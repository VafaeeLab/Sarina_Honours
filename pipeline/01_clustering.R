# =============================================================================
# 01_clustering.R
# Outcome-independent clustering audit + frozen primary dietary profiles
#
# PURPOSE
#   1. Read the frozen diet-preprocessing RDS from 00_preprocess.R.
#   2. Reproduce the thesis clustering decision audit across Ward.D2, K-means,
#      GMM, fuzzy c-means (FCM), and UMAP + K-means.
#   3. Tune each method WITHOUT using HCC/case status.
#   4. Compare frozen candidates using repeated-subsample Jaccard stability.
#   5. Freeze the configured primary clustering solution and save participant
#      assignments for downstream Stage 3/4 analyses.
#
# INPUT
#   RDS_PATHS$preprocess
#
# OUTPUTS
#   RDS_PATHS$clustering
#   intermediate/01_clustering_audit.rds
#   results/01_clustering/*.csv, *.png, *.pdf
#   figures/01_clustering/primary_profile_heatmap.png/.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(psych)
  library(cluster)
  library(mclust)
  library(e1071)
  library(uwot)
  library(gridExtra)
  library(grid)
  library(scales)
})

source("pipeline/config.R")
source("pipeline/utils.R")

set.seed(MASTER_SEED)

# =============================================================================
# 0. CONFIGURATION ALIASES + OUTPUT LOCATIONS
# =============================================================================
K_MAX               <- as.integer(CLUSTER_K_MAX)
PC_MAX_CAP          <- as.integer(CLUSTER_PC_MAX_CAP)
PA_ITER             <- as.integer(CLUSTER_PA_ITER)
GAP_B               <- as.integer(CLUSTER_GAP_B)
TUNE_RESAMPLES      <- as.integer(CLUSTER_TUNE_RESAMPLES)
FINAL_RESAMPLES     <- as.integer(CLUSTER_FINAL_RESAMPLES)
SUBSAMPLE_FRACTION  <- as.numeric(CLUSTER_SUBSAMPLE_FRACTION)
FCM_M_VALUES        <- as.numeric(CLUSTER_FCM_M_VALUES)
UMAP_NEIGHBORS      <- as.integer(CLUSTER_UMAP_NEIGHBORS)
UMAP_MIN_DIST       <- as.numeric(CLUSTER_UMAP_MIN_DIST)

OUTPUT_DIR <- file.path(DIRS$results, "01_clustering")
FIGURE_DIR <- file.path(DIRS$figures, "01_clustering")
AUDIT_RDS  <- file.path(DIRS$intermediate, "01_clustering_audit.rds")

ensure_dir(OUTPUT_DIR)
ensure_dir(FIGURE_DIR)
ensure_parent_dir(AUDIT_RDS)
ensure_parent_dir(RDS_PATHS$clustering)

if (!CLUSTER_PRIMARY_SELECTION %in% c("configured", "stability")) {
  stop(
    "CLUSTER_PRIMARY_SELECTION must be 'configured' or 'stability'.",
    call. = FALSE
  )
}

if (K_MAX < 2L) stop("CLUSTER_K_MAX must be >= 2.", call. = FALSE)
if (PC_MAX_CAP < 2L) stop("CLUSTER_PC_MAX_CAP must be >= 2.", call. = FALSE)
if (PA_ITER < 20L) stop("CLUSTER_PA_ITER must be >= 20.", call. = FALSE)
if (GAP_B < 20L) stop("CLUSTER_GAP_B must be >= 20.", call. = FALSE)
if (TUNE_RESAMPLES < 10L) stop("CLUSTER_TUNE_RESAMPLES must be >= 10.", call. = FALSE)
if (FINAL_RESAMPLES < 20L) stop("CLUSTER_FINAL_RESAMPLES must be >= 20.", call. = FALSE)
if (!is.finite(SUBSAMPLE_FRACTION) || SUBSAMPLE_FRACTION <= 0.5 || SUBSAMPLE_FRACTION >= 1) {
  stop("CLUSTER_SUBSAMPLE_FRACTION must be between 0.5 and 1.", call. = FALSE)
}

# =============================================================================
# 1. READ PREPROCESSING OUTPUT
# =============================================================================

preprocess <- read_rds_checked(
  RDS_PATHS$preprocess,
  "preprocessed diet RDS"
)

if (is.null(preprocess$adj_scaled) || is.null(preprocess$diet_final)) {
  stop(
    "00_preprocessed_diet.rds must contain adj_scaled and diet_final.",
    call. = FALSE
  )
}

adj_scaled <- preprocess$adj_scaled
diet_final <- preprocess$diet_final

X <- as.matrix(adj_scaled)
storage.mode(X) <- "double"

if (!is.numeric(X) || any(!is.finite(X))) {
  stop("adj_scaled must be a finite numeric matrix with no NA/Inf values.", call. = FALSE)
}

if (nrow(X) != nrow(diet_final)) {
  stop("nrow(adj_scaled) must equal nrow(diet_final).", call. = FALSE)
}

if (ncol(X) < 2L) {
  stop("At least two food-group variables are required for clustering.", call. = FALSE)
}

if (nrow(X) < 30L) {
  stop("Too few participants for this clustering audit workflow.", call. = FALSE)
}

assert_columns(diet_final, "case_status", "preprocessed diet_final")

participant_id <- if ("id" %in% names(diet_final)) {
  normalise_id(diet_final$id)
} else if ("study_id" %in% names(diet_final)) {
  normalise_id(diet_final$study_id)
} else if (!is.null(rownames(X))) {
  normalise_id(rownames(X))
} else {
  as.character(seq_len(nrow(X)))
}

if (anyNA(participant_id) || anyDuplicated(participant_id)) {
  stop("Clustering participant IDs must be unique and non-missing.", call. = FALSE)
}

rownames(X) <- participant_id

n_obs <- nrow(X)
n_fg  <- ncol(X)

cat("\n============================================================\n")
cat("01 CLUSTERING — OUTCOME-INDEPENDENT DECISION AUDIT\n")
cat("============================================================\n")
cat("Participants:", n_obs, "\n")
cat("Food-group variables:", n_fg, "\n")
cat("Primary-selection mode:", CLUSTER_PRIMARY_SELECTION, "\n")
if (CLUSTER_PRIMARY_SELECTION == "configured") {
  cat("Configured primary method:", PRIMARY_CLUSTER_METHOD, "\n")
}
cat("Disease status is NOT used in tuning/selection.\n")

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

safe_silhouette <- function(labels, d) {
  labels <- as.integer(factor(labels))
  if (length(unique(labels)) < 2L) return(NA_real_)
  if (any(table(labels) < 2L)) return(NA_real_)
  tryCatch(
    mean(cluster::silhouette(labels, d)[, 3]),
    error = function(e) NA_real_
  )
}

run_gap <- function(data, FUNcluster, label, K.max = K_MAX,
                    B = GAP_B, seed = MASTER_SEED) {
  cat("  Gap:", label, "... ")
  set.seed(seed)
  g <- cluster::clusGap(
    x = as.matrix(data),
    FUNcluster = FUNcluster,
    K.max = K.max,
    B = B,
    d.power = 2,
    spaceH0 = "scaledPCA",
    verbose = FALSE
  )
  best_k <- cluster::maxSE(
    f = g$Tab[, "gap"],
    SE.f = g$Tab[, "SE.sim"],
    method = "Tibs2001SEmax"
  )
  details <- tibble(
    k = seq_len(K.max),
    gap = as.numeric(g$Tab[, "gap"]),
    se = as.numeric(g$Tab[, "SE.sim"]),
    selected = seq_len(K.max) == best_k
  )
  cat("selected k =", best_k, "\n")
  list(object = g, best_k = as.integer(best_k), details = details)
}

best_match_jaccard <- function(reference_labels, candidate_labels) {
  reference_labels <- as.character(reference_labels)
  candidate_labels <- as.character(candidate_labels)
  
  if (length(reference_labels) != length(candidate_labels)) {
    stop("Jaccard label vectors must have equal lengths.")
  }
  if (anyNA(reference_labels) || anyNA(candidate_labels)) {
    return(rep(NA_real_, length(unique(reference_labels))))
  }
  
  ref_levels <- unique(reference_labels)
  cand_levels <- unique(candidate_labels)
  
  vapply(ref_levels, function(r) {
    ref_set <- which(reference_labels == r)
    if (length(ref_set) == 0L) return(NA_real_)
    
    scores <- vapply(cand_levels, function(cn) {
      cand_set <- which(candidate_labels == cn)
      union_n <- length(union(ref_set, cand_set))
      if (union_n == 0L) return(NA_real_)
      length(intersect(ref_set, cand_set)) / union_n
    }, numeric(1))
    
    if (all(is.na(scores))) NA_real_ else max(scores, na.rm = TRUE)
  }, numeric(1))
}

resample_stability <- function(X, full_labels, fit_fun,
                               B = FINAL_RESAMPLES,
                               fraction = SUBSAMPLE_FRACTION,
                               seed = MASTER_SEED) {
  if (length(unique(full_labels)) < 2L) {
    return(list(
      mean_jaccard = NA_real_,
      min_cluster_jaccard = NA_real_,
      cluster_jaccard = NA_real_,
      success_rate = NA_real_,
      valid_runs = 0L,
      B = B
    ))
  }
  
  n <- nrow(X)
  size <- max(length(unique(full_labels)) * 3L,
              floor(fraction * n))
  size <- min(size, n - 1L)
  
  ref_levels <- unique(as.character(full_labels))
  jaccard_matrix <- matrix(
    NA_real_, nrow = B, ncol = length(ref_levels),
    dimnames = list(NULL, ref_levels)
  )
  success <- logical(B)
  
  for (b in seq_len(B)) {
    set.seed(seed + b)
    idx <- sort(sample.int(n, size = size, replace = FALSE))
    
    fitted <- tryCatch(
      fit_fun(X[idx, , drop = FALSE], seed = seed + 10000L + b),
      error = function(e) NULL
    )
    
    if (is.null(fitted) || is.null(fitted$labels) ||
        length(fitted$labels) != length(idx) || anyNA(fitted$labels) ||
        length(unique(fitted$labels)) < 2L) {
      next
    }
    
    j <- best_match_jaccard(full_labels[idx], fitted$labels)
    names(j) <- unique(as.character(full_labels[idx]))
    shared_names <- intersect(colnames(jaccard_matrix), names(j))
    jaccard_matrix[b, shared_names] <- j[shared_names]
    success[b] <- TRUE
  }
  
  cluster_means <- colMeans(jaccard_matrix[success, , drop = FALSE],
                            na.rm = TRUE)
  cluster_means[!is.finite(cluster_means)] <- NA_real_
  
  list(
    mean_jaccard = if (all(is.na(cluster_means))) NA_real_ else
      mean(cluster_means, na.rm = TRUE),
    min_cluster_jaccard = if (all(is.na(cluster_means))) NA_real_ else
      min(cluster_means, na.rm = TRUE),
    cluster_jaccard = cluster_means,
    success_rate = mean(success),
    valid_runs = sum(success),
    B = B,
    replicate_jaccard = jaccard_matrix
  )
}

fisher_after_freeze <- function(labels, status, seed = MASTER_SEED) {
  if (length(unique(labels)) < 2L) return(NA_real_)
  tab <- table(factor(labels), factor(status))
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(NA_real_)
  
  exact_result <- tryCatch(
    stats::fisher.test(tab),
    error = function(e) NULL
  )
  if (!is.null(exact_result)) return(exact_result$p.value)
  
  set.seed(seed)
  stats::fisher.test(tab, simulate.p.value = TRUE, B = 100000)$p.value
}
# Calculate the fuzzy c-means objective for a fitted model:
# sum_i sum_j u_ij^m ||x_i - v_j||^2
fcm_objective_value <- function(data, fit, m) {
  X_local <- as.matrix(data)
  U <- as.matrix(fit$membership)
  V <- as.matrix(fit$centers)
  
  if (
    nrow(U) != nrow(X_local) ||
    ncol(U) != nrow(V) ||
    anyNA(U) || anyNA(V) ||
    any(!is.finite(U)) || any(!is.finite(V))
  ) {
    return(Inf)
  }
  
  squared_distances <- vapply(
    seq_len(nrow(V)),
    function(j) {
      rowSums(
        (
          X_local - matrix(
            V[j, ],
            nrow = nrow(X_local),
            ncol = ncol(X_local),
            byrow = TRUE
          )
        )^2
      )
    },
    numeric(nrow(X_local))
  )
  
  objective <- sum((U^m) * as.matrix(squared_distances))
  if (!is.finite(objective)) Inf else objective
}

# Xie-Beni validity index:
#   numerator   = fuzzy within-cluster objective
#   denominator = n times the minimum squared distance between centres
# Lower values indicate more compact clusters with better separated centres.
xie_beni_index <- function(data, fit, m, tolerance = 1e-10) {
  X_local <- as.matrix(data)
  U <- as.matrix(fit$membership)
  V <- as.matrix(fit$centers)
  
  if (
    nrow(V) < 2L ||
    nrow(U) != nrow(X_local) ||
    ncol(U) != nrow(V) ||
    anyNA(U) || anyNA(V) ||
    any(!is.finite(U)) || any(!is.finite(V))
  ) {
    return(Inf)
  }
  
  numerator <- fcm_objective_value(X_local, fit, m)
  
  centre_distances_squared <- as.matrix(stats::dist(V))^2
  diag(centre_distances_squared) <- Inf
  minimum_centre_distance <- min(centre_distances_squared, na.rm = TRUE)
  
  if (
    !is.finite(numerator) ||
    !is.finite(minimum_centre_distance) ||
    minimum_centre_distance <= tolerance
  ) {
    return(Inf)
  }
  
  numerator / (nrow(X_local) * minimum_centre_distance)
}

# Fit a fixed k x m FCM configuration from several random starts.
# The lowest-objective valid fit is retained before Xie-Beni is calculated.
fit_best_fcm <- function(
    data,
    k,
    m,
    n_starts = FCM_N_STARTS,
    seed = MASTER_SEED,
    degeneracy_tolerance = 1e-4
) {
  X_local <- as.matrix(data)
  
  best_fit <- NULL
  best_objective <- Inf
  successful_starts <- 0L
  valid_starts <- 0L
  
  for (start_id in seq_len(n_starts)) {
    set.seed(seed + start_id)
    
    fit <- tryCatch(
      e1071::cmeans(
        x = X_local,
        centers = k,
        m = m,
        iter.max = 500,
        dist = "euclidean",
        method = "cmeans"
      ),
      error = function(e) NULL
    )
    
    if (is.null(fit)) next
    successful_starts <- successful_starts + 1L
    
    U <- as.matrix(fit$membership)
    hard_labels <- as.integer(fit$cluster)
    cluster_sizes <- tabulate(hard_labels, nbins = k)
    
    if (
      nrow(U) != nrow(X_local) ||
      ncol(U) != k ||
      anyNA(U) || any(!is.finite(U)) ||
      length(unique(hard_labels)) != k ||
      any(cluster_sizes == 0L)
    ) {
      next
    }
    
    partition_coefficient <- sum(U^2) / nrow(U)
    mean_max_membership <- mean(apply(U, 1, max))
    chance_membership <- 1 / k
    
    if (
      partition_coefficient <= chance_membership + degeneracy_tolerance ||
      mean_max_membership <= chance_membership + degeneracy_tolerance
    ) {
      next
    }
    
    xb <- xie_beni_index(X_local, fit, m)
    if (!is.finite(xb)) next
    
    valid_starts <- valid_starts + 1L
    objective <- fcm_objective_value(X_local, fit, m)
    
    if (is.finite(objective) && objective < best_objective) {
      best_fit <- fit
      best_objective <- objective
    }
  }
  
  if (is.null(best_fit)) return(NULL)
  
  list(
    fit = best_fit,
    objective = best_objective,
    xie_beni = xie_beni_index(X_local, best_fit, m),
    successful_starts = successful_starts,
    valid_starts = valid_starts
  )
}

extract_best_icl <- function(data, G = 1:K_MAX) {
  icl_obj <- tryCatch(
    mclust::mclustICL(data, G = G),
    error = function(e) NULL
  )
  if (is.null(icl_obj)) return(NULL)
  
  icl_mat <- as.matrix(icl_obj)
  if (all(!is.finite(icl_mat))) return(NULL)
  
  idx <- which(icl_mat == max(icl_mat, na.rm = TRUE), arr.ind = TRUE)[1, ]
  best_G <- suppressWarnings(as.integer(rownames(icl_mat)[idx[1]]))
  best_model <- colnames(icl_mat)[idx[2]]
  best_value <- as.numeric(icl_mat[idx[1], idx[2]])
  
  if (!is.finite(best_G)) return(NULL)
  
  fit <- tryCatch({
    if (best_G == 1L) {
      mclust::Mclust(data, G = 1L, verbose = FALSE)
    } else {
      mclust::Mclust(data, G = best_G,
                     modelNames = best_model, verbose = FALSE)
    }
  }, error = function(e) NULL)
  
  if (is.null(fit)) return(NULL)
  
  list(
    G = best_G,
    model = best_model,
    ICL = best_value,
    fit = fit,
    labels = fit$classification,
    uncertainty = if (best_G > 1L) mean(fit$uncertainty) else 0
  )
}

format_value <- function(x, digits = 3L) {
  if (length(x) == 0L || is.na(x) || !is.finite(x)) return("—")
  formatC(x, format = "f", digits = digits)
}

safe_filename <- function(x) {
  gsub("[^A-Za-z0-9]+", "_", x)
}

# =============================================================================
# PHASE 1 — PCA RANGE FROM PARALLEL ANALYSIS
# =============================================================================

cat("\n--- PHASE 1: PCA CANDIDATE RANGE ---\n")

set.seed(MASTER_SEED)
pa <- psych::fa.parallel(
  X,
  fa = "pc",
  n.iter = PA_ITER,
  plot = FALSE,
  show.legend = FALSE,
  error.bars = FALSE
)

pa_central <- suppressWarnings(as.integer(pa$ncomp))
if (!is.finite(pa_central) || pa_central < 1L) pa_central <- 1L
pa_central <- min(pa_central, n_fg)

q_max <- min(n_fg, PC_MAX_CAP, max(2L, pa_central + 2L))
pc_range <- seq.int(2L, q_max)

pca_full <- prcomp(X, center = FALSE, scale. = FALSE)

# =============================================================================
# NUMBER OF PCs REQUIRED TO EXPLAIN AT LEAST 80% OF VARIANCE
# =============================================================================

VARIANCE_THRESHOLD <- 80

# Variance explained by each individual PC
pc_variance_percent <- 100 * (
  pca_full$sdev^2 /
    sum(pca_full$sdev^2)
)

# Cumulative variance explained
pc_cumulative_variance_percent <- cumsum(
  pc_variance_percent
)

# Minimum number of PCs required to reach at least 80%
n_pcs_for_80 <- which(
  pc_cumulative_variance_percent >= VARIANCE_THRESHOLD
)[1]

# Full PCA variance table
pc_variance_full <- tibble::tibble(
  PC = seq_along(pc_variance_percent),
  
  PC_label = paste0(
    "PC",
    seq_along(pc_variance_percent)
  ),
  
  individual_variance_percent = round(
    pc_variance_percent,
    2
  ),
  
  cumulative_variance_percent = round(
    pc_cumulative_variance_percent,
    2
  ),
  
  reaches_80_percent = seq_along(
    pc_variance_percent
  ) == n_pcs_for_80
)

cat("\n============================================================\n")
cat("PCA VARIANCE SUMMARY\n")
cat("============================================================\n")

print(
  pc_variance_full,
  n = Inf,
  width = Inf
)

cat(
  "\nMinimum PCs required for at least",
  VARIANCE_THRESHOLD,
  "% cumulative variance:",
  n_pcs_for_80,
  "PCs\n"
)

cat(
  "Cumulative variance at PC",
  n_pcs_for_80,
  ":",
  round(
    pc_cumulative_variance_percent[n_pcs_for_80],
    2
  ),
  "%\n"
)

if (n_pcs_for_80 > 1L) {
  cat(
    "Cumulative variance at PC",
    n_pcs_for_80 - 1L,
    ":",
    round(
      pc_cumulative_variance_percent[n_pcs_for_80 - 1L],
      2
    ),
    "%\n"
  )
}

# Save all PC variance results
readr::write_csv(
  pc_variance_full,
  file.path(
    OUTPUT_DIR,
    "pca_variance_all_components.csv"
  )
)

variance <- summary(pca_full)$importance
pc1_pct <- 100 * variance[2, 1]
pc2_pct <- 100 * variance[2, 2]

cat("Parallel analysis retained:", pa_central, "PC(s)\n")
cat("Candidate PCA range:", paste(pc_range, collapse = ", "), "\n")

pc_variance_table <- tibble(
  n_pcs = pc_range,
  cumulative_variance_pct = 100 * variance[3, pc_range]
)
readr::write_csv(pc_variance_table,
                 file.path(OUTPUT_DIR, "pc_variance_candidates.csv"))

# =============================================================================
# PHASE 2 — WARD.D2: k SELECTED BY GAP INCLUDING k=1
# =============================================================================

cat("\n--- PHASE 2: WARD.D2 ---\n")

gap_ward <- run_gap(
  X,
  FUNcluster = function(x, k) {
    hc <- hclust(dist(x, method = "euclidean"), method = "ward.D2")
    list(cluster = cutree(hc, k = k))
  },
  label = "Ward.D2 on full matrix"
)

ward_k <- gap_ward$best_k
ward_gap_value <- gap_ward$details$gap[gap_ward$details$k == ward_k]
ward_labels <- if (ward_k == 1L) {
  rep(1L, n_obs)
} else {
  cutree(hclust(dist(X), method = "ward.D2"), k = ward_k)
}

ward_fit_fun <- function(data, seed) {
  hc <- hclust(dist(data), method = "ward.D2")
  list(labels = cutree(hc, k = ward_k))
}

readr::write_csv(gap_ward$details,
                 file.path(OUTPUT_DIR, "ward_gap_details.csv"))

# =============================================================================
# PHASE 3 — K-MEANS: JOINT PC RANGE × k VIA GAP, PC COUNT VIA STABILITY
# =============================================================================

cat("\n--- PHASE 3: K-MEANS PCA GRID ---\n")

kmeans_candidates <- vector("list", length(pc_range))

for (ii in seq_along(pc_range)) {
  q <- pc_range[ii]
  pc_data <- pca_full$x[, seq_len(q), drop = FALSE]
  
  gap_q <- run_gap(
    pc_data,
    FUNcluster = function(x, k) {
      set.seed(MASTER_SEED + q + k)
      list(cluster = kmeans(x, centers = k, nstart = 50,
                            iter.max = 300)$cluster)
    },
    label = paste0("K-means on PC1-", q),
    seed = MASTER_SEED + q
  )
  
  k_q <- gap_q$best_k
  full_labels <- if (k_q == 1L) {
    rep(1L, n_obs)
  } else {
    set.seed(MASTER_SEED + q)
    kmeans(pc_data, centers = k_q, nstart = 100,
           iter.max = 300)$cluster
  }
  
  tune_stability <- if (k_q == 1L) {
    list(mean_jaccard = NA_real_, min_cluster_jaccard = NA_real_,
         success_rate = NA_real_)
  } else {
    fit_fun_q <- local({
      q_local <- q
      k_local <- k_q
      function(data, seed) {
        pc <- prcomp(data, center = FALSE, scale. = FALSE)$x[
          , seq_len(q_local), drop = FALSE
        ]
        set.seed(seed)
        km <- kmeans(pc, centers = k_local, nstart = 50, iter.max = 300)
        list(labels = km$cluster)
      }
    })
    resample_stability(
      X = X,
      full_labels = full_labels,
      fit_fun = fit_fun_q,
      B = TUNE_RESAMPLES,
      fraction = SUBSAMPLE_FRACTION,
      seed = MASTER_SEED + 100L * q
    )
  }
  
  kmeans_candidates[[ii]] <- tibble(
    n_pcs = q,
    best_k = k_q,
    gap_value = gap_q$details$gap[gap_q$details$k == k_q],
    tune_mean_jaccard = tune_stability$mean_jaccard,
    tune_min_jaccard = tune_stability$min_cluster_jaccard,
    tune_success_rate = tune_stability$success_rate,
    labels = list(full_labels),
    gap_details = list(gap_q$details)
  )
}

kmeans_candidates <- bind_rows(kmeans_candidates)

kmeans_multi <- kmeans_candidates %>% filter(best_k > 1L)
if (nrow(kmeans_multi) == 0L) {
  selected_km <- kmeans_candidates %>%
    mutate(distance_from_pa = abs(n_pcs - pa_central)) %>%
    arrange(distance_from_pa, n_pcs) %>%
    slice(1)
} else {
  selected_km <- kmeans_multi %>%
    arrange(desc(tune_min_jaccard), desc(tune_mean_jaccard),
            best_k, n_pcs) %>%
    slice(1)
}

km_q <- selected_km$n_pcs
km_k <- selected_km$best_k
km_labels <- selected_km$labels[[1]]
km_gap_value <- selected_km$gap_value

km_fit_fun <- local({
  q_local <- km_q
  k_local <- km_k
  function(data, seed) {
    pc <- prcomp(data, center = FALSE, scale. = FALSE)$x[
      , seq_len(q_local), drop = FALSE
    ]
    if (k_local == 1L) return(list(labels = rep(1L, nrow(data))))
    set.seed(seed)
    fit <- kmeans(pc, centers = k_local, nstart = 100, iter.max = 300)
    list(labels = fit$cluster)
  }
})

readr::write_csv(
  kmeans_candidates %>% select(-labels, -gap_details),
  file.path(OUTPUT_DIR, "kmeans_pc_candidates.csv")
)

# =============================================================================
# PHASE 4 — GMM: ICL SELECTS G WITHIN PC COUNT; PC COUNT VIA STABILITY
# =============================================================================

cat("\n--- PHASE 4: GMM PCA GRID ---\n")

gmm_candidates <- vector("list", length(pc_range))

for (ii in seq_along(pc_range)) {
  q <- pc_range[ii]
  pc_data <- pca_full$x[, seq_len(q), drop = FALSE]
  best <- extract_best_icl(pc_data, G = 1:K_MAX)
  
  if (is.null(best)) {
    gmm_candidates[[ii]] <- tibble(
      n_pcs = q,
      best_k = NA_integer_,
      model = "FAILED",
      ICL = NA_real_,
      mean_uncertainty = NA_real_,
      tune_mean_jaccard = NA_real_,
      tune_min_jaccard = NA_real_,
      tune_success_rate = 0,
      labels = list(rep(NA_integer_, n_obs))
    )
    next
  }
  
  tune_stability <- if (best$G == 1L) {
    list(mean_jaccard = NA_real_, min_cluster_jaccard = NA_real_,
         success_rate = NA_real_)
  } else {
    fit_fun_q <- local({
      q_local <- q
      k_local <- best$G
      model_local <- best$model
      function(data, seed) {
        pc <- prcomp(data, center = FALSE, scale. = FALSE)$x[
          , seq_len(q_local), drop = FALSE
        ]
        fit <- tryCatch(
          mclust::Mclust(pc, G = k_local,
                         modelNames = model_local, verbose = FALSE),
          error = function(e) NULL
        )
        if (is.null(fit)) return(NULL)
        list(labels = fit$classification)
      }
    })
    resample_stability(
      X = X,
      full_labels = best$labels,
      fit_fun = fit_fun_q,
      B = TUNE_RESAMPLES,
      fraction = SUBSAMPLE_FRACTION,
      seed = MASTER_SEED + 200L * q
    )
  }
  
  gmm_candidates[[ii]] <- tibble(
    n_pcs = q,
    best_k = best$G,
    model = best$model,
    ICL = best$ICL,
    mean_uncertainty = best$uncertainty,
    tune_mean_jaccard = tune_stability$mean_jaccard,
    tune_min_jaccard = tune_stability$min_cluster_jaccard,
    tune_success_rate = tune_stability$success_rate,
    labels = list(best$labels)
  )
}

gmm_candidates <- bind_rows(gmm_candidates)
gmm_valid <- gmm_candidates %>% filter(!is.na(best_k))
if (nrow(gmm_valid) == 0L) stop("All GMM candidate fits failed.")

gmm_multi <- gmm_valid %>% filter(best_k > 1L)
if (nrow(gmm_multi) == 0L) {
  selected_gmm <- gmm_valid %>%
    mutate(distance_from_pa = abs(n_pcs - pa_central)) %>%
    arrange(distance_from_pa, n_pcs) %>%
    slice(1)
} else {
  selected_gmm <- gmm_multi %>%
    arrange(desc(tune_min_jaccard), desc(tune_mean_jaccard),
            best_k, n_pcs) %>%
    slice(1)
}

gmm_q <- selected_gmm$n_pcs
gmm_k <- selected_gmm$best_k
gmm_model <- selected_gmm$model
gmm_labels <- selected_gmm$labels[[1]]
gmm_icl_value <- selected_gmm$ICL

gmm_fit_fun <- local({
  q_local <- gmm_q
  k_local <- gmm_k
  model_local <- gmm_model
  function(data, seed) {
    pc <- prcomp(data, center = FALSE, scale. = FALSE)$x[
      , seq_len(q_local), drop = FALSE
    ]
    if (k_local == 1L) return(list(labels = rep(1L, nrow(data))))
    fit <- tryCatch(
      mclust::Mclust(pc, G = k_local,
                     modelNames = model_local, verbose = FALSE),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    list(labels = fit$classification)
  }
})

readr::write_csv(
  gmm_candidates %>% select(-labels),
  file.path(OUTPUT_DIR, "gmm_pc_candidates.csv")
)

# =============================================================================
# PHASE 5 — FCM: JOINT k × m, SELECTED BY MINIMUM XIE-BENI
# =============================================================================

cat("\n--- PHASE 5: FCM GRID ---\n")

fcm_grid <- tidyr::crossing(
  k = 2:K_MAX,
  m = FCM_M_VALUES
)

fcm_results <- vector("list", nrow(fcm_grid))

for (ii in seq_len(nrow(fcm_grid))) {
  k <- fcm_grid$k[ii]
  m <- fcm_grid$m[ii]
  
  cat("  FCM: k =", k, "| m =", m, "... ")
  
  fitted <- fit_best_fcm(
    data = X,
    k = k,
    m = m,
    n_starts = FCM_N_STARTS,
    seed = MASTER_SEED + ii * 10000L
  )
  
  if (is.null(fitted)) {
    cat("no valid fit\n")
    
    fcm_results[[ii]] <- tibble(
      k = k,
      m = m,
      xie_beni = NA_real_,
      objective = NA_real_,
      partition_coefficient = NA_real_,
      partition_entropy = NA_real_,
      mean_max_membership = NA_real_,
      minimum_cluster_size = 0L,
      successful_starts = 0L,
      valid_starts = 0L,
      labels = list(rep(NA_integer_, n_obs)),
      membership = list(NULL),
      fit = list(NULL)
    )
    next
  }
  
  fit <- fitted$fit
  U <- as.matrix(fit$membership)
  hard_labels <- as.integer(fit$cluster)
  cluster_sizes <- tabulate(hard_labels, nbins = k)
  
  fcm_results[[ii]] <- tibble(
    k = k,
    m = m,
    xie_beni = fitted$xie_beni,
    objective = fitted$objective,
    partition_coefficient = sum(U^2) / nrow(U),
    partition_entropy = -sum(U * log(U + 1e-12)) / nrow(U),
    mean_max_membership = mean(apply(U, 1, max)),
    minimum_cluster_size = min(cluster_sizes),
    successful_starts = fitted$successful_starts,
    valid_starts = fitted$valid_starts,
    labels = list(hard_labels),
    membership = list(U),
    fit = list(fit)
  )
  
  cat("XB =", round(fitted$xie_beni, 4), "\n")
}

fcm_results <- bind_rows(fcm_results)

# Lower Xie-Beni is better. Fewer clusters and m closest to 2 are used only
# as deterministic tie-breakers when criterion values are equal.
selected_fcm <- fcm_results %>%
  filter(is.finite(xie_beni)) %>%
  arrange(xie_beni, k, abs(m - 2)) %>%
  slice(1)

if (nrow(selected_fcm) == 0L) {
  stop(
    "All FCM candidates failed or collapsed. Inspect fcm_grid.csv and the input matrix."
  )
}

fcm_k <- selected_fcm$k
fcm_m <- selected_fcm$m
fcm_labels <- selected_fcm$labels[[1]]
fcm_membership <- selected_fcm$membership[[1]]
fcm_selected_fit <- selected_fcm$fit[[1]]
fcm_xb_value <- selected_fcm$xie_beni

cat(
  "Selected FCM configuration:",
  "k =", fcm_k,
  "| m =", fcm_m,
  "| Xie-Beni =", round(fcm_xb_value, 4),
  "| valid starts =", selected_fcm$valid_starts,
  "of", FCM_N_STARTS,
  "\n"
)

fcm_fit_fun <- local({
  k_local <- fcm_k
  m_local <- fcm_m
  n_starts_local <- FCM_N_STARTS
  
  function(data, seed) {
    fitted <- fit_best_fcm(
      data = data,
      k = k_local,
      m = m_local,
      n_starts = n_starts_local,
      seed = seed
    )
    
    if (is.null(fitted)) return(NULL)
    list(labels = fitted$fit$cluster)
  }
})

readr::write_csv(
  fcm_results %>% select(-labels, -membership, -fit),
  file.path(OUTPUT_DIR, "fcm_grid.csv")
)

# =============================================================================
# PHASE 6 — UMAP+K-MEANS: k BY NATIVE SILHOUETTE; PARAMETERS BY STABILITY
# =============================================================================

cat("\n--- PHASE 6: UMAP PARAMETER GRID ---\n")

umap_grid <- tidyr::crossing(
  n_neighbors = UMAP_NEIGHBORS,
  min_dist = UMAP_MIN_DIST
)
umap_results <- vector("list", nrow(umap_grid))

for (ii in seq_len(nrow(umap_grid))) {
  nn <- umap_grid$n_neighbors[ii]
  md <- umap_grid$min_dist[ii]
  
  set.seed(MASTER_SEED + ii)
  emb <- uwot::umap(
    X,
    n_neighbors = nn,
    min_dist = md,
    n_components = 2,
    metric = "euclidean",
    seed = MASTER_SEED + ii,
    verbose = FALSE
  )
  
  d_emb <- dist(emb)
  k_stats <- map_dfr(2:K_MAX, function(k) {
    set.seed(MASTER_SEED + 100L * ii + k)
    fit <- kmeans(emb, centers = k, nstart = 100, iter.max = 300)
    tibble(
      k = k,
      silhouette = safe_silhouette(fit$cluster, d_emb),
      labels = list(fit$cluster)
    )
  })
  
  best_k_row <- k_stats %>%
    arrange(desc(silhouette), k) %>%
    slice(1)
  
  k_selected <- best_k_row$k
  labels_selected <- best_k_row$labels[[1]]
  
  fit_fun_params <- local({
    nn_local <- nn
    md_local <- md
    k_local <- k_selected
    function(data, seed) {
      emb_sub <- uwot::umap(
        data,
        n_neighbors = min(nn_local, nrow(data) - 1L),
        min_dist = md_local,
        n_components = 2,
        metric = "euclidean",
        seed = seed,
        verbose = FALSE
      )
      set.seed(seed + 1L)
      km <- kmeans(emb_sub, centers = k_local,
                   nstart = 50, iter.max = 300)
      list(labels = km$cluster)
    }
  })
  
  tune_stability <- resample_stability(
    X = X,
    full_labels = labels_selected,
    fit_fun = fit_fun_params,
    B = TUNE_RESAMPLES,
    fraction = SUBSAMPLE_FRACTION,
    seed = MASTER_SEED + 1000L * ii
  )
  
  umap_results[[ii]] <- tibble(
    n_neighbors = nn,
    min_dist = md,
    best_k = k_selected,
    silhouette = best_k_row$silhouette,
    tune_mean_jaccard = tune_stability$mean_jaccard,
    tune_min_jaccard = tune_stability$min_cluster_jaccard,
    tune_success_rate = tune_stability$success_rate,
    labels = list(labels_selected)
  )
}

umap_results <- bind_rows(umap_results)
selected_umap <- umap_results %>%
  arrange(desc(tune_min_jaccard), desc(tune_mean_jaccard),
          desc(silhouette), best_k,
          abs(n_neighbors - 15), abs(min_dist - 0.15)) %>%
  slice(1)

umap_nn <- selected_umap$n_neighbors
umap_md <- selected_umap$min_dist
umap_k <- selected_umap$best_k
umap_labels <- selected_umap$labels[[1]]
umap_sil_value <- selected_umap$silhouette

umap_fit_fun <- local({
  nn_local <- umap_nn
  md_local <- umap_md
  k_local <- umap_k
  function(data, seed) {
    emb <- uwot::umap(
      data,
      n_neighbors = min(nn_local, nrow(data) - 1L),
      min_dist = md_local,
      n_components = 2,
      metric = "euclidean",
      seed = seed,
      verbose = FALSE
    )
    set.seed(seed + 1L)
    fit <- kmeans(emb, centers = k_local, nstart = 100, iter.max = 300)
    list(labels = fit$cluster)
  }
})

readr::write_csv(
  umap_results %>% select(-labels),
  file.path(OUTPUT_DIR, "umap_grid.csv")
)

# =============================================================================
# PHASE 7 — FINAL STABILITY OF EACH FROZEN CANDIDATE
# =============================================================================

cat("\n--- PHASE 7: FINAL SUBSAMPLE STABILITY ---\n")
cat("B =", FINAL_RESAMPLES,
    "|", SUBSAMPLE_FRACTION * 100, "% subsamples without replacement\n")
cat("PCA/UMAP and clustering are refitted in every subsample.\n")

candidate_objects <- list(
  "Ward.D2" = list(
    method = "Ward.D2",
    space = "Full 18D",
    best_k = ward_k,
    criterion = "Gap (1-SE)",
    criterion_value = ward_gap_value,
    labels = ward_labels,
    fit_fun = ward_fit_fun
  ),
  "K-means" = list(
    method = "K-means",
    space = paste0("PCA ", km_q, "D"),
    best_k = km_k,
    criterion = "Gap; PCs→stability",
    criterion_value = km_gap_value,
    labels = km_labels,
    fit_fun = km_fit_fun
  ),
  "GMM" = list(
    method = "GMM",
    space = paste0("PCA ", gmm_q, "D"),
    best_k = gmm_k,
    criterion = "ICL; PCs→stability",
    criterion_value = gmm_icl_value,
    labels = gmm_labels,
    fit_fun = gmm_fit_fun,
    model = gmm_model
  ),
  "FCM" = list(
    method = "FCM",
    space = paste0("Full 18D (m=", fcm_m, ")"),
    best_k = fcm_k,
    criterion = "Xie-Beni (lower; k>=2)",
    criterion_value = fcm_xb_value,
    labels = fcm_labels,
    membership = fcm_membership,
    fit_fun = fcm_fit_fun,
    m = fcm_m,
    n_starts = FCM_N_STARTS
  ),
  "UMAP+Km" = list(
    method = "UMAP+Km",
    space = paste0("UMAP 2D (", umap_nn, "/", umap_md, ")"),
    best_k = umap_k,
    criterion = "Silh.; params→stability (k≥2)",
    criterion_value = umap_sil_value,
    labels = umap_labels,
    fit_fun = umap_fit_fun,
    n_neighbors = umap_nn,
    min_dist = umap_md
  )
)

for (method_name in names(candidate_objects)) {
  cat("  ", method_name, "... ", sep = "")
  obj <- candidate_objects[[method_name]]
  
  stability <- if (obj$best_k == 1L) {
    resample_stability(X, obj$labels, obj$fit_fun,
                       B = FINAL_RESAMPLES)
  } else {
    resample_stability(
      X = X,
      full_labels = obj$labels,
      fit_fun = obj$fit_fun,
      B = FINAL_RESAMPLES,
      fraction = SUBSAMPLE_FRACTION,
      seed = MASTER_SEED + match(method_name, names(candidate_objects)) * 10000L
    )
  }
  
  candidate_objects[[method_name]]$stability <- stability
  candidate_objects[[method_name]]$fisher_p <- fisher_after_freeze(
    obj$labels, diet_final$case_status,
    seed = MASTER_SEED + match(method_name, names(candidate_objects))
  )
  
  if (obj$best_k == 1L) {
    cat("k=1; stability/Fisher not applicable\n")
  } else {
    cat(
      "mean Jaccard =", round(stability$mean_jaccard, 3),
      "| success =", round(100 * stability$success_rate, 1), "%\n"
    )
  }
}

# =============================================================================
# PHASE 8 — COMPACT DECISION AUDIT TABLE
# =============================================================================

cat("\n--- PHASE 8: DECISION AUDIT TABLE ---\n")

audit_table <- map_dfr(candidate_objects, function(obj) {
  criterion_digits <- if (obj$method == "GMM") 1L else 3L
  tibble(
    Method = obj$method,
    Space = obj$space,
    `Best k` = obj$best_k,
    `Criterion used` = obj$criterion,
    `Criterion value` = format_value(obj$criterion_value, criterion_digits),
    `Stability (Jaccard)` = format_value(obj$stability$mean_jaccard, 3L),
    `Fisher p` = format_value(obj$fisher_p, 3L)
  )
})

print(audit_table, n = Inf, width = Inf)
readr::write_csv(audit_table,
                 file.path(OUTPUT_DIR, "decision_audit.csv"))

# Detailed stability output for checking failed runs and weak clusters.
stability_details <- map_dfr(candidate_objects, function(obj) {
  cl_j <- obj$stability$cluster_jaccard
  tibble(
    Method = obj$method,
    Best_k = obj$best_k,
    Mean_Jaccard = obj$stability$mean_jaccard,
    Minimum_cluster_Jaccard = obj$stability$min_cluster_jaccard,
    Successful_resamples = obj$stability$valid_runs,
    Requested_resamples = obj$stability$B,
    Success_rate = obj$stability$success_rate,
    Cluster_Jaccards = if (all(is.na(cl_j))) "NA" else
      paste(names(cl_j), formatC(cl_j, digits = 3, format = "f"),
            sep = ":", collapse = "; ")
  )
})
readr::write_csv(stability_details,
                 file.path(OUTPUT_DIR, "stability_details.csv"))

# Cross-method ARI is an agreement diagnostic, not a selection criterion.
method_names <- names(candidate_objects)
ari_matrix <- matrix(
  NA_real_, nrow = length(method_names), ncol = length(method_names),
  dimnames = list(method_names, method_names)
)
for (i in seq_along(method_names)) {
  for (j in seq_along(method_names)) {
    li <- candidate_objects[[method_names[i]]]$labels
    lj <- candidate_objects[[method_names[j]]]$labels
    if (!anyNA(li) && !anyNA(lj)) {
      ari_matrix[i, j] <- mclust::adjustedRandIndex(li, lj)
    }
  }
}
readr::write_csv(
  as.data.frame(ari_matrix) %>% rownames_to_column("Method"),
  file.path(OUTPUT_DIR, "cross_method_ari.csv")
)

# =============================================================================
# PHASE 9 — SLIDE-READY TABLE FIGURE
# =============================================================================

navy <- "#1A1A2E"
row_fills <- rep(c("white", "#F3F5F7"), length.out = nrow(audit_table))

table_theme <- gridExtra::ttheme_minimal(
  core = list(
    fg_params = list(fontsize = 8.2, col = "#202020"),
    bg_params = list(fill = row_fills, col = NA)
  ),
  colhead = list(
    fg_params = list(fontsize = 8.2, fontface = "bold", col = navy),
    bg_params = list(fill = "#DCE6F1", col = NA)
  ),
  padding = unit(c(2.4, 2.4), "mm")
)

table_grob <- gridExtra::tableGrob(
  as.data.frame(audit_table),
  rows = NULL,
  theme = table_theme
)

title_grob <- textGrob(
  "Clustering decision audit",
  x = 0, hjust = 0,
  gp = gpar(fontsize = 11.5, fontface = "bold", col = navy)
)

footer_grob <- textGrob(
  paste0(
    "Criteria are method-specific and cannot be compared across rows. ",
    "Higher values are preferred for Gap, ICL and silhouette; lower values ",
    "are preferred for Xie-Beni. Stability is the mean best-matched cluster Jaccard ",
    "across ", FINAL_RESAMPLES, " repeated ",
    SUBSAMPLE_FRACTION * 100, "% subsamples. ",
    "Fisher p is post-selection and was not used for model selection. ",
    "FCM and UMAP candidate grids begin at k=2."
  ),
  x = 0,
  hjust = 0,
  gp = gpar(
    fontsize = 7.2,
    col = "grey35",
    fontface = "italic"
  )
)

slide_table_grob <- arrangeGrob(
  title_grob,
  table_grob,
  footer_grob,
  ncol = 1,
  heights = unit(c(0.48, 1, 0.74), c("cm", "null", "cm")),
  padding = unit(1.5, "mm")
)

pdf(file.path(OUTPUT_DIR, "slide5_decision_audit_table.pdf"),
    width = 7.1, height = 3.2)
grid.draw(slide_table_grob)
dev.off()

png(file.path(OUTPUT_DIR, "slide5_decision_audit_table.png"),
    width = 2130, height = 960, res = 300)
grid.draw(slide_table_grob)
dev.off()

# =============================================================================
# PHASE 10 — CANDIDATE PCA SCATTERS (NO FINAL METHOD SELECTED)
# =============================================================================

cat("\n--- PHASE 10: CANDIDATE PCA PROJECTIONS ---\n")

pc_plot_base <- tibble(
  PC1 = pca_full$x[, 1],
  PC2 = pca_full$x[, 2],
  case_status = factor(diet_final$case_status),
  participant_id = participant_id
)

for (method_name in names(candidate_objects)) {
  obj <- candidate_objects[[method_name]]
  if (obj$best_k == 1L || anyNA(obj$labels)) next
  
  plot_df <- pc_plot_base %>%
    mutate(cluster = factor(obj$labels))
  
  p <- ggplot(plot_df, aes(PC1, PC2)) +
    geom_point(
      aes(fill = cluster, shape = case_status),
      size = 2.6, alpha = 0.82, colour = "white", stroke = 0.4
    ) +
    scale_shape_manual(values = c("Cirrhosis" = 21, "HCC" = 24)) +
    guides(
      fill = guide_legend(
        override.aes = list(shape = 22, colour = "grey30", stroke = 0.3, size = 4)
      ),
      shape = guide_legend(
        override.aes = list(fill = "grey50", colour = "grey30", stroke = 0.3, size = 4)
      )
    ) +
    labs(
      title = paste0("PCA projection of the ", method_name, " candidate"),
      subtitle = paste0(
        "Audit candidate only — not yet selected as the final method | k = ",
        obj$best_k
      ),
      x = paste0("PC1 (", round(pc1_pct, 1), "% variance)"),
      y = paste0("PC2 (", round(pc2_pct, 1), "% variance)"),
      fill = "Candidate cluster",
      shape = "Disease status",
      caption = "Clusters were fitted in the method-specific space shown in the audit table; this PCA plot is only a common 2D projection."
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", colour = navy),
      plot.subtitle = element_text(size = 9.5, colour = "grey35"),
      plot.caption = element_text(size = 8, colour = "grey45", hjust = 0),
      panel.grid.minor = element_blank(),
      legend.position = "right"
    )
  
  stub <- safe_filename(method_name)
  ggsave(
    file.path(OUTPUT_DIR, paste0("candidate_pca_scatter_", stub, ".pdf")),
    p, width = 7.0, height = 5.1
  )
  ggsave(
    file.path(OUTPUT_DIR, paste0("candidate_pca_scatter_", stub, ".png")),
    p, width = 7.0, height = 5.1, dpi = 300
  )
}

# Save the complete audit object before stopping.
audit_object <- list(
  settings = list(
    seed = MASTER_SEED,
    K_MAX = K_MAX,
    pc_range = pc_range,
    pa_central = pa_central,
    GAP_B = GAP_B,
    tune_resamples = TUNE_RESAMPLES,
    final_resamples = FINAL_RESAMPLES,
    subsample_fraction = SUBSAMPLE_FRACTION,
    fcm_n_starts = FCM_N_STARTS
  ),
  audit_table = audit_table,
  stability_details = stability_details,
  cross_method_ari = ari_matrix,
  candidates = candidate_objects,
  tuning = list(
    ward_gap = gap_ward$details,
    kmeans = kmeans_candidates,
    gmm = gmm_candidates,
    fcm = fcm_results,
    umap = umap_results
  ),
  pca = pca_full
)
saveRDS(audit_object,
        file.path(OUTPUT_DIR, "clustering_audit_object.rds"))

# Save a pipeline-level copy of the full audit object.
save_rds_checked(audit_object, AUDIT_RDS)

cat("\nAudit outputs written to:", normalizePath(OUTPUT_DIR, winslash = "/"), "\n")
cat("Audit RDS saved to:", AUDIT_RDS, "\n")

# =============================================================================
# 11. COMMON CROSS-METHOD STABILITY COMPARISON
# =============================================================================

method_stability <- purrr::imap_dfr(candidate_objects, function(obj, method_name) {
  tibble::tibble(
    method = method_name,
    best_k = as.integer(obj$best_k),
    space = obj$space,
    native_criterion = obj$criterion,
    native_criterion_value = as.numeric(obj$criterion_value),
    mean_jaccard = as.numeric(obj$stability$mean_jaccard),
    minimum_cluster_jaccard = as.numeric(obj$stability$min_cluster_jaccard),
    resample_success_rate = as.numeric(obj$stability$success_rate),
    successful_resamples = as.integer(obj$stability$valid_runs),
    requested_resamples = as.integer(obj$stability$B),
    fisher_p_post_freeze = as.numeric(obj$fisher_p)
  )
}) %>%
  dplyr::mutate(
    eligible_multicluster = best_k > 1L &
      is.finite(minimum_cluster_jaccard) &
      is.finite(mean_jaccard)
  )

readr::write_csv(
  method_stability,
  file.path(OUTPUT_DIR, "cross_method_stability_comparison.csv")
)

cat("\n=== COMMON STABILITY COMPARISON ===\n")
print(method_stability, n = Inf, width = Inf)

# =============================================================================
# 12. FREEZE PRIMARY CLUSTERING METHOD
# =============================================================================
# Two supported rules:
#   configured -> preserve the method explicitly set in config.R after audit
#   stability  -> choose among frozen multicluster candidates by the common
#                 stability metric only: minimum cluster Jaccard, then mean
#                 Jaccard, then success rate. Native criteria have already
#                 selected the best configuration WITHIN each method.

if (CLUSTER_PRIMARY_SELECTION == "configured") {
  primary_method <- PRIMARY_CLUSTER_METHOD
  
  if (!primary_method %in% names(candidate_objects)) {
    stop(
      "PRIMARY_CLUSTER_METHOD must be one of: ",
      paste(names(candidate_objects), collapse = ", "),
      call. = FALSE
    )
  }
} else {
  eligible <- method_stability %>%
    dplyr::filter(eligible_multicluster) %>%
    dplyr::arrange(
      dplyr::desc(minimum_cluster_jaccard),
      dplyr::desc(mean_jaccard),
      dplyr::desc(resample_success_rate),
      best_k,
      method
    )
  
  if (nrow(eligible) == 0L) {
    stop(
      "No multicluster candidate had finite Jaccard stability; cannot use ",
      "CLUSTER_PRIMARY_SELECTION = 'stability'.",
      call. = FALSE
    )
  }
  
  primary_method <- eligible$method[[1]]
}

primary_obj <- candidate_objects[[primary_method]]

if (is.null(primary_obj) || is.null(primary_obj$labels)) {
  stop("Selected primary candidate has no cluster labels.", call. = FALSE)
}

if (!is.finite(primary_obj$best_k) || primary_obj$best_k < 2L) {
  stop(
    "Selected primary method produced k < 2. A one-cluster solution cannot ",
    "define dietary profiles.",
    call. = FALSE
  )
}

if (length(primary_obj$labels) != n_obs || anyNA(primary_obj$labels)) {
  stop("Primary cluster labels do not align with the preprocessing rows.", call. = FALSE)
}

if (
  primary_method == "FCM" &&
  isTRUE(CLUSTER_ASSERT_EXPECTED_FCM)
) {
  selected_k <- as.integer(primary_obj$best_k)
  selected_m <- as.numeric(primary_obj$m)
  
  if (selected_k != as.integer(FCM_K)) {
    stop(
      "FCM audit no longer selects the expected k = ", FCM_K,
      "; it selected k = ", selected_k,
      ". Review the audit before changing the frozen thesis solution.",
      call. = FALSE
    )
  }
  
  if (!is.finite(selected_m) || abs(selected_m - as.numeric(FCM_M)) > 1e-12) {
    stop(
      "FCM audit no longer selects the expected m = ", FCM_M,
      "; it selected m = ", selected_m,
      ". Review the audit before changing the frozen thesis solution.",
      call. = FALSE
    )
  }
}

raw_labels <- as.character(primary_obj$labels)
label_levels <- sort(unique(raw_labels))
hard_cluster <- match(raw_labels, label_levels)

if (length(unique(hard_cluster)) != as.integer(primary_obj$best_k)) {
  stop("Primary hard labels do not contain the expected number of clusters.", call. = FALSE)
}

# Maximum fuzzy membership is only defined for FCM here.
max_membership <- rep(NA_real_, n_obs)
primary_membership <- NULL

if (primary_method == "FCM" && !is.null(primary_obj$membership)) {
  primary_membership <- as.matrix(primary_obj$membership)
  
  if (nrow(primary_membership) == n_obs) {
    max_membership <- apply(primary_membership, 1, max)
    rownames(primary_membership) <- participant_id
  } else {
    warning(
      "FCM membership matrix does not align with participants; maximum ",
      "membership will be stored as NA.",
      call. = FALSE
    )
    primary_membership <- NULL
  }
}

cat("\n============================================================\n")
cat("PRIMARY CLUSTERING SOLUTION\n")
cat("============================================================\n")
cat("Method:", primary_method, "\n")
cat("k:", primary_obj$best_k, "\n")
cat("Space:", primary_obj$space, "\n")
cat("Native criterion:", primary_obj$criterion, "\n")
cat("Mean Jaccard:", round(primary_obj$stability$mean_jaccard, 3), "\n")
cat("Minimum-cluster Jaccard:", round(primary_obj$stability$min_cluster_jaccard, 3), "\n")

# =============================================================================
# 13. OUTCOME-INDEPENDENT PROFILE NAMING
# =============================================================================

pretty_food_labels <- function(x) {
  key <- stringr::str_remove(x, "_adj$")
  
  label_map <- c(
    veg = "Vegetables",
    legumes = "Legumes",
    fruit = "Fruit",
    grains_cereals = "Grains & cereals",
    dairy_and_alternatives = "Dairy & alternatives",
    red_meat = "Red meat",
    processed_meat = "Processed meat",
    poultry = "Poultry",
    seafood = "Seafood",
    eggs = "Eggs",
    nuts = "Nuts",
    discretionary_sweet = "Discretionary - sweet",
    discretionary_savoury = "Discretionary - savoury",
    alcohol = "Alcohol",
    saturated = "Saturated fats & oils",
    unsaturated = "Unsaturated fats & oils",
    healthy_drinks = "Water, tea & coffee",
    soy_plant = "Soy & plant-based foods"
  )
  
  out <- unname(label_map[key])
  missing <- is.na(out)
  out[missing] <- key[missing] %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_to_sentence()
  out
}

participant_scores <- tibble::as_tibble(X, .name_repair = "minimal") %>%
  dplyr::mutate(cluster = hard_cluster, .before = 1)

profile_means_raw <- participant_scores %>%
  tidyr::pivot_longer(
    cols = -cluster,
    names_to = "food_group",
    values_to = "z_score"
  ) %>%
  dplyr::group_by(cluster, food_group) %>%
  dplyr::summarise(
    mean_z = mean(z_score),
    .groups = "drop"
  )

diverse_positive <- c(
  "veg_adj", "fruit_adj", "legumes_adj", "nuts_adj", "seafood_adj",
  "eggs_adj", "healthy_drinks_adj", "unsaturated_adj", "soy_plant_adj",
  "grains_cereals_adj"
)

diverse_negative <- c(
  "discretionary_sweet_adj", "discretionary_savoury_adj",
  "processed_meat_adj", "alcohol_adj", "saturated_adj"
)

signature_foods <- intersect(
  c(diverse_positive, diverse_negative),
  unique(profile_means_raw$food_group)
)

diversity_scores <- profile_means_raw %>%
  dplyr::filter(food_group %in% signature_foods) %>%
  dplyr::mutate(
    direction = dplyr::if_else(food_group %in% diverse_positive, 1, -1)
  ) %>%
  dplyr::group_by(cluster) %>%
  dplyr::summarise(
    diversity_signature = mean(direction * mean_z),
    .groups = "drop"
  )

cluster_ids <- sort(unique(hard_cluster))

if (length(cluster_ids) == 2L && length(signature_foods) >= 6L) {
  food_diverse_cluster <- diversity_scores %>%
    dplyr::slice_max(diversity_signature, n = 1, with_ties = FALSE) %>%
    dplyr::pull(cluster)
  
  profile_map <- stats::setNames(
    rep("Food-limited", 2L),
    as.character(cluster_ids)
  )
  profile_map[as.character(food_diverse_cluster)] <- "Food-diverse"
} else {
  profile_map <- stats::setNames(
    paste0("Diet-profile-", seq_along(cluster_ids)),
    as.character(cluster_ids)
  )
}

profile_means <- profile_means_raw %>%
  dplyr::mutate(
    diet_profile = unname(profile_map[as.character(cluster)]),
    food_label = pretty_food_labels(food_group)
  )

cluster_name_mapping <- tibble::tibble(
  cluster = cluster_ids,
  diet_profile = unname(profile_map[as.character(cluster_ids)])
) %>%
  dplyr::left_join(diversity_scores, by = "cluster")

# =============================================================================
# 14. PARTICIPANT ASSIGNMENTS + PROFILE SUMMARIES
# =============================================================================

assignment_output <- tibble::tibble(
  study_id = participant_id,
  case_status = as.character(diet_final$case_status),
  clustering_method = primary_method,
  cluster = hard_cluster,
  diet_profile = unname(profile_map[as.character(hard_cluster)]),
  max_membership = max_membership
)

if (primary_method == "FCM") {
  assignment_output <- assignment_output %>%
    dplyr::mutate(
      fcm_cluster = cluster,
      fcm_profile = diet_profile,
      fcm_max_membership = max_membership
    )
}

profile_summary <- assignment_output %>%
  dplyr::group_by(diet_profile) %>%
  dplyr::summarise(
    n_total = dplyr::n(),
    n_HCC = sum(case_status == "HCC", na.rm = TRUE),
    n_Cirrhosis = sum(case_status == "Cirrhosis", na.rm = TRUE),
    HCC_percent = 100 * n_HCC / n_total,
    Cirrhosis_percent = 100 * n_Cirrhosis / n_total,
    mean_max_membership = if (all(is.na(max_membership))) {
      NA_real_
    } else {
      mean(max_membership, na.rm = TRUE)
    },
    .groups = "drop"
  )

profile_characteristics <- profile_means %>%
  dplyr::group_by(diet_profile) %>%
  dplyr::summarise(
    highest_mean_z_foods = paste(
      utils::head(food_label[order(mean_z, decreasing = TRUE)], 5),
      collapse = ", "
    ),
    lowest_mean_z_foods = paste(
      utils::head(food_label[order(mean_z, decreasing = FALSE)], 5),
      collapse = ", "
    ),
    .groups = "drop"
  )

profile_summary <- profile_summary %>%
  dplyr::left_join(profile_characteristics, by = "diet_profile")

# Post-freeze disease-status association for descriptive context only.
case_table <- table(
  assignment_output$diet_profile,
  assignment_output$case_status
)

primary_fisher <- if (
  nrow(case_table) >= 2L &&
  ncol(case_table) >= 2L &&
  all(rowSums(case_table) > 0L) &&
  all(colSums(case_table) > 0L)
) {
  stats::fisher.test(case_table)
} else {
  NULL
}

# =============================================================================
# 15. PRIMARY PROFILE HEATMAP
# =============================================================================

if (all(c("Food-limited", "Food-diverse") %in% profile_means$diet_profile)) {
  profile_levels_for_plot <- c("Food-limited", "Food-diverse")
  food_order <- profile_means %>%
    dplyr::filter(diet_profile == "Food-limited") %>%
    dplyr::arrange(mean_z) %>%
    dplyr::pull(food_label)
} else {
  profile_levels_for_plot <- cluster_name_mapping$diet_profile
  reference_profile <- profile_levels_for_plot[[1]]
  food_order <- profile_means %>%
    dplyr::filter(diet_profile == reference_profile) %>%
    dplyr::arrange(mean_z) %>%
    dplyr::pull(food_label)
}

heatmap_data <- profile_means %>%
  dplyr::mutate(
    diet_profile = factor(diet_profile, levels = profile_levels_for_plot),
    food_label = factor(food_label, levels = rev(food_order)),
    cell_label = sprintf("%+.2f", mean_z)
  )

z_limit <- max(abs(heatmap_data$mean_z), na.rm = TRUE)
z_limit <- max(0.5, ceiling(z_limit * 10) / 10)

p_primary_heatmap <- ggplot2::ggplot(
  heatmap_data,
  ggplot2::aes(x = diet_profile, y = food_label, fill = mean_z)
) +
  ggplot2::geom_tile(colour = "white", linewidth = 1.0) +
  ggplot2::geom_text(
    ggplot2::aes(label = cell_label),
    fontface = "bold",
    size = 3.6
  ) +
  ggplot2::scale_fill_gradient2(
    low = "#7F9CCB",
    mid = "#F7F7F7",
    high = "#C94D4D",
    midpoint = 0,
    limits = c(-z_limit, z_limit),
    oob = scales::squish,
    name = "Mean z"
  ) +
  ggplot2::labs(
    title = paste0("Food-group profiles: ", primary_method),
    subtitle = paste0("Frozen primary solution; k = ", primary_obj$best_k),
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 13) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    panel.grid = ggplot2::element_blank(),
    axis.text.x = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(FIGURE_DIR, "primary_profile_heatmap.png"),
  p_primary_heatmap,
  width = 9.2,
  height = 10.0,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIGURE_DIR, "primary_profile_heatmap.pdf"),
  p_primary_heatmap,
  width = 9.2,
  height = 10.0,
  bg = "white"
)

# =============================================================================
# 15B. PCA SCATTER PLOT — CLUSTER ASSIGNMENTS IN PC1/PC2
# =============================================================================
# Colour = dietary profile, Shape = disease status (circle/triangle).
# PCA is run on the same adj_scaled matrix used for clustering (mean-centred,
# not variance-scaled to match the clustering input space).
# =============================================================================

cluster_pca <- stats::prcomp(X, center = TRUE, scale. = FALSE)
cluster_pc_var <- cluster_pca$sdev^2 / sum(cluster_pca$sdev^2)

pca_scatter_data <- tibble::tibble(
  study_id = participant_id,
  PC1 = cluster_pca$x[, 1],
  PC2 = cluster_pca$x[, 2]
) %>%
  dplyr::left_join(
    assignment_output %>%
      dplyr::select(study_id, diet_profile, case_status, max_membership),
    by = "study_id"
  ) %>%
  dplyr::mutate(
    diet_profile = factor(
      diet_profile,
      levels = profile_levels_for_plot
    ),
    case_status = factor(
      case_status,
      levels = c("Cirrhosis", "HCC")
    )
  )

# For FCM: use max_membership as alpha (faded = ambiguous assignment).
# For non-FCM methods: max_membership is NA, so use fixed alpha.
is_fcm <- primary_method == "FCM" &&
  !all(is.na(pca_scatter_data$max_membership))

cluster_colours <- c(
  "Food-limited" = "#7F9CCB",
  "Food-diverse" = "#C94D4D"
)

# Fall back to a generic palette when profile names are non-standard.
if (!all(profile_levels_for_plot %in% names(cluster_colours))) {
  n_profiles <- length(profile_levels_for_plot)
  cluster_colours <- stats::setNames(
    scales::hue_pal()(n_profiles),
    profile_levels_for_plot
  )
}

if (is_fcm) {
  p_pca_scatter <- ggplot2::ggplot(
    pca_scatter_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      colour = diet_profile,
      shape = case_status,
      alpha = max_membership
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_point(size = 3.0) +
    ggplot2::scale_alpha_continuous(
      name = "Maximum FCM\nmembership",
      range = c(0.15, 1.0),
      limits = c(0.50, 1.0),
      breaks = seq(0.50, 1.00, by = 0.10),
      labels = sprintf("%.2f", seq(0.50, 1.00, by = 0.10)),
      guide = ggplot2::guide_legend(
        override.aes = list(
          shape = 16,
          colour = "grey20",
          size = 3.5
        )
      )
    ) +
    ggplot2::scale_colour_manual(
      values = cluster_colours,
      name = "Dietary profile"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Cirrhosis" = 16, "HCC" = 17),
      name = "Disease status"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, shape = 16)
      ),
      shape = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, colour = "grey30")
      )
    )
} else {
  p_pca_scatter <- ggplot2::ggplot(
    pca_scatter_data,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      colour = diet_profile,
      shape = case_status
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_point(size = 3.0, alpha = 0.80) +
    ggplot2::scale_colour_manual(
      values = cluster_colours,
      name = "Dietary profile"
    ) +
    ggplot2::scale_shape_manual(
      values = c("Cirrhosis" = 16, "HCC" = 17),
      name = "Disease status"
    ) +
    ggplot2::guides(
      colour = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, shape = 16)
      ),
      shape = ggplot2::guide_legend(
        override.aes = list(size = 4, alpha = 1, colour = "grey30")
      )
    )
}

fcm_subtitle <- if (is_fcm) {
  paste0(
    " | Opacity = max FCM membership (",
    sprintf("%.0f", 100 * min(pca_scatter_data$max_membership, na.rm = TRUE)),
    "%-",
    sprintf("%.0f", 100 * max(pca_scatter_data$max_membership, na.rm = TRUE)),
    "%)"
  )
} else {
  ""
}

p_pca_scatter <- p_pca_scatter +
  ggplot2::labs(
    title = paste0("Dietary-pattern clustering in PCA space (", primary_method, ")"),
    subtitle = paste0(
      "k = ", primary_obj$best_k,
      " | n = ", nrow(pca_scatter_data),
      " | Mean-centred food-group z-scores",
      fcm_subtitle
    ),
    x = paste0("PC1 (", round(100 * cluster_pc_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(100 * cluster_pc_var[2], 1), "% variance)")
  ) +
  ggplot2::theme_minimal(base_size = 12) +
  ggplot2::theme(
    panel.grid.minor = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "right"
  )

ggplot2::ggsave(
  file.path(FIGURE_DIR, "primary_cluster_pca_scatter.png"),
  p_pca_scatter,
  width = 9.0,
  height = 7.0,
  dpi = 300,
  bg = "white"
)

ggplot2::ggsave(
  file.path(FIGURE_DIR, "primary_cluster_pca_scatter.pdf"),
  p_pca_scatter,
  width = 9.0,
  height = 7.0,
  bg = "white"
)

# =============================================================================
# 16. EXPORT FROZEN PRIMARY OUTPUTS
# =============================================================================

readr::write_csv(
  assignment_output,
  file.path(OUTPUT_DIR, "participant_assignments.csv")
)

readr::write_csv(
  profile_means %>%
    dplyr::select(cluster, diet_profile, food_group, food_label, mean_z) %>%
    dplyr::arrange(diet_profile, mean_z),
  file.path(OUTPUT_DIR, "primary_profile_mean_z_scores.csv")
)

readr::write_csv(
  profile_summary,
  file.path(OUTPUT_DIR, "primary_profile_summary.csv")
)

readr::write_csv(
  cluster_name_mapping,
  file.path(OUTPUT_DIR, "cluster_name_mapping.csv")
)

primary_selection_record <- method_stability %>%
  dplyr::filter(method == primary_method) %>%
  dplyr::mutate(
    selection_mode = CLUSTER_PRIMARY_SELECTION,
    configured_method = if (CLUSTER_PRIMARY_SELECTION == "configured") {
      PRIMARY_CLUSTER_METHOD
    } else {
      NA_character_
    }
  )

readr::write_csv(
  primary_selection_record,
  file.path(OUTPUT_DIR, "primary_selection_record.csv")
)

clustering_output <- list(
  assignments = assignment_output,
  profile_means = profile_means,
  profile_summary = profile_summary,
  cluster_name_mapping = cluster_name_mapping,
  primary_method = primary_method,
  primary_candidate = list(
    method = primary_method,
    space = primary_obj$space,
    best_k = primary_obj$best_k,
    criterion = primary_obj$criterion,
    criterion_value = primary_obj$criterion_value,
    stability = primary_obj$stability,
    fisher_p_post_freeze = primary_obj$fisher_p,
    m = primary_obj$m %||% NULL,
    n_neighbors = primary_obj$n_neighbors %||% NULL,
    min_dist = primary_obj$min_dist %||% NULL
  ),
  membership = primary_membership,
  audit = list(
    audit_rds = AUDIT_RDS,
    decision_audit = audit_table,
    method_stability = method_stability,
    stability_details = stability_details,
    cross_method_ari = ari_matrix,
    pca_parallel_analysis_n = pa_central,
    pca_candidate_range = pc_range,
    pca_variance = pc_variance_full
  ),
  post_freeze = list(
    case_table = case_table,
    fisher_test = primary_fisher
  ),
  settings = list(
    seed = MASTER_SEED,
    primary_selection = CLUSTER_PRIMARY_SELECTION,
    primary_method = primary_method,
    K_MAX = K_MAX,
    PC_MAX_CAP = PC_MAX_CAP,
    PA_ITER = PA_ITER,
    GAP_B = GAP_B,
    tune_resamples = TUNE_RESAMPLES,
    final_resamples = FINAL_RESAMPLES,
    subsample_fraction = SUBSAMPLE_FRACTION,
    fcm_m_values = FCM_M_VALUES,
    fcm_n_starts = FCM_N_STARTS,
    umap_neighbors = UMAP_NEIGHBORS,
    umap_min_dist = UMAP_MIN_DIST
  )
)

save_rds_checked(
  clustering_output,
  RDS_PATHS$clustering
)

cat("\n============================================================\n")
cat("CLUSTERING COMPLETE\n")
cat("============================================================\n")
cat("Primary method:", primary_method, "\n")
cat("Primary k:", primary_obj$best_k, "\n")
if (primary_method == "FCM") {
  cat("Primary FCM m:", primary_obj$m, "\n")
}
cat("Assignments:", nrow(assignment_output), "participants\n")
cat("\nProfile summary:\n")
print(profile_summary, width = Inf)
cat("\nProfile x disease-status table (post-freeze only):\n")
print(case_table)
if (!is.null(primary_fisher)) {
  cat("Post-freeze Fisher p =", signif(primary_fisher$p.value, 4), "\n")
}
cat("\nSaved pipeline RDS:", RDS_PATHS$clustering, "\n")
cat("Saved audit RDS:", AUDIT_RDS, "\n")
cat("Results directory:", normalizePath(OUTPUT_DIR, winslash = "/"), "\n")
cat("Figures directory:", normalizePath(FIGURE_DIR, winslash = "/"), "\n")
cat("============================================================\n")