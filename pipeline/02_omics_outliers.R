# =============================================================================
# 02_stage1_screen.R
# Omics preprocessing + Stage 1 multivariate outlier screening
#
# PURPOSE
#     * Visit 1 only
#     * Cohort A = preprocessed diet cohort with usable Visit-1 data in all
#       three omics panels
#     * Drop features with >30% missingness in Cohort A
#     * Force log1p for panels listed in FORCE_LOG_PANELS (cytokines)
#     * Median-impute after transformation
#     * Remove zero-variance features
#     * Mean-centred, NOT variance-scaled PCA
#     * Retain PCs using parallel analysis
#     * Hotelling T-squared screening at the configured confidence level
#     * NO Q-residual exclusion
#
# INPUTS
#   RDS_PATHS$preprocess
#   Raw omics files defined in OMICS_PANELS
#
# OUTPUT
#   RDS_PATHS$stage1
#
# NOTE
#   Processed Cohort A data are used for within-diet-profile Conditions A/B.
#   Processed Cohort B data retain the full Visit-1 all-omics reference cohort
#   for downstream Condition C. The same feature set and log rule are used,
#   but median imputation is performed within each cohort, matching the
#   existing analysis logic.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(psych)
})

source("pipeline/config.R")
source("pipeline/utils.R")

set.seed(MASTER_SEED)
ensure_parent_dir(RDS_PATHS$stage1)

required_stage1_settings <- c(
  "STAGE1_MAX_FEATURE_MISSING",
  "STAGE1_PA_ITER",
  "STAGE1_HOTELLING_CONFIDENCE"
)

missing_stage1_settings <- required_stage1_settings[
  !vapply(
    required_stage1_settings,
    exists,
    logical(1),
    inherits = TRUE
  )
]

if (length(missing_stage1_settings) > 0L) {
  stop(
    "Missing Stage 1 setting(s) in config.R: ",
    paste(missing_stage1_settings, collapse = ", "),
    ". Add: STAGE1_MAX_FEATURE_MISSING <- 0.30; ",
    "STAGE1_PA_ITER <- 1000L; ",
    "STAGE1_HOTELLING_CONFIDENCE <- 0.99.",
    call. = FALSE
  )
}

if (
  !is.numeric(STAGE1_MAX_FEATURE_MISSING) ||
  length(STAGE1_MAX_FEATURE_MISSING) != 1L ||
  STAGE1_MAX_FEATURE_MISSING < 0 ||
  STAGE1_MAX_FEATURE_MISSING >= 1
) {
  stop("STAGE1_MAX_FEATURE_MISSING must be in [0, 1).", call. = FALSE)
}

if (
  !is.numeric(STAGE1_HOTELLING_CONFIDENCE) ||
  length(STAGE1_HOTELLING_CONFIDENCE) != 1L ||
  STAGE1_HOTELLING_CONFIDENCE <= 0 ||
  STAGE1_HOTELLING_CONFIDENCE >= 1
) {
  stop("STAGE1_HOTELLING_CONFIDENCE must be in (0, 1).", call. = FALSE)
}

if (STAGE1_PA_ITER < 100L) {
  warning(
    "STAGE1_PA_ITER is below 100; parallel analysis may be unstable.",
    call. = FALSE
  )
}

# =============================================================================
# 1. LOAD THE DIET COHORT
# =============================================================================

preprocess <- read_rds_checked(
  RDS_PATHS$preprocess,
  "preprocessed diet RDS"
)

diet_final <- preprocess$diet_final %||%
  preprocess$data$diet_final

if (is.null(diet_final)) {
  stop(
    "Preprocessing RDS does not contain diet_final.",
    call. = FALSE
  )
}

diet_id_col <- if ("id" %in% names(diet_final)) {
  "id"
} else if ("study_id" %in% names(diet_final)) {
  "study_id"
} else {
  stop(
    "diet_final must contain either 'id' or 'study_id'.",
    call. = FALSE
  )
}

assert_columns(
  diet_final,
  c(diet_id_col, "case_status"),
  "diet_final"
)

diet_final <- diet_final %>%
  dplyr::mutate(
    .stage1_id = normalise_id(.data[[diet_id_col]])
  )

if (anyNA(diet_final$.stage1_id)) {
  stop("diet_final contains missing participant IDs.", call. = FALSE)
}

if (anyDuplicated(diet_final$.stage1_id)) {
  stop("diet_final contains duplicated participant IDs.", call. = FALSE)
}

diet_ids <- diet_final$.stage1_id

# =============================================================================
# 2. LOAD ALL OMICS PANELS STRUCTURALLY
# =============================================================================

panel_names <- names(OMICS_PANELS)

if (!setequal(
  panel_names,
  c("cytokines", "small_molecules", "lipoproteins")
)) {
  warning(
    "OMICS_PANELS differs from the expected three Aim 2 panels: ",
    paste(panel_names, collapse = ", "),
    call. = FALSE
  )
}

panel_loaded <- stats::setNames(
  lapply(
    panel_names,
    function(panel_name) {
      obj <- load_omics_panel(
        panel_name = panel_name,
        visit = OMICS_VISIT
      )
      
      assert_unique_ids(
        obj$data,
        id_col = "study_id",
        label = paste0(obj$label, " Visit ", OMICS_VISIT)
      )
      
      obj
    }
  ),
  panel_names
)

# Raw Visit-1 IDs are also collected before all-missing panel rows are removed.
# This preserves the full reference cohort definition used by Condition C.
get_raw_visit_ids <- function(panel_name) {
  spec <- get_panel_spec(panel_name)
  
  read_ids <- function(path, id_col) {
    assert_file_exists(path)
    
    dat <- readr::read_csv(
      path,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
    
    assert_columns(
      dat,
      c(id_col, spec$visit_col),
      paste0(spec$label, " raw Visit-1 ID source")
    )
    
    visit_value <- suppressWarnings(
      as.integer(dat[[spec$visit_col]])
    )
    
    ids <- dat[[id_col]][visit_value == OMICS_VISIT]
    ids <- normalise_id(ids)
    unique(ids[!is.na(ids)])
  }
  
  unique(
    c(
      read_ids(spec$cirrhosis_file, spec$cirrhosis_id_col),
      read_ids(spec$hcc_file, spec$hcc_id_col)
    )
  )
}

raw_visit_id_sets <- stats::setNames(
  lapply(panel_names, get_raw_visit_ids),
  panel_names
)

usable_visit_id_sets <- lapply(
  panel_loaded,
  function(obj) unique(obj$data$study_id)
)

# Cohort B: raw Visit-1 presence in all 3 omics panels.
cohort_b_ids <- Reduce(
  intersect,
  raw_visit_id_sets
)

# Cohort A: LER-excluded diet cohort AND a usable (not entirely missing)
# Visit-1 profile in every omics panel.
usable_all_three <- Reduce(
  intersect,
  usable_visit_id_sets
)

cohort_a_ids <- diet_ids[
  diet_ids %in% usable_all_three
]

if (length(cohort_a_ids) < 10L) {
  stop(
    "Cohort A has fewer than 10 participants after matching diet and omics.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("STAGE 1 COHORTS\n")
cat("============================================================\n")
cat("Diet cohort from 00_preprocess:", length(diet_ids), "\n")
cat("Cohort A (diet + usable all 3 omics):", length(cohort_a_ids), "\n")
cat("Cohort B (raw Visit-1 presence in all 3 omics):", length(cohort_b_ids), "\n")

cohort_a_status <- diet_final %>%
  dplyr::filter(.stage1_id %in% cohort_a_ids) %>%
  dplyr::count(case_status, name = "n")

print(cohort_a_status)

# =============================================================================
# 3. HELPERS FOR FROZEN PANEL PREPROCESSING
# =============================================================================

replace_nonfinite_with_na <- function(data, feature_cols) {
  out <- data
  
  for (feature in feature_cols) {
    x <- suppressWarnings(as.numeric(out[[feature]]))
    x[!is.finite(x)] <- NA_real_
    out[[feature]] <- x
  }
  
  out
}

median_impute_with_record <- function(data, feature_cols, label) {
  assert_columns(data, feature_cols, label)
  
  out <- data
  
  medians <- vapply(
    feature_cols,
    function(feature) {
      x <- suppressWarnings(as.numeric(out[[feature]]))
      
      if (all(is.na(x))) {
        stop(
          label,
          ": retained feature '",
          feature,
          "' is entirely missing.",
          call. = FALSE
        )
      }
      
      stats::median(x, na.rm = TRUE)
    },
    numeric(1)
  )
  
  for (feature in feature_cols) {
    x <- suppressWarnings(as.numeric(out[[feature]]))
    x[is.na(x)] <- medians[[feature]]
    out[[feature]] <- x
  }
  
  list(
    data = out,
    medians = medians
  )
}

subset_panel_in_id_order <- function(data, ids) {
  ids <- normalise_id(ids)
  ids_present <- ids[ids %in% data$study_id]
  
  idx <- match(ids_present, data$study_id)
  data[idx, , drop = FALSE]
}

apply_panel_transform <- function(data, feature_cols, panel_name) {
  transformed <- apply_configured_log1p(
    data = data,
    feature_cols = feature_cols,
    panel_name = panel_name
  )
  
  transformed
}

run_parallel_analysis <- function(X, panel_name, seed) {
  set.seed(seed)
  
  pa <- suppressWarnings(
    psych::fa.parallel(
      X,
      fa = "pc",
      n.iter = STAGE1_PA_ITER,
      cor = "cov",
      plot = FALSE,
      show.legend = FALSE,
      error.bars = FALSE
    )
  )
  
  n_components <- suppressWarnings(
    as.integer(pa$ncomp)
  )
  
  max_allowed <- min(
    ncol(X),
    nrow(X) - 2L
  )
  
  if (!is.finite(n_components) || n_components < 1L) {
    warning(
      panel_name,
      ": parallel analysis retained <1 PC; using 1 PC.",
      call. = FALSE
    )
    n_components <- 1L
  }
  
  n_components <- min(
    n_components,
    max_allowed
  )
  
  if (n_components < 1L) {
    stop(
      panel_name,
      ": too few observations for Hotelling T-squared.",
      call. = FALSE
    )
  }
  
  list(
    n_components = n_components,
    object = pa
  )
}

hotelling_t2_from_pca <- function(pca, n_components, confidence) {
  scores <- pca$x[
    ,
    seq_len(n_components),
    drop = FALSE
  ]
  
  eigenvalues <- pca$sdev[
    seq_len(n_components)
  ]^2
  
  if (
    any(!is.finite(eigenvalues)) ||
    any(eigenvalues <= sqrt(.Machine$double.eps))
  ) {
    stop(
      "Non-positive PCA eigenvalue encountered in retained components.",
      call. = FALSE
    )
  }
  
  contributions <- sweep(
    scores^2,
    MARGIN = 2L,
    STATS = eigenvalues,
    FUN = "/"
  )
  
  t2 <- rowSums(contributions)
  
  n <- nrow(scores)
  A <- n_components
  
  if (n <= A) {
    stop(
      "Hotelling T-squared requires n > retained PCs.",
      call. = FALSE
    )
  }
  
  # F-based PCA Hotelling T-squared control limit.
  # The (n + 1) / n factor corresponds to the standard individual-observation
  # PCA monitoring limit.
  cutoff <- (
    A * (n + 1) * (n - 1) /
      (n * (n - A))
  ) * stats::qf(
    confidence,
    df1 = A,
    df2 = n - A
  )
  
  list(
    t2 = t2,
    cutoff = cutoff,
    contributions = contributions
  )
}

# =============================================================================
# 4. PROCESS EACH PANEL + RUN STAGE 1 T-SQUARED SCREEN
# =============================================================================

panel_results <- list()

for (panel_index in seq_along(panel_names)) {
  panel_name <- panel_names[[panel_index]]
  panel_obj  <- panel_loaded[[panel_name]]
  
  cat("\n============================================================\n")
  cat("PANEL:", panel_obj$label, "\n")
  cat("============================================================\n")
  
  # ---------------------------------------------------------------------------
  # 4a. Cohort A raw panel
  # ---------------------------------------------------------------------------
  
  panel_a_raw <- subset_panel_in_id_order(
    panel_obj$data,
    cohort_a_ids
  )
  
  if (nrow(panel_a_raw) != length(cohort_a_ids)) {
    stop(
      panel_obj$label,
      ": not all Cohort A participants were found after panel loading.",
      call. = FALSE
    )
  }
  
  all_features <- panel_obj$feature_cols
  
  panel_a_raw <- replace_nonfinite_with_na(
    panel_a_raw,
    all_features
  )
  
  feature_missingness <- tibble::tibble(
    feature = all_features,
    missing_fraction = vapply(
      all_features,
      function(feature) {
        mean(is.na(panel_a_raw[[feature]]))
      },
      numeric(1)
    )
  )
  
  features_dropped_missing <- feature_missingness %>%
    dplyr::filter(
      missing_fraction > STAGE1_MAX_FEATURE_MISSING
    ) %>%
    dplyr::pull(feature)
  
  candidate_features <- setdiff(
    all_features,
    features_dropped_missing
  )
  
  if (length(candidate_features) < 2L) {
    stop(
      panel_obj$label,
      ": fewer than 2 features remain after missingness filtering.",
      call. = FALSE
    )
  }
  
  cat("Raw features:", length(all_features), "\n")
  cat(
    "Dropped for >",
    STAGE1_MAX_FEATURE_MISSING * 100,
    "% missingness:",
    length(features_dropped_missing),
    "\n"
  )
  
  # ---------------------------------------------------------------------------
  # 4b. Apply configured transformation BEFORE imputation
  # ---------------------------------------------------------------------------
  
  transform_a <- apply_panel_transform(
    data = panel_a_raw,
    feature_cols = candidate_features,
    panel_name = panel_name
  )
  
  panel_a_transformed <- transform_a$data
  log_applied <- transform_a$log_applied
  
  # ---------------------------------------------------------------------------
  # 4c. Median-impute Cohort A
  # ---------------------------------------------------------------------------
  
  imputed_a <- median_impute_with_record(
    data = panel_a_transformed,
    feature_cols = candidate_features,
    label = paste0(panel_obj$label, " Cohort A")
  )
  
  panel_a_imputed <- imputed_a$data
  
  # Remove features with invalid / zero variance AFTER the frozen
  # transform and Cohort A imputation.
  feature_sds <- vapply(
    panel_a_imputed[, candidate_features, drop = FALSE],
    stats::sd,
    numeric(1)
  )
  
  features_dropped_zero_variance <- names(feature_sds)[
    !is.finite(feature_sds) |
      feature_sds <= sqrt(.Machine$double.eps)
  ]
  
  retained_features <- setdiff(
    candidate_features,
    features_dropped_zero_variance
  )
  
  if (length(retained_features) < 2L) {
    stop(
      panel_obj$label,
      ": fewer than 2 non-zero-variance features remain.",
      call. = FALSE
    )
  }
  
  features_dropped <- setdiff(
    all_features,
    retained_features
  )
  
  processed_a <- panel_a_imputed %>%
    dplyr::select(
      study_id,
      case_status,
      dplyr::all_of(retained_features)
    )
  
  if (
    any(
      !is.finite(
        as.matrix(
          processed_a[, retained_features, drop = FALSE]
        )
      )
    )
  ) {
    stop(
      panel_obj$label,
      ": non-finite values remain in processed Cohort A data.",
      call. = FALSE
    )
  }
  
  cat("Log1p applied:", log_applied, "\n")
  cat(
    "Dropped for zero variance:",
    length(features_dropped_zero_variance),
    "\n"
  )
  cat("Retained features:", length(retained_features), "\n")
  
  # ---------------------------------------------------------------------------
  # 4d. Freeze corresponding Cohort B processed data
  # ---------------------------------------------------------------------------
  
  panel_b_raw <- subset_panel_in_id_order(
    panel_obj$data,
    cohort_b_ids
  )
  
  panel_b_raw <- replace_nonfinite_with_na(
    panel_b_raw,
    retained_features
  )
  
  transform_b <- apply_panel_transform(
    data = panel_b_raw,
    feature_cols = retained_features,
    panel_name = panel_name
  )
  
  if (!identical(transform_b$log_applied, log_applied)) {
    stop(
      panel_obj$label,
      ": inconsistent log decision between Cohort A and Cohort B.",
      call. = FALSE
    )
  }
  
  imputed_b <- median_impute_with_record(
    data = transform_b$data,
    feature_cols = retained_features,
    label = paste0(panel_obj$label, " Cohort B")
  )
  
  processed_b <- imputed_b$data %>%
    dplyr::select(
      study_id,
      case_status,
      dplyr::all_of(retained_features)
    )
  
  # ---------------------------------------------------------------------------
  # 4e. Mean-centred PCA on Cohort A (NO variance scaling)
  # ---------------------------------------------------------------------------
  
  X_a <- as.matrix(
    processed_a[, retained_features, drop = FALSE]
  )
  storage.mode(X_a) <- "double"
  rownames(X_a) <- processed_a$study_id
  
  pca <- stats::prcomp(
    X_a,
    center = TRUE,
    scale. = FALSE
  )
  
  pa <- run_parallel_analysis(
    X = X_a,
    panel_name = panel_obj$label,
    seed = MASTER_SEED + 1000L * panel_index
  )
  
  n_components <- pa$n_components
  
  # ---------------------------------------------------------------------------
  # 4f. Hotelling T-squared at configured confidence; NO Q exclusion
  # ---------------------------------------------------------------------------
  
  t2_result <- hotelling_t2_from_pca(
    pca = pca,
    n_components = n_components,
    confidence = STAGE1_HOTELLING_CONFIDENCE
  )
  
  flagged <- t2_result$t2 > t2_result$cutoff
  
  flagged_ids <- processed_a$study_id[flagged]
  
  t2_table <- tibble::tibble(
    study_id = processed_a$study_id,
    case_status = processed_a$case_status,
    hotelling_t2 = as.numeric(t2_result$t2),
    t2_cutoff = as.numeric(t2_result$cutoff),
    flagged = flagged
  ) %>%
    dplyr::arrange(dplyr::desc(hotelling_t2))
  
  contribution_table <- tibble::as_tibble(
    t2_result$contributions,
    .name_repair = "minimal"
  )
  
  names(contribution_table) <- paste0(
    "PC",
    seq_len(n_components)
  )
  
  contribution_table <- contribution_table %>%
    dplyr::mutate(
      study_id = processed_a$study_id,
      .before = 1
    )
  
  pca_loadings <- tibble::as_tibble(
    pca$rotation[
      ,
      seq_len(n_components),
      drop = FALSE
    ],
    rownames = "feature"
  )
  
  pca_variance <- tibble::tibble(
    component = seq_along(pca$sdev),
    eigenvalue = pca$sdev^2,
    variance_fraction = (
      pca$sdev^2 /
        sum(pca$sdev^2)
    ),
    cumulative_variance_fraction = cumsum(
      pca$sdev^2 /
        sum(pca$sdev^2)
    ),
    retained = seq_along(pca$sdev) <= n_components
  )
  
  cat("Parallel-analysis PCs retained:", n_components, "\n")
  cat(
    "Hotelling T2 confidence:",
    STAGE1_HOTELLING_CONFIDENCE,
    "\n"
  )
  cat("Hotelling T2 cutoff:", signif(t2_result$cutoff, 5), "\n")
  cat("Flagged participants:", length(flagged_ids), "\n")
  
  if (length(flagged_ids) > 0L) {
    cat("  ", paste(flagged_ids, collapse = ", "), "\n", sep = "")
  }
  
  # ---------------------------------------------------------------------------
  # 4g. Freeze panel audit result
  # ---------------------------------------------------------------------------
  
  panel_results[[panel_name]] <- list(
    panel_name = panel_name,
    label = panel_obj$label,
    
    # Compatibility fields expected by the existing biomarker script.
    features_dropped = features_dropped,
    log_applied = log_applied,
    flagged_ids = flagged_ids,
    
    # Explicit frozen recipe.
    all_features = all_features,
    retained_features = retained_features,
    features_dropped_missing = features_dropped_missing,
    features_dropped_zero_variance = features_dropped_zero_variance,
    feature_missingness = feature_missingness,
    
    imputation_medians_cohort_a = imputed_a$medians[
      retained_features
    ],
    imputation_medians_cohort_b = imputed_b$medians[
      retained_features
    ],
    
    processed_cohort_a = processed_a,
    processed_cohort_b = processed_b,
    
    pca = list(
      model = pca,
      n_components_parallel = n_components,
      loadings = pca_loadings,
      variance = pca_variance,
      centered = TRUE,
      scaled = FALSE
    ),
    
    hotelling = list(
      confidence = STAGE1_HOTELLING_CONFIDENCE,
      cutoff = t2_result$cutoff,
      table = t2_table,
      pc_contributions = contribution_table,
      q_residual_exclusion = FALSE
    ),
    
    counts = list(
      n_cohort_a = nrow(processed_a),
      n_cohort_b_usable = nrow(processed_b),
      n_features_raw = length(all_features),
      n_features_retained = length(retained_features),
      n_flagged = length(flagged_ids)
    )
  )
}

# =============================================================================
# 5. SAVE THE STAGE 1 AUDIT OBJECT
# =============================================================================

stage1_output <- list(
  cohort_a_ids = cohort_a_ids,
  cohort_b_ids = cohort_b_ids,
  
  cohort = list(
    diet_ids = diet_ids,
    usable_all_three_ids = usable_all_three,
    cohort_a_ids = cohort_a_ids,
    cohort_b_ids = cohort_b_ids,
    cohort_a_status = cohort_a_status,
    raw_visit_id_sets = raw_visit_id_sets,
    usable_visit_id_sets = usable_visit_id_sets
  ),
  
  panel_results = panel_results,
  
  settings = list(
    seed = MASTER_SEED,
    omics_visit = OMICS_VISIT,
    max_feature_missing = STAGE1_MAX_FEATURE_MISSING,
    imputation = OMICS_IMPUTATION,
    force_log_panels = FORCE_LOG_PANELS,
    pca_center = TRUE,
    pca_scale = FALSE,
    parallel_analysis_iterations = STAGE1_PA_ITER,
    hotelling_confidence = STAGE1_HOTELLING_CONFIDENCE,
    q_residual_exclusion = FALSE
  )
)

save_rds_checked(
  stage1_output,
  RDS_PATHS$stage1
)

cat("\n============================================================\n")
cat("STAGE 1 SCREEN COMPLETE\n")
cat("Saved:", RDS_PATHS$stage1, "\n")
cat("============================================================\n")

summary_table <- purrr::imap_dfr(
  panel_results,
  function(x, panel_name) {
    tibble::tibble(
      panel = panel_name,
      n_cohort_a = x$counts$n_cohort_a,
      n_cohort_b_usable = x$counts$n_cohort_b_usable,
      n_features_raw = x$counts$n_features_raw,
      n_features_retained = x$counts$n_features_retained,
      log1p = x$log_applied,
      n_pcs = x$pca$n_components_parallel,
      n_t2_flagged = x$counts$n_flagged
    )
  }
)

print(summary_table, n = Inf, width = Inf)