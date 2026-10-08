# =============================================================================
# 04_common.R
# Shared setup for all 04x deep-dive sub-scripts.
#
# Each 04x script sources this file AFTER sourcing config.R and utils.R.
# This file loads frozen upstream objects, defines shared helpers, and sets
# output directories. It does NOT perform any analysis.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

if (!requireNamespace("vegan", quietly = TRUE)) {
  stop(
    "Package 'vegan' is required for PERMANOVA/PERMDISP. Install it with ",
    "install.packages('vegan').",
    call. = FALSE
  )
}

# =============================================================================
# 1. OPTIONAL CONFIG SETTINGS WITH DEFAULTS
# =============================================================================

DEEP_DIVE_PERMUTATIONS <- if (exists(
  "DEEP_DIVE_PERMUTATIONS",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_PERMUTATIONS)
} else {
  9999L
}

DEEP_DIVE_STANDARDIZE_DISTANCE <- if (exists(
  "DEEP_DIVE_STANDARDIZE_DISTANCE",
  inherits = TRUE
)) {
  isTRUE(DEEP_DIVE_STANDARDIZE_DISTANCE)
} else {
  TRUE
}

DEEP_DIVE_TOP_LOADINGS <- if (exists(
  "DEEP_DIVE_TOP_LOADINGS",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_TOP_LOADINGS)
} else {
  10L
}

if (!is.finite(DEEP_DIVE_PERMUTATIONS) || DEEP_DIVE_PERMUTATIONS < 99L) {
  stop("DEEP_DIVE_PERMUTATIONS must be at least 99.", call. = FALSE)
}

if (!is.finite(DEEP_DIVE_TOP_LOADINGS) || DEEP_DIVE_TOP_LOADINGS < 1L) {
  stop("DEEP_DIVE_TOP_LOADINGS must be >= 1.", call. = FALSE)
}

DEEP_DIVE_PLOT_TOP_N_CLIFFS <- if (exists(
  "DEEP_DIVE_PLOT_TOP_N_CLIFFS",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_PLOT_TOP_N_CLIFFS)
} else {
  20L
}

DEEP_DIVE_PLOT_TOP_N_MAD <- if (exists(
  "DEEP_DIVE_PLOT_TOP_N_MAD",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_PLOT_TOP_N_MAD)
} else {
  20L
}

if (!is.finite(DEEP_DIVE_PLOT_TOP_N_CLIFFS) || DEEP_DIVE_PLOT_TOP_N_CLIFFS < 1L) {
  stop("DEEP_DIVE_PLOT_TOP_N_CLIFFS must be >= 1.", call. = FALSE)
}

if (!is.finite(DEEP_DIVE_PLOT_TOP_N_MAD) || DEEP_DIVE_PLOT_TOP_N_MAD < 1L) {
  stop("DEEP_DIVE_PLOT_TOP_N_MAD must be >= 1.", call. = FALSE)
}

# =============================================================================
# 2. OUTPUT DIRECTORIES
# =============================================================================

DEEP_DIVE_RESULTS_DIR <- file.path(DIRS$results, "04_deep_dive")
DEEP_DIVE_FIGURE_DIR <- file.path(DIRS$figures, "04_deep_dive")
ensure_dir(DEEP_DIVE_RESULTS_DIR)
ensure_dir(DEEP_DIVE_FIGURE_DIR)

# =============================================================================
# 3. LOAD FROZEN INPUTS
# =============================================================================

stage1 <- read_rds_checked(
  RDS_PATHS$stage1,
  "Stage 1 omics screen RDS"
)

biomarker <- read_rds_checked(
  RDS_PATHS$biomarker,
  "biomarker RDS"
)

assignments <- biomarker$clustering_assignments
clinical <- biomarker$clinical

if (!is.data.frame(assignments)) {
  stop(
    "03_biomarker.rds does not contain clustering_assignments.",
    call. = FALSE
  )
}

if (!is.data.frame(clinical)) {
  stop(
    "03_biomarker.rds does not contain clinical covariates.",
    call. = FALSE
  )
}

assert_columns(
  assignments,
  c("study_id", "fcm_profile"),
  "clustering assignments"
)

assert_columns(
  clinical,
  c("study_id", "age", "bmi"),
  "clinical covariates"
)

assignments <- assignments %>%
  dplyr::mutate(
    study_id = normalise_id(study_id),
    fcm_profile = as.character(fcm_profile)
  )

clinical <- clinical %>%
  dplyr::mutate(study_id = normalise_id(study_id))

profile_levels <- sort(unique(assignments$fcm_profile))

# =============================================================================
# 4. PROFILE LABELS
# =============================================================================

DIVERSE_PROFILE <- "Food-diverse"
LIMITED_PROFILE <- "Food-limited"

required_profiles <- c(LIMITED_PROFILE, DIVERSE_PROFILE)
missing_profiles <- setdiff(required_profiles, profile_levels)

if (length(missing_profiles) > 0L) {
  stop(
    "Aim 2 profile-level analyses require profile labels 'Food-limited' and ",
    "'Food-diverse'. Missing: ",
    paste(missing_profiles, collapse = ", "),
    call. = FALSE
  )
}

# =============================================================================
# 5. SHARED HELPERS
# =============================================================================

safe_plot_stub <- function(x) {
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- gsub("^_+|_+$", "", x)
  x
}

save_plot_pair <- function(
    plot_object,
    stub,
    width = 9,
    height = 6,
    dpi = 300
) {
  png_path <- file.path(
    DEEP_DIVE_FIGURE_DIR,
    paste0(stub, ".png")
  )
  
  pdf_path <- file.path(
    DEEP_DIVE_FIGURE_DIR,
    paste0(stub, ".pdf")
  )
  
  ggplot2::ggsave(
    filename = png_path,
    plot = plot_object,
    width = width,
    height = height,
    dpi = dpi,
    bg = "white"
  )
  
  ggplot2::ggsave(
    filename = pdf_path,
    plot = plot_object,
    width = width,
    height = height,
    bg = "white"
  )
  
  list(png = png_path, pdf = pdf_path)
}

get_panel_data <- function(
    panel_name,
    condition = c("A_excluded", "B_full", "C_full_cohort"),
    profile = NULL
) {
  condition <- match.arg(condition)
  
  panel <- stage1$panel_results[[panel_name]]
  
  if (is.null(panel)) {
    stop("Unknown Stage 1 panel: ", panel_name, call. = FALSE)
  }
  
  feature_cols <- panel$retained_features
  
  if (condition == "C_full_cohort") {
    dat <- panel$processed_cohort_b %>%
      dplyr::mutate(study_id = normalise_id(study_id))
    
    return(
      list(
        data = dat,
        features = feature_cols,
        label = panel$label %||% panel_name,
        flags_removed = character(0)
      )
    )
  }
  
  dat <- panel$processed_cohort_a %>%
    dplyr::mutate(study_id = normalise_id(study_id)) %>%
    dplyr::left_join(
      assignments %>% dplyr::select(study_id, fcm_profile),
      by = "study_id"
    )
  
  if (anyNA(dat$fcm_profile)) {
    missing_ids <- dat$study_id[is.na(dat$fcm_profile)]
    stop(
      panel_name,
      ": Cohort A participant(s) lack profile assignments: ",
      paste(head(missing_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  
  flags <- normalise_id(panel$flagged_ids %||% character(0))
  
  if (condition == "A_excluded") {
    dat <- dat %>%
      dplyr::filter(!study_id %in% flags)
  }
  
  if (!is.null(profile)) {
    dat <- dat %>%
      dplyr::filter(fcm_profile == profile)
  }
  
  list(
    data = dat,
    features = feature_cols,
    label = panel$label %||% panel_name,
    flags_removed = if (condition == "A_excluded") flags else character(0)
  )
}

finite_numeric_matrix <- function(data, feature_cols, label) {
  assert_columns(data, feature_cols, label)
  
  X <- as.matrix(data[, feature_cols, drop = FALSE])
  storage.mode(X) <- "double"
  
  if (any(!is.finite(X))) {
    stop(label, ": feature matrix contains NA/Inf.", call. = FALSE)
  }
  
  X
}

standardize_for_distance <- function(X) {
  if (!DEEP_DIVE_STANDARDIZE_DISTANCE) {
    return(X)
  }
  
  sds <- apply(X, 2, stats::sd)
  keep <- is.finite(sds) & sds > sqrt(.Machine$double.eps)
  
  if (sum(keep) < 2L) {
    stop(
      "Fewer than two non-zero-variance features remain for distance analysis.",
      call. = FALSE
    )
  }
  
  X <- X[, keep, drop = FALSE]
  scale(X, center = TRUE, scale = TRUE)
}

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  if (length(x) == 0L || length(y) == 0L) {
    return(NA_real_)
  }
  
  comparisons <- outer(
    x,
    y,
    FUN = function(a, b) sign(a - b)
  )
  
  sum(comparisons) / (length(x) * length(y))
}

cliffs_delta_magnitude <- function(delta) {
  abs_delta <- abs(delta)
  
  dplyr::case_when(
    is.na(abs_delta) ~ NA_character_,
    abs_delta < 0.147 ~ "negligible",
    abs_delta < 0.330 ~ "small",
    abs_delta < 0.474 ~ "medium",
    TRUE ~ "large"
  )
}

run_spearman <- function(x, y) {
  complete <- is.finite(x) & is.finite(y)
  
  x <- x[complete]
  y <- y[complete]
  
  if (
    length(x) < 5L ||
    stats::sd(x) <= sqrt(.Machine$double.eps) ||
    stats::sd(y) <= sqrt(.Machine$double.eps)
  ) {
    return(
      list(
        n = length(x),
        rho = NA_real_,
        p = NA_real_
      )
    )
  }
  
  test <- tryCatch(
    suppressWarnings(
      stats::cor.test(
        x,
        y,
        method = "spearman",
        exact = FALSE
      )
    ),
    error = function(e) NULL
  )
  
  list(
    n = length(x),
    rho = if (!is.null(test)) {
      unname(test$estimate)
    } else {
      NA_real_
    },
    p = if (!is.null(test)) {
      test$p.value
    } else {
      NA_real_
    }
  )
}

format_fdr_star <- function(fdr) {
  dplyr::case_when(
    !is.na(fdr) & fdr < 0.05 ~ "**",
    !is.na(fdr) & fdr < FDR_THRESHOLD ~ "*",
    TRUE ~ ""
  )
}

pretty_nutrient_label <- function(x) {
  label_map <- c(
    retinol_eq = "Retinol equivalents",
    thiamin = "Thiamin",
    riboflavin = "Riboflavin",
    niacin_eq = "Niacin equivalents",
    total_folates = "Dietary folate equivalents",
    vit_b6 = "Vitamin B6",
    vit_b12 = "Vitamin B12",
    vit_c = "Vitamin C",
    vit_e = "Vitamin E",
    calcium = "Calcium",
    iodine = "Iodine",
    iron = "Iron",
    magnesium = "Magnesium",
    phosphorus = "Phosphorus",
    potassium = "Potassium",
    selenium = "Selenium",
    sodium = "Sodium",
    zinc = "Zinc"
  )
  
  out <- unname(label_map[x])
  missing <- is.na(out)
  out[missing] <- x[missing] %>%
    stringr::str_replace_all("_", " ") %>%
    stringr::str_squish() %>%
    stringr::str_to_sentence()
  out
}

get_cytokine_z_data <- function(
    condition = c("B_full", "A_excluded")
) {
  condition <- match.arg(condition)
  
  obj <- get_panel_data(
    panel_name = "cytokines",
    condition = condition
  )
  
  dat <- obj$data
  feature_cols <- obj$features
  
  X <- finite_numeric_matrix(
    dat,
    feature_cols,
    paste0("cytokines ", condition)
  )
  
  sds <- apply(X, 2, stats::sd)
  keep <- is.finite(sds) & sds > sqrt(.Machine$double.eps)
  
  X_z <- scale(
    X[, keep, drop = FALSE],
    center = TRUE,
    scale = TRUE
  )
  
  z_data <- tibble::as_tibble(
    X_z,
    .name_repair = "minimal"
  ) %>%
    dplyr::mutate(
      study_id = dat$study_id,
      case_status = dat$case_status,
      fcm_profile = dat$fcm_profile,
      .before = 1
    )
  
  list(
    data = z_data,
    features = colnames(X_z),
    label = obj$label
  )
}

cat("\n============================================================\n")
cat("04 DEEP DIVE — COMMON SETUP LOADED\n")
cat("============================================================\n")
cat("Profiles:", paste(profile_levels, collapse = ", "), "\n")
cat("PERMANOVA permutations:", DEEP_DIVE_PERMUTATIONS, "\n")
cat(
  "Standardise features for distance analyses:",
  DEEP_DIVE_STANDARDIZE_DISTANCE,
  "\n"
)

# =============================================================================
# End 04_deep_dive_common.R
# =============================================================================