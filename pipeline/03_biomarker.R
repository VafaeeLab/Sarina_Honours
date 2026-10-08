# =============================================================================
# 03_biomarker.R
# Frozen-omics biomarker analysis within diet-defined FCM profiles
#
# PURPOSE
#   Run the primary/sensitivity/reference biomarker analyses using ONLY frozen
#   outputs from earlier pipeline stages. This script does not re-load or
#   re-transform raw omics data.
#
# INPUTS
#   RDS_PATHS$stage1      -- frozen omics preprocessing + T2 flags
#   RDS_PATHS$clustering  -- frozen FCM participant assignments
#   Clinical files defined in CLINICAL_SPEC
#
# CONDITIONS
#   A_excluded      -- within dietary profile; panel-specific Stage 1 T2 flags
#                      excluded (primary)
#   B_full          -- within dietary profile; all Cohort A participants retained
#                      (outlier sensitivity)
#   C_full_cohort   -- Cohort B, unstratified reference analysis
#
# ANALYSES
#   1. Wilcoxon rank-sum: HCC vs Cirrhosis, feature by feature
#   2. Benjamini-Hochberg FDR within each panel x profile x condition analysis
#   3. Cliff's delta + magnitude + median difference/effect direction
#   4. For features with FDR < FDR_THRESHOLD:
#        HCC ~ feature_z + age_z
#        HCC ~ feature_z + bmi_z
#
# OUTPUTS
#   RDS_PATHS$biomarker
#   results/03_biomarker/*.csv
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
})

source("pipeline/config.R")
source("pipeline/utils.R")

set.seed(MASTER_SEED)

BIOMARKER_OUTPUT_DIR <- file.path(DIRS$results, "03_biomarker")
ensure_dir(BIOMARKER_OUTPUT_DIR)
ensure_parent_dir(RDS_PATHS$biomarker)

# =============================================================================
# 1. LOAD FROZEN UPSTREAM OBJECTS
# =============================================================================

stage1 <- read_rds_checked(
  RDS_PATHS$stage1,
  "Stage 1 omics screen RDS"
)

clustering <- read_rds_checked(
  RDS_PATHS$clustering,
  "clustering RDS"
)

if (is.null(stage1$panel_results) || length(stage1$panel_results) == 0L) {
  stop("Stage 1 RDS does not contain panel_results.", call. = FALSE)
}

# -----------------------------------------------------------------------------
# Extract frozen clustering assignments.
#
# Canonical 01_clustering.rds contract should contain $assignments, but several
# aliases are accepted so this script can survive the transition from the old
# clustering code.
# -----------------------------------------------------------------------------

extract_clustering_assignments <- function(x) {
  candidates <- list(
    x$assignments,
    x$assignment_output,
    x$fcm_assignments,
    x$fcm$assignments,
    x$output$assignments
  )
  
  candidates <- candidates[
    vapply(candidates, is.data.frame, logical(1))
  ]
  
  if (length(candidates) == 0L) {
    stop(
      "01_clustering.rds must contain a participant assignment table, ",
      "preferably clustering$assignments.",
      call. = FALSE
    )
  }
  
  dat <- candidates[[1]]
  
  id_col <- intersect(
    c("study_id", "id", "participant_id"),
    names(dat)
  )[1]
  
  profile_col <- intersect(
    c("fcm_profile", "profile", "diet_profile"),
    names(dat)
  )[1]
  
  cluster_col <- intersect(
    c("fcm_cluster", "cluster", "diet_cluster"),
    names(dat)
  )[1]
  
  if (is.na(id_col)) {
    stop("Clustering assignments contain no participant-ID column.", call. = FALSE)
  }
  
  if (is.na(profile_col)) {
    stop(
      "Clustering assignments contain no profile label column. Expected ",
      "'fcm_profile' (preferred), 'profile', or 'diet_profile'.",
      call. = FALSE
    )
  }
  
  out <- tibble::tibble(
    study_id = normalise_id(dat[[id_col]]),
    fcm_profile = as.character(dat[[profile_col]])
  )
  
  if (!is.na(cluster_col)) {
    out$fcm_cluster <- dat[[cluster_col]]
  } else {
    out$fcm_cluster <- NA_integer_
  }
  
  membership_col <- intersect(
    c("fcm_max_membership", "max_membership"),
    names(dat)
  )[1]
  
  if (!is.na(membership_col)) {
    out$fcm_max_membership <- suppressWarnings(
      as.numeric(dat[[membership_col]])
    )
  } else {
    out$fcm_max_membership <- NA_real_
  }
  
  if (anyNA(out$study_id)) {
    stop("Clustering assignments contain missing participant IDs.", call. = FALSE)
  }
  
  if (anyDuplicated(out$study_id)) {
    stop("Clustering assignments contain duplicated participant IDs.", call. = FALSE)
  }
  
  if (anyNA(out$fcm_profile) || any(out$fcm_profile == "")) {
    stop("Clustering assignments contain missing profile labels.", call. = FALSE)
  }
  
  out
}

assignments <- extract_clustering_assignments(clustering)
profile_levels <- sort(unique(assignments$fcm_profile))

if (length(profile_levels) < 2L) {
  stop(
    "At least two dietary profiles are required for the within-profile analysis.",
    call. = FALSE
  )
}

cat("\n============================================================\n")
cat("03 BIOMARKER ANALYSIS\n")
cat("============================================================\n")
cat("Dietary profiles:", paste(profile_levels, collapse = ", "), "\n")
cat("FDR threshold:", FDR_THRESHOLD, "\n")

# =============================================================================
# 2. LOAD CLINICAL COVARIATES
# =============================================================================

load_clinical_covariates <- function() {
  assert_file_exists(
    CLINICAL_SPEC$cirrhosis_file,
    "Cirrhosis clinical file"
  )
  assert_file_exists(
    CLINICAL_SPEC$hcc_file,
    "HCC clinical file"
  )
  
  read_one <- function(path, case_status) {
    dat <- readr::read_csv(
      path,
      show_col_types = FALSE,
      name_repair = "minimal"
    )
    
    assert_columns(
      dat,
      c(
        CLINICAL_SPEC$id_col,
        CLINICAL_SPEC$age_col,
        CLINICAL_SPEC$bmi_col
      ),
      paste(case_status, "clinical data")
    )
    
    tibble::tibble(
      study_id = normalise_id(dat[[CLINICAL_SPEC$id_col]]),
      case_status = case_status,
      age = suppressWarnings(as.numeric(dat[[CLINICAL_SPEC$age_col]])),
      bmi = suppressWarnings(as.numeric(dat[[CLINICAL_SPEC$bmi_col]]))
    )
  }
  
  clinical <- dplyr::bind_rows(
    read_one(CLINICAL_SPEC$cirrhosis_file, "Cirrhosis"),
    read_one(CLINICAL_SPEC$hcc_file, "HCC")
  )
  
  if (anyNA(clinical$study_id)) {
    stop("Clinical data contain missing participant IDs.", call. = FALSE)
  }
  
  duplicate_ids <- clinical$study_id[duplicated(clinical$study_id)] %>% unique()
  if (length(duplicate_ids) > 0L) {
    stop(
      "Clinical data contain duplicated participant IDs: ",
      paste(head(duplicate_ids, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  
  clinical
}

clinical <- load_clinical_covariates()

cat("Clinical participants:", nrow(clinical), "\n")
cat("Age missing:", sum(is.na(clinical$age)), "\n")
cat("BMI missing:", sum(is.na(clinical$bmi)), "\n")

# =============================================================================
# 3. ANALYSIS HELPERS
# =============================================================================

cliffs_delta <- function(x, y) {
  x <- x[is.finite(x)]
  y <- y[is.finite(y)]
  
  nx <- length(x)
  ny <- length(y)
  
  if (nx == 0L || ny == 0L) {
    return(NA_real_)
  }
  
  comparisons <- outer(
    x,
    y,
    FUN = function(a, b) sign(a - b)
  )
  
  sum(comparisons) / (nx * ny)
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

run_univariate <- function(data, feature_cols, analysis_label) {
  assert_columns(
    data,
    c("study_id", "case_status", feature_cols),
    analysis_label
  )
  
  case_status <- as.character(data$case_status)
  valid_status <- case_status %in% c("HCC", "Cirrhosis")
  
  if (!all(valid_status)) {
    warning(
      analysis_label,
      ": dropping participants with unexpected case_status values.",
      call. = FALSE
    )
    data <- data[valid_status, , drop = FALSE]
  }
  
  n_hcc <- sum(data$case_status == "HCC")
  n_cirr <- sum(data$case_status == "Cirrhosis")
  
  if (n_hcc < 3L || n_cirr < 3L) {
    warning(
      analysis_label,
      ": too few participants for Wilcoxon testing (HCC=", n_hcc,
      ", Cirrhosis=", n_cirr, ").",
      call. = FALSE
    )
    return(NULL)
  }
  
  results <- purrr::map_dfr(
    feature_cols,
    function(feature) {
      x_hcc <- suppressWarnings(
        as.numeric(data[[feature]][data$case_status == "HCC"])
      )
      x_cirr <- suppressWarnings(
        as.numeric(data[[feature]][data$case_status == "Cirrhosis"])
      )
      
      x_hcc <- x_hcc[is.finite(x_hcc)]
      x_cirr <- x_cirr[is.finite(x_cirr)]
      
      if (length(x_hcc) < 3L || length(x_cirr) < 3L) {
        return(
          tibble::tibble(
            feature = feature,
            n_hcc = length(x_hcc),
            n_cirrhosis = length(x_cirr),
            median_hcc = if (length(x_hcc) > 0L) median(x_hcc) else NA_real_,
            median_cirrhosis = if (length(x_cirr) > 0L) median(x_cirr) else NA_real_,
            median_diff = NA_real_,
            cliffs_delta = NA_real_,
            wilcoxon_p = NA_real_
          )
        )
      }
      
      wt <- tryCatch(
        stats::wilcox.test(
          x_hcc,
          x_cirr,
          exact = FALSE,
          correct = TRUE
        ),
        error = function(e) NULL
      )
      
      med_hcc <- stats::median(x_hcc)
      med_cirr <- stats::median(x_cirr)
      
      tibble::tibble(
        feature = feature,
        n_hcc = length(x_hcc),
        n_cirrhosis = length(x_cirr),
        median_hcc = med_hcc,
        median_cirrhosis = med_cirr,
        median_diff = med_hcc - med_cirr,
        cliffs_delta = cliffs_delta(x_hcc, x_cirr),
        wilcoxon_p = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }
  )
  
  results %>%
    dplyr::mutate(
      fdr_p = stats::p.adjust(wilcoxon_p, method = "BH"),
      direction = dplyr::case_when(
        median_diff > 0 ~ "HCC > Cirrhosis",
        median_diff < 0 ~ "HCC < Cirrhosis",
        is.na(median_diff) ~ NA_character_,
        TRUE ~ "No median difference"
      ),
      cliffs_magnitude = cliffs_delta_magnitude(cliffs_delta),
      fdr_significant = !is.na(fdr_p) & fdr_p < FDR_THRESHOLD
    ) %>%
    dplyr::arrange(wilcoxon_p)
}

fit_adjusted_logistic <- function(model_df, covariate_name) {
  model_df <- model_df %>%
    dplyr::filter(
      is.finite(feature_val),
      is.finite(covar_val),
      !is.na(hcc_binary)
    )
  
  n_complete <- nrow(model_df)
  n_hcc <- sum(model_df$hcc_binary == 1L)
  n_cirr <- sum(model_df$hcc_binary == 0L)
  
  empty_result <- function() {
    tibble::tibble(
      covariate = covariate_name,
      n_complete = n_complete,
      n_hcc = n_hcc,
      n_cirrhosis = n_cirr,
      feature_or = NA_real_,
      feature_ci_low = NA_real_,
      feature_ci_high = NA_real_,
      feature_p = NA_real_,
      covariate_or = NA_real_,
      covariate_p = NA_real_,
      model_converged = FALSE
    )
  }
  
  if (
    n_complete < 10L ||
    n_hcc < 3L ||
    n_cirr < 3L ||
    stats::sd(model_df$feature_val) <= sqrt(.Machine$double.eps) ||
    stats::sd(model_df$covar_val) <= sqrt(.Machine$double.eps)
  ) {
    return(empty_result())
  }
  
  model_df <- model_df %>%
    dplyr::mutate(
      feature_z = as.numeric(scale(feature_val)),
      covar_z = as.numeric(scale(covar_val))
    )
  
  fit <- tryCatch(
    suppressWarnings(
      stats::glm(
        hcc_binary ~ feature_z + covar_z,
        data = model_df,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(empty_result())
  }
  
  coef_table <- summary(fit)$coefficients
  
  if (!all(c("feature_z", "covar_z") %in% rownames(coef_table))) {
    return(empty_result())
  }
  
  feature_est <- coef_table["feature_z", "Estimate"]
  feature_se <- coef_table["feature_z", "Std. Error"]
  
  tibble::tibble(
    covariate = covariate_name,
    n_complete = n_complete,
    n_hcc = n_hcc,
    n_cirrhosis = n_cirr,
    feature_or = exp(feature_est),
    feature_ci_low = exp(feature_est - 1.96 * feature_se),
    feature_ci_high = exp(feature_est + 1.96 * feature_se),
    feature_p = coef_table["feature_z", "Pr(>|z|)"],
    covariate_or = exp(coef_table["covar_z", "Estimate"]),
    covariate_p = coef_table["covar_z", "Pr(>|z|)"],
    model_converged = isTRUE(fit$converged)
  )
}

run_confounder_adjustment <- function(
    data,
    feature_cols,
    clinical_data,
    univariate_results,
    fdr_cutoff
) {
  if (is.null(univariate_results) || nrow(univariate_results) == 0L) {
    return(tibble::tibble())
  }
  
  sig_features <- univariate_results %>%
    dplyr::filter(!is.na(fdr_p), fdr_p < fdr_cutoff) %>%
    dplyr::pull(feature) %>%
    unique()
  
  sig_features <- intersect(sig_features, feature_cols)
  
  if (length(sig_features) == 0L) {
    return(tibble::tibble())
  }
  
  analysis_data <- data %>%
    dplyr::select(
      study_id,
      case_status,
      dplyr::all_of(sig_features)
    ) %>%
    dplyr::mutate(
      hcc_binary = as.integer(case_status == "HCC")
    ) %>%
    dplyr::left_join(
      clinical_data %>% dplyr::select(study_id, age, bmi),
      by = "study_id"
    )
  
  purrr::map_dfr(
    sig_features,
    function(feature) {
      purrr::map_dfr(
        c("age", "bmi"),
        function(covariate) {
          model_df <- tibble::tibble(
            hcc_binary = analysis_data$hcc_binary,
            feature_val = suppressWarnings(
              as.numeric(analysis_data[[feature]])
            ),
            covar_val = suppressWarnings(
              as.numeric(analysis_data[[covariate]])
            )
          )
          
          fit_adjusted_logistic(
            model_df,
            covariate_name = covariate
          ) %>%
            dplyr::mutate(feature = feature, .before = 1)
        }
      )
    }
  )
}

add_analysis_labels <- function(
    data,
    panel_name,
    panel_label,
    profile,
    condition
) {
  data %>%
    dplyr::mutate(
      panel_name = panel_name,
      panel = panel_label,
      profile = profile,
      condition = condition,
      .before = 1
    )
}

# =============================================================================
# 4. RUN CONDITIONS A, B AND C
# =============================================================================

all_univariate <- list()
all_adjusted <- list()
analysis_counts <- list()

for (panel_name in names(stage1$panel_results)) {
  panel <- stage1$panel_results[[panel_name]]
  panel_label <- panel$label %||% panel_name
  feature_cols <- panel$retained_features
  
  if (
    is.null(panel$processed_cohort_a) ||
    is.null(panel$processed_cohort_b) ||
    length(feature_cols) == 0L
  ) {
    stop(
      "Stage 1 panel '", panel_name,
      "' is missing frozen processed data or retained features.",
      call. = FALSE
    )
  }
  
  panel_a <- panel$processed_cohort_a %>%
    dplyr::mutate(study_id = normalise_id(study_id)) %>%
    dplyr::left_join(
      assignments %>% dplyr::select(study_id, fcm_profile),
      by = "study_id"
    )
  
  unmatched <- panel_a %>%
    dplyr::filter(is.na(fcm_profile)) %>%
    dplyr::pull(study_id)
  
  if (length(unmatched) > 0L) {
    stop(
      panel_label,
      ": Cohort A participant(s) lack clustering assignments: ",
      paste(head(unmatched, 10L), collapse = ", "),
      call. = FALSE
    )
  }
  
  flags <- normalise_id(panel$flagged_ids %||% character(0))
  
  cat("\n------------------------------------------------------------\n")
  cat("PANEL:", panel_label, "\n")
  cat("Retained features:", length(feature_cols), "\n")
  cat("Stage 1 flagged participants:", length(flags), "\n")
  
  # ---------------------------------------------------------------------------
  # Conditions A/B: within each frozen dietary profile in Cohort A
  # ---------------------------------------------------------------------------
  
  for (profile in profile_levels) {
    profile_data <- panel_a %>%
      dplyr::filter(fcm_profile == profile)
    
    condition_sets <- list(
      A_excluded = profile_data %>%
        dplyr::filter(!study_id %in% flags),
      B_full = profile_data
    )
    
    for (condition in names(condition_sets)) {
      dat <- condition_sets[[condition]]
      analysis_label <- paste(panel_name, profile, condition, sep = "__")
      
      cat(
        profile, "|", condition,
        "| n =", nrow(dat),
        "| HCC =", sum(dat$case_status == "HCC"),
        "| Cirrhosis =", sum(dat$case_status == "Cirrhosis"),
        "\n"
      )
      
      analysis_counts[[analysis_label]] <- tibble::tibble(
        panel_name = panel_name,
        panel = panel_label,
        profile = profile,
        condition = condition,
        n_total = nrow(dat),
        n_hcc = sum(dat$case_status == "HCC"),
        n_cirrhosis = sum(dat$case_status == "Cirrhosis"),
        n_stage1_flags_removed = if (condition == "A_excluded") {
          sum(profile_data$study_id %in% flags)
        } else {
          0L
        }
      )
      
      uni <- run_univariate(
        dat,
        feature_cols,
        analysis_label
      )
      
      if (!is.null(uni)) {
        uni <- add_analysis_labels(
          uni,
          panel_name,
          panel_label,
          profile,
          condition
        )
        
        all_univariate[[analysis_label]] <- uni
        
        adj <- run_confounder_adjustment(
          data = dat,
          feature_cols = feature_cols,
          clinical_data = clinical,
          univariate_results = uni,
          fdr_cutoff = FDR_THRESHOLD
        )
        
        if (nrow(adj) > 0L) {
          all_adjusted[[analysis_label]] <- add_analysis_labels(
            adj,
            panel_name,
            panel_label,
            profile,
            condition
          )
        }
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Condition C: frozen Cohort B, unstratified reference
  # ---------------------------------------------------------------------------
  
  panel_c <- panel$processed_cohort_b %>%
    dplyr::mutate(study_id = normalise_id(study_id))
  
  condition <- "C_full_cohort"
  profile <- "All (unstratified)"
  analysis_label <- paste(panel_name, "all", condition, sep = "__")
  
  cat(
    profile, "|", condition,
    "| n =", nrow(panel_c),
    "| HCC =", sum(panel_c$case_status == "HCC"),
    "| Cirrhosis =", sum(panel_c$case_status == "Cirrhosis"),
    "\n"
  )
  
  analysis_counts[[analysis_label]] <- tibble::tibble(
    panel_name = panel_name,
    panel = panel_label,
    profile = profile,
    condition = condition,
    n_total = nrow(panel_c),
    n_hcc = sum(panel_c$case_status == "HCC"),
    n_cirrhosis = sum(panel_c$case_status == "Cirrhosis"),
    n_stage1_flags_removed = 0L
  )
  
  uni_c <- run_univariate(
    panel_c,
    feature_cols,
    analysis_label
  )
  
  if (!is.null(uni_c)) {
    uni_c <- add_analysis_labels(
      uni_c,
      panel_name,
      panel_label,
      profile,
      condition
    )
    
    all_univariate[[analysis_label]] <- uni_c
    
    adj_c <- run_confounder_adjustment(
      data = panel_c,
      feature_cols = feature_cols,
      clinical_data = clinical,
      univariate_results = uni_c,
      fdr_cutoff = FDR_THRESHOLD
    )
    
    if (nrow(adj_c) > 0L) {
      all_adjusted[[analysis_label]] <- add_analysis_labels(
        adj_c,
        panel_name,
        panel_label,
        profile,
        condition
      )
    }
  }
}

# =============================================================================
# 5. COMBINE AND SUMMARISE
# =============================================================================

univariate_all <- dplyr::bind_rows(all_univariate)
adjusted_all <- dplyr::bind_rows(all_adjusted)
analysis_counts <- dplyr::bind_rows(analysis_counts)

if (nrow(univariate_all) == 0L) {
  stop("No univariate biomarker analyses completed successfully.", call. = FALSE)
}

summary_table <- univariate_all %>%
  dplyr::group_by(
    panel_name,
    panel,
    profile,
    condition
  ) %>%
  dplyr::summarise(
    n_hcc = dplyr::first(n_hcc),
    n_cirrhosis = dplyr::first(n_cirrhosis),
    n_features = dplyr::n(),
    n_nominal_p05 = sum(wilcoxon_p < 0.05, na.rm = TRUE),
    n_fdr_p10 = sum(fdr_p < 0.10, na.rm = TRUE),
    n_fdr_p05 = sum(fdr_p < 0.05, na.rm = TRUE),
    .groups = "drop"
  )

significant_features <- univariate_all %>%
  dplyr::filter(!is.na(fdr_p), fdr_p < FDR_THRESHOLD) %>%
  dplyr::arrange(panel_name, feature, profile, condition, fdr_p)

# A/B concordance by dietary profile, with Condition C attached as an
# unstratified reference for the same panel x feature.
ab_concordance <- univariate_all %>%
  dplyr::filter(condition %in% c("A_excluded", "B_full")) %>%
  dplyr::mutate(significant = fdr_p < FDR_THRESHOLD) %>%
  dplyr::select(
    panel_name,
    panel,
    profile,
    feature,
    condition,
    significant,
    fdr_p
  ) %>%
  tidyr::pivot_wider(
    names_from = condition,
    values_from = c(significant, fdr_p),
    names_sep = "__"
  )

c_reference <- univariate_all %>%
  dplyr::filter(condition == "C_full_cohort") %>%
  dplyr::transmute(
    panel_name,
    feature,
    significant__C_full_cohort = fdr_p < FDR_THRESHOLD,
    fdr_p__C_full_cohort = fdr_p
  )

feature_concordance <- ab_concordance %>%
  dplyr::left_join(
    c_reference,
    by = c("panel_name", "feature")
  ) %>%
  dplyr::arrange(panel_name, feature, profile)

# =============================================================================
# 6. EXPORT CSV RESULTS
# =============================================================================

readr::write_csv(
  univariate_all,
  file.path(BIOMARKER_OUTPUT_DIR, "univariate_results_all.csv")
)

readr::write_csv(
  analysis_counts,
  file.path(BIOMARKER_OUTPUT_DIR, "analysis_sample_counts.csv")
)

readr::write_csv(
  summary_table,
  file.path(BIOMARKER_OUTPUT_DIR, "stage3_summary.csv")
)

readr::write_csv(
  significant_features,
  file.path(BIOMARKER_OUTPUT_DIR, "significant_features_fdr.csv")
)

readr::write_csv(
  feature_concordance,
  file.path(BIOMARKER_OUTPUT_DIR, "feature_concordance_across_conditions.csv")
)

if (nrow(adjusted_all) > 0L) {
  readr::write_csv(
    adjusted_all,
    file.path(BIOMARKER_OUTPUT_DIR, "confounder_adjusted_all.csv")
  )
}

for (panel_name in unique(univariate_all$panel_name)) {
  panel_results <- univariate_all %>%
    dplyr::filter(.data$panel_name == panel_name)
  
  readr::write_csv(
    panel_results,
    file.path(
      BIOMARKER_OUTPUT_DIR,
      paste0("univariate_results_", panel_name, ".csv")
    )
  )
}

# =============================================================================
# 7. SAVE FROZEN BIOMARKER OBJECT
# =============================================================================

biomarker_output <- list(
  univariate = univariate_all,
  confounder_adjusted = adjusted_all,
  significant_features = significant_features,
  concordance = feature_concordance,
  summary = summary_table,
  analysis_counts = analysis_counts,
  
  # Stored so 04_deep_dive.R does not need to reload these sources.
  clustering_assignments = assignments,
  clinical = clinical,
  
  settings = list(
    seed = MASTER_SEED,
    fdr_threshold = FDR_THRESHOLD,
    conditions = BIOMARKER_CONDITIONS,
    univariate_test = "Wilcoxon rank-sum",
    multiplicity = "Benjamini-Hochberg within each panel/profile/condition analysis",
    effect_size = "Cliff's delta (positive = HCC tends higher than Cirrhosis)",
    cliffs_delta_magnitude_cutoffs = c(
      negligible = 0.147,
      small = 0.330,
      medium = 0.474
    ),
    adjusted_models = c(
      "HCC ~ feature_z + age_z",
      "HCC ~ feature_z + bmi_z"
    ),
    omics_preprocessing_source = RDS_PATHS$stage1,
    clustering_source = RDS_PATHS$clustering
  )
)

save_rds_checked(
  biomarker_output,
  RDS_PATHS$biomarker
)

# =============================================================================
# 8. CONSOLE SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("03 BIOMARKER ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Saved RDS:", RDS_PATHS$biomarker, "\n")
cat("CSV directory:", BIOMARKER_OUTPUT_DIR, "\n\n")

print(summary_table, n = Inf, width = Inf)

cat("\nTotal FDR <", FDR_THRESHOLD, "results:", nrow(significant_features), "\n")
cat("============================================================\n")
