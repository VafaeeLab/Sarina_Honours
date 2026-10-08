# =============================================================================
# 04e_interaction.R
# Exploratory biomarker x dietary-profile interaction models
#
# PURPOSE
#   For features reaching FDR < threshold in at least one within-profile
#   Condition A or B biomarker analysis:
#     HCC ~ feature_z * fcm_profile + age_z + bmi_z
#   BH-FDR is applied across all interaction tests within each panel/condition.
#
# INPUTS
#   RDS_PATHS$stage1
#   RDS_PATHS$biomarker
#
# OUTPUTS
#   intermediate/04e_interaction.rds
#   results/04_deep_dive/interaction_candidate_features.csv
#   results/04_deep_dive/profile_interaction_models.csv
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")
source("pipeline/04_COMMON.R")

set.seed(MASTER_SEED)

# =============================================================================
# 13. FEATURE x DIET-PROFILE INTERACTION ANALYSES
# =============================================================================
# These are explicitly exploratory. Candidate features are restricted to those
# reaching FDR < FDR_THRESHOLD in at least one within-profile Condition A or B
# analysis in 03_biomarker.R.
# =============================================================================

candidate_interaction_features <- biomarker$univariate %>%
  dplyr::filter(
    condition %in% c("A_excluded", "B_full"),
    !is.na(fdr_p),
    fdr_p < FDR_THRESHOLD
  ) %>%
  dplyr::distinct(panel_name, feature) %>%
  dplyr::arrange(panel_name, feature)

fit_profile_interaction <- function(
    panel_name,
    feature,
    condition = c("A_excluded", "B_full")
) {
  condition <- match.arg(condition)
  
  obj <- get_panel_data(
    panel_name = panel_name,
    condition = condition
  )
  
  if (!feature %in% obj$features) {
    return(tibble::tibble())
  }
  
  dat <- obj$data %>%
    dplyr::select(
      study_id,
      case_status,
      fcm_profile,
      feature_val = dplyr::all_of(feature)
    ) %>%
    dplyr::left_join(
      clinical %>% dplyr::select(study_id, age, bmi),
      by = "study_id"
    ) %>%
    dplyr::filter(
      case_status %in% c("Cirrhosis", "HCC"),
      !is.na(fcm_profile),
      is.finite(feature_val),
      is.finite(age),
      is.finite(bmi)
    )
  
  empty <- function() {
    tibble::tibble(
      panel_name = panel_name,
      panel = obj$label,
      condition = condition,
      feature = feature,
      interaction_term = NA_character_,
      n_complete = nrow(dat),
      interaction_log_or_ratio = NA_real_,
      interaction_or_ratio = NA_real_,
      interaction_ci_low = NA_real_,
      interaction_ci_high = NA_real_,
      interaction_p = NA_real_,
      model_converged = FALSE
    )
  }
  
  if (
    nrow(dat) < 20L ||
    length(unique(dat$fcm_profile)) < 2L ||
    sum(dat$case_status == "HCC") < 5L ||
    sum(dat$case_status == "Cirrhosis") < 5L ||
    stats::sd(dat$feature_val) <= sqrt(.Machine$double.eps) ||
    stats::sd(dat$age) <= sqrt(.Machine$double.eps) ||
    stats::sd(dat$bmi) <= sqrt(.Machine$double.eps)
  ) {
    return(empty())
  }
  
  dat <- dat %>%
    dplyr::mutate(
      hcc_binary = as.integer(case_status == "HCC"),
      fcm_profile = factor(fcm_profile, levels = profile_levels),
      feature_z = as.numeric(scale(feature_val)),
      age_z = as.numeric(scale(age)),
      bmi_z = as.numeric(scale(bmi))
    )
  
  fit <- tryCatch(
    suppressWarnings(
      stats::glm(
        hcc_binary ~ feature_z * fcm_profile + age_z + bmi_z,
        data = dat,
        family = stats::binomial(link = "logit")
      )
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(empty())
  }
  
  coef_table <- summary(fit)$coefficients
  interaction_rows <- grep(
    "feature_z:fcm_profile|fcm_profile.*:feature_z",
    rownames(coef_table),
    value = TRUE
  )
  
  if (length(interaction_rows) == 0L) {
    return(empty())
  }
  
  purrr::map_dfr(
    interaction_rows,
    function(term) {
      estimate <- coef_table[term, "Estimate"]
      se <- coef_table[term, "Std. Error"]
      
      tibble::tibble(
        panel_name = panel_name,
        panel = obj$label,
        condition = condition,
        feature = feature,
        interaction_term = term,
        n_complete = nrow(dat),
        interaction_log_or_ratio = estimate,
        interaction_or_ratio = exp(estimate),
        interaction_ci_low = exp(estimate - 1.96 * se),
        interaction_ci_high = exp(estimate + 1.96 * se),
        interaction_p = coef_table[term, "Pr(>|z|)"],
        model_converged = isTRUE(fit$converged)
      )
    }
  )
}

interaction_results <- list()

if (nrow(candidate_interaction_features) > 0L) {
  interaction_results <- purrr::pmap_dfr(
    candidate_interaction_features,
    function(panel_name, feature) {
      dplyr::bind_rows(
        fit_profile_interaction(
          panel_name,
          feature,
          condition = "A_excluded"
        ),
        fit_profile_interaction(
          panel_name,
          feature,
          condition = "B_full"
        )
      )
    }
  )
  
  if (nrow(interaction_results) > 0L) {
    interaction_results <- interaction_results %>%
      dplyr::group_by(panel_name, condition) %>%
      dplyr::mutate(
        interaction_fdr = stats::p.adjust(
          interaction_p,
          method = "BH"
        )
      ) %>%
      dplyr::ungroup() %>%
      dplyr::arrange(
        panel_name,
        condition,
        interaction_p
      )
  }
} else {
  interaction_results <- tibble::tibble()
}

readr::write_csv(
  candidate_interaction_features,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "interaction_candidate_features.csv"
  )
)

if (nrow(interaction_results) > 0L) {
  readr::write_csv(
    interaction_results,
    file.path(
      DEEP_DIVE_RESULTS_DIR,
      "profile_interaction_models.csv"
    )
  )
}

# =============================================================================
# SAVE
# =============================================================================

interaction_output <- list(
  interaction_candidates = candidate_interaction_features,
  profile_interactions = interaction_results
)

save_rds_checked(
  interaction_output,
  file.path(DIRS$intermediate, "04e_interaction.rds")
)

cat("\n============================================================\n")
cat("04e INTERACTION MODELS COMPLETE\n")
cat("============================================================\n")
cat(
  "Interaction candidate features:",
  nrow(candidate_interaction_features),
  "\n"
)
cat("Interaction result rows:", nrow(interaction_results), "\n")
cat("============================================================\n")