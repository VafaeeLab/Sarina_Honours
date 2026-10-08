# =============================================================================
# 04b_permanova.R
# Multivariate HCC-vs-Cirrhosis and Food-diverse-vs-Food-limited tests
#
# PURPOSE
#   1. HCC-vs-Cirrhosis PERMANOVA/PERMDISP within dietary profiles + reference,
#      with age/BMI-adjusted PERMANOVA sensitivity analysis.
#   2. Food-diverse vs Food-limited PERMANOVA/PERMDISP across each omics panel.
#
# INPUTS
#   RDS_PATHS$stage1
#   RDS_PATHS$biomarker
#
# OUTPUTS
#   intermediate/04b_permanova.rds
#   results/04_deep_dive/permanova_permdisp_results.csv
#   results/04_deep_dive/aim2_diet_profile_permanova_permdisp.csv
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")
source("code/04_COMMON.R")

set.seed(MASTER_SEED)

# =============================================================================
# 1. HCC-VS-CIRRHOSIS PERMANOVA + PERMDISP
# =============================================================================

run_multivariate_test <- function(
    panel_name,
    condition,
    profile = NULL,
    seed_offset = 0L
) {
  obj <- get_panel_data(
    panel_name = panel_name,
    condition = condition,
    profile = profile
  )
  
  dat <- obj$data
  feature_cols <- obj$features
  
  profile_label <- if (is.null(profile)) {
    "All (unstratified)"
  } else {
    profile
  }
  
  n_hcc <- sum(dat$case_status == "HCC")
  n_cirr <- sum(dat$case_status == "Cirrhosis")
  
  base_result <- tibble::tibble(
    panel_name = panel_name,
    panel = obj$label,
    profile = profile_label,
    condition = condition,
    n_total = nrow(dat),
    n_hcc = n_hcc,
    n_cirrhosis = n_cirr,
    n_features = length(feature_cols),
    standardized_for_distance = DEEP_DIVE_STANDARDIZE_DISTANCE,
    permanova_F = NA_real_,
    permanova_R2 = NA_real_,
    permanova_p = NA_real_,
    permdisp_F = NA_real_,
    permdisp_p = NA_real_,
    adjusted_n = NA_integer_,
    adjusted_permanova_F = NA_real_,
    adjusted_permanova_R2 = NA_real_,
    adjusted_permanova_p = NA_real_
  )
  
  if (nrow(dat) < 6L || n_hcc < 3L || n_cirr < 3L) {
    warning(
      panel_name, " / ", profile_label, " / ", condition,
      ": too few participants for multivariate testing.",
      call. = FALSE
    )
    return(base_result)
  }
  
  X <- finite_numeric_matrix(
    dat,
    feature_cols,
    paste(panel_name, profile_label, condition)
  )
  
  X_distance <- standardize_for_distance(X)
  distance <- stats::dist(X_distance, method = "euclidean")
  
  meta <- dat %>%
    dplyr::transmute(
      study_id,
      case_status = factor(
        case_status,
        levels = c("Cirrhosis", "HCC")
      )
    )
  
  set.seed(MASTER_SEED + seed_offset)
  permanova <- tryCatch(
    vegan::adonis2(
      distance ~ case_status,
      data = meta,
      permutations = DEEP_DIVE_PERMUTATIONS,
      by = "margin"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(permanova) && "case_status" %in% rownames(permanova)) {
    base_result$permanova_F <- permanova["case_status", "F"]
    base_result$permanova_R2 <- permanova["case_status", "R2"]
    base_result$permanova_p <- permanova["case_status", "Pr(>F)"]
  }
  
  dispersion <- tryCatch(
    vegan::betadisper(
      distance,
      group = meta$case_status,
      type = "median"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(dispersion)) {
    set.seed(MASTER_SEED + seed_offset + 50000L)
    dispersion_test <- tryCatch(
      vegan::permutest(
        dispersion,
        permutations = DEEP_DIVE_PERMUTATIONS
      ),
      error = function(e) NULL
    )
    
    if (!is.null(dispersion_test)) {
      tab <- dispersion_test$tab
      if (nrow(tab) >= 1L) {
        base_result$permdisp_F <- tab[1, "F"]
        base_result$permdisp_p <- tab[1, "Pr(>F)"]
      }
    }
  }
  
  # ---------------------------------------------------------------------------
  # Adjusted PERMANOVA: disease status + age + BMI, tested marginally.
  # ---------------------------------------------------------------------------
  
  adjusted_meta <- dat %>%
    dplyr::select(study_id, case_status) %>%
    dplyr::left_join(
      clinical %>% dplyr::select(study_id, age, bmi),
      by = "study_id"
    ) %>%
    dplyr::mutate(
      row_index = dplyr::row_number()
    ) %>%
    dplyr::filter(
      is.finite(age),
      is.finite(bmi)
    )
  
  base_result$adjusted_n <- nrow(adjusted_meta)
  
  if (
    nrow(adjusted_meta) >= 10L &&
    sum(adjusted_meta$case_status == "HCC") >= 3L &&
    sum(adjusted_meta$case_status == "Cirrhosis") >= 3L &&
    stats::sd(adjusted_meta$age) > sqrt(.Machine$double.eps) &&
    stats::sd(adjusted_meta$bmi) > sqrt(.Machine$double.eps)
  ) {
    idx <- adjusted_meta$row_index
    X_adjusted <- X_distance[idx, , drop = FALSE]
    
    adjusted_meta <- adjusted_meta %>%
      dplyr::mutate(
        case_status = factor(
          case_status,
          levels = c("Cirrhosis", "HCC")
        ),
        age_z = as.numeric(scale(age)),
        bmi_z = as.numeric(scale(bmi))
      )
    
    adjusted_distance <- stats::dist(
      X_adjusted,
      method = "euclidean"
    )
    
    set.seed(MASTER_SEED + seed_offset + 100000L)
    adjusted_permanova <- tryCatch(
      vegan::adonis2(
        adjusted_distance ~ case_status + age_z + bmi_z,
        data = adjusted_meta,
        permutations = DEEP_DIVE_PERMUTATIONS,
        by = "margin"
      ),
      error = function(e) NULL
    )
    
    if (
      !is.null(adjusted_permanova) &&
      "case_status" %in% rownames(adjusted_permanova)
    ) {
      base_result$adjusted_permanova_F <- adjusted_permanova[
        "case_status",
        "F"
      ]
      base_result$adjusted_permanova_R2 <- adjusted_permanova[
        "case_status",
        "R2"
      ]
      base_result$adjusted_permanova_p <- adjusted_permanova[
        "case_status",
        "Pr(>F)"
      ]
    }
  }
  
  base_result
}

multivariate_results <- list()
test_index <- 0L

for (panel_name in names(stage1$panel_results)) {
  for (profile in profile_levels) {
    for (condition in c("A_excluded", "B_full")) {
      test_index <- test_index + 1L
      key <- paste(panel_name, profile, condition, sep = "__")
      
      multivariate_results[[key]] <- run_multivariate_test(
        panel_name = panel_name,
        condition = condition,
        profile = profile,
        seed_offset = 1000L * test_index
      )
    }
  }
  
  test_index <- test_index + 1L
  key <- paste(panel_name, "all", "C_full_cohort", sep = "__")
  
  multivariate_results[[key]] <- run_multivariate_test(
    panel_name = panel_name,
    condition = "C_full_cohort",
    profile = NULL,
    seed_offset = 1000L * test_index
  )
}

multivariate_results <- dplyr::bind_rows(multivariate_results) %>%
  dplyr::mutate(
    permanova_fdr = stats::p.adjust(permanova_p, method = "BH"),
    permdisp_fdr = stats::p.adjust(permdisp_p, method = "BH"),
    adjusted_permanova_fdr = stats::p.adjust(
      adjusted_permanova_p,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(panel_name, condition, profile)

readr::write_csv(
  multivariate_results,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "permanova_permdisp_results.csv"
  )
)

# =============================================================================
# 2. DIET-PROFILE PERMANOVA + PERMDISP
# =============================================================================

run_diet_profile_multivariate <- function(
    panel_name,
    condition = c("B_full", "A_excluded"),
    seed_offset = 0L
) {
  condition <- match.arg(condition)
  
  obj <- get_panel_data(
    panel_name = panel_name,
    condition = condition
  )
  
  dat <- obj$data %>%
    dplyr::filter(
      fcm_profile %in% c(LIMITED_PROFILE, DIVERSE_PROFILE)
    )
  
  feature_cols <- obj$features
  
  base_result <- tibble::tibble(
    panel_name = panel_name,
    panel = obj$label,
    condition = condition,
    n_total = nrow(dat),
    n_food_diverse = sum(dat$fcm_profile == DIVERSE_PROFILE),
    n_food_limited = sum(dat$fcm_profile == LIMITED_PROFILE),
    n_features = length(feature_cols),
    permanova_F = NA_real_,
    permanova_R2 = NA_real_,
    permanova_p = NA_real_,
    permdisp_F = NA_real_,
    permdisp_p = NA_real_,
    mean_distance_food_diverse = NA_real_,
    mean_distance_food_limited = NA_real_
  )
  
  if (
    base_result$n_food_diverse < 3L ||
    base_result$n_food_limited < 3L
  ) {
    return(base_result)
  }
  
  X <- finite_numeric_matrix(
    dat,
    feature_cols,
    paste0(panel_name, " diet-profile multivariate")
  )
  
  X_distance <- standardize_for_distance(X)
  distance <- stats::dist(X_distance, method = "euclidean")
  
  meta <- dat %>%
    dplyr::transmute(
      fcm_profile = factor(
        fcm_profile,
        levels = c(LIMITED_PROFILE, DIVERSE_PROFILE)
      )
    )
  
  set.seed(MASTER_SEED + seed_offset)
  permanova <- tryCatch(
    vegan::adonis2(
      distance ~ fcm_profile,
      data = meta,
      permutations = DEEP_DIVE_PERMUTATIONS,
      by = "margin"
    ),
    error = function(e) NULL
  )
  
  if (!is.null(permanova) && "fcm_profile" %in% rownames(permanova)) {
    base_result$permanova_F <- permanova["fcm_profile", "F"]
    base_result$permanova_R2 <- permanova["fcm_profile", "R2"]
    base_result$permanova_p <- permanova["fcm_profile", "Pr(>F)"]
  }
  
  dispersion <- tryCatch(
    vegan::betadisper(
      distance,
      group = meta$fcm_profile,
      type = "median",
      bias.adjust = TRUE
    ),
    error = function(e) NULL
  )
  
  if (!is.null(dispersion)) {
    set.seed(MASTER_SEED + seed_offset + 50000L)
    dispersion_test <- tryCatch(
      vegan::permutest(
        dispersion,
        permutations = DEEP_DIVE_PERMUTATIONS
      ),
      error = function(e) NULL
    )
    
    if (!is.null(dispersion_test)) {
      tab <- dispersion_test$tab
      
      if (nrow(tab) >= 1L) {
        base_result$permdisp_F <- tab[1, "F"]
        base_result$permdisp_p <- tab[1, "Pr(>F)"]
      }
    }
    
    distances <- tibble::tibble(
      profile = as.character(meta$fcm_profile),
      distance_to_centroid = as.numeric(dispersion$distances)
    )
    
    base_result$mean_distance_food_diverse <- mean(
      distances$distance_to_centroid[
        distances$profile == DIVERSE_PROFILE
      ],
      na.rm = TRUE
    )
    
    base_result$mean_distance_food_limited <- mean(
      distances$distance_to_centroid[
        distances$profile == LIMITED_PROFILE
      ],
      na.rm = TRUE
    )
  }
  
  base_result
}

diet_profile_multivariate <- list()
diet_profile_test_index <- 0L

for (condition in c("B_full", "A_excluded")) {
  for (panel_name in names(stage1$panel_results)) {
    diet_profile_test_index <- diet_profile_test_index + 1L
    
    key <- paste(panel_name, condition, sep = "__")
    
    diet_profile_multivariate[[key]] <- run_diet_profile_multivariate(
      panel_name = panel_name,
      condition = condition,
      seed_offset = 200000L + 1000L * diet_profile_test_index
    )
  }
}

diet_profile_multivariate <- dplyr::bind_rows(
  diet_profile_multivariate
) %>%
  dplyr::group_by(condition) %>%
  dplyr::mutate(
    permanova_fdr_across_panels = stats::p.adjust(
      permanova_p,
      method = "BH"
    ),
    permdisp_fdr_across_panels = stats::p.adjust(
      permdisp_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  diet_profile_multivariate,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_diet_profile_permanova_permdisp.csv"
  )
)

# =============================================================================
# 3. SAVE
# =============================================================================

permanova_output <- list(
  multivariate = multivariate_results,
  diet_profile_multivariate = diet_profile_multivariate
)

save_rds_checked(
  permanova_output,
  file.path(DIRS$intermediate, "04b_permanova.rds")
)

cat("\n============================================================\n")
cat("04b PERMANOVA/PERMDISP COMPLETE\n")
cat("============================================================\n")
cat("\nHCC PERMANOVA/PERMDISP summary:\n")
print(
  multivariate_results %>%
    dplyr::select(
      panel, profile, condition, n_total,
      permanova_R2, permanova_p, permanova_fdr,
      permdisp_p, adjusted_permanova_p, adjusted_permanova_fdr
    ),
  n = Inf, width = Inf
)
cat("\nDiet-profile PERMANOVA/PERMDISP summary:\n")
print(
  diet_profile_multivariate %>%
    dplyr::select(
      panel, condition, n_total,
      permanova_R2, permanova_p, permanova_fdr_across_panels,
      permdisp_p
    ),
  n = Inf, width = Inf
)
cat("============================================================\n")