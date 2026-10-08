# =============================================================================
# 04d_micronutrient.R
# Micronutrient deep dive: energy adjustment, D0, D1, D2
#
# PURPOSE
#   D0: 18 Willett-energy-adjusted micronutrients by dietary profile.
#   D1: micronutrients vs cytokine PC1/PC2 (Spearman + adjusted regression).
#   D2: micronutrients vs exploratory priority cytokines (Spearman + adjusted).
#
# INPUTS
#   RDS_PATHS$stage1
#   RDS_PATHS$biomarker
#   intermediate/04c_cytokine_dispersion.rds  (priority cytokines, PCA scores)
#   Diet workbook (DIET_SPEC$file)
#
# OUTPUTS
#   intermediate/04d_micronutrient.rds
#   results/04_deep_dive/D0_micronutrient_energy_adjustment_diagnostics.csv
#   results/04_deep_dive/D0_Nutrients_By_Diet_Profile.csv
#   results/04_deep_dive/D1_Nutrients_PC.csv
#   results/04_deep_dive/D1_Adjusted.csv
#   results/04_deep_dive/D2_Nutrients_Cyt.csv
#   results/04_deep_dive/D2_Adjusted.csv
#   figures/04_deep_dive/D0_*.png/.pdf
#   figures/04_deep_dive/D1_*.png/.pdf
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")
source("pipeline/04_COMMON.R")

set.seed(MASTER_SEED)

deep_dive_plot_paths <- list()

# Load cytokine outputs from 04c.
cytokine_04c <- read_rds_checked(
  file.path(DIRS$intermediate, "04c_cytokine_dispersion.rds"),
  "04c cytokine dispersion RDS"
)

cytokine_pca <- cytokine_04c$cytokine_pca$model
cytokine_pc_variance <- cytokine_04c$cytokine_pca$variance
cytokine_pc_scores <- cytokine_04c$cytokine_pca$scores
cytokine_pc_loadings <- cytokine_04c$cytokine_pca$loadings
priority_cytokines <- cytokine_04c$priority_cytokines
profile_feature_results <- cytokine_04c$profile_feature_screen

# Get cytokine z-data for D2.
cytokine_z <- get_cytokine_z_data("B_full")
cytokine_z_data <- cytokine_z$data
cytokine_features <- cytokine_z$features

# =============================================================================
# 9. MICRONUTRIENT PREPARATION: WILLETT ENERGY RESIDUAL METHOD
# =============================================================================
# The 18 micronutrients carried forward from the prior Aim 2 analysis are:
# retinol equivalents, thiamin, riboflavin, niacin equivalents, total folates,
# B6, B12, vitamin C, vitamin E, calcium, iodine, iron, magnesium, phosphorus,
# potassium, selenium, sodium and zinc.
# =============================================================================

if (!requireNamespace("readxl", quietly = TRUE)) {
  stop(
    "Package 'readxl' is required for the micronutrient deep dive. ",
    "Install it with install.packages('readxl').",
    call. = FALSE
  )
}

DEEP_DIVE_MICRONUTRIENTS <- if (exists(
  "DEEP_DIVE_MICRONUTRIENTS",
  inherits = TRUE
)) {
  DEEP_DIVE_MICRONUTRIENTS
} else {
  c(
    "retinol_eq",
    "thiamin",
    "riboflavin",
    "niacin_eq",
    "total_folates",
    "vit_b6",
    "vit_b12",
    "vit_c",
    "vit_e",
    "calcium",
    "iodine",
    "iron",
    "magnesium",
    "phosphorus",
    "potassium",
    "selenium",
    "sodium",
    "zinc"
  )
}

assert_file_exists(
  DIET_SPEC$file,
  "diet workbook for micronutrient analysis"
)

load_diet_nutrients <- function() {
  read_one <- function(sheet_name, case_status) {
    dat <- readxl::read_excel(
      DIET_SPEC$file,
      sheet = sheet_name,
      .name_repair = "minimal"
    )
    
    assert_columns(
      dat,
      c(
        DIET_SPEC$id_col,
        DIET_SPEC$energy_col,
        DEEP_DIVE_MICRONUTRIENTS
      ),
      paste0("Diet sheet ", sheet_name)
    )
    
    out <- tibble::tibble(
      study_id = normalise_id(dat[[DIET_SPEC$id_col]]),
      case_status = case_status,
      total_energy_kj = suppressWarnings(
        as.numeric(dat[[DIET_SPEC$energy_col]])
      )
    )
    
    for (nutrient in DEEP_DIVE_MICRONUTRIENTS) {
      out[[nutrient]] <- suppressWarnings(
        as.numeric(dat[[nutrient]])
      )
    }
    
    out
  }
  
  dplyr::bind_rows(
    read_one(
      DIET_SPEC$sheets[["Cirrhosis"]],
      "Cirrhosis"
    ),
    read_one(
      DIET_SPEC$sheets[["HCC"]],
      "HCC"
    )
  ) %>%
    dplyr::filter(
      study_id %in% assignments$study_id
    ) %>%
    dplyr::mutate(
      total_energy_kcal = total_energy_kj / 4.184
    ) %>%
    dplyr::left_join(
      assignments %>%
        dplyr::select(study_id, fcm_profile),
      by = "study_id"
    )
}

nutrient_raw <- load_diet_nutrients()

if (anyDuplicated(nutrient_raw$study_id)) {
  stop(
    "Micronutrient dataset contains duplicated participant IDs after combining sheets.",
    call. = FALSE
  )
}

willett_adjust_one <- function(
    nutrient,
    energy
) {
  complete <- is.finite(nutrient) &
    is.finite(energy) &
    energy > 0
  
  adjusted <- rep(NA_real_, length(nutrient))
  
  if (
    sum(complete) < 10L ||
    stats::sd(nutrient[complete]) <=
    sqrt(.Machine$double.eps)
  ) {
    return(
      list(
        adjusted = adjusted,
        n_complete = sum(complete),
        r_squared = NA_real_,
        p_energy = NA_real_,
        predicted_at_mean_energy = NA_real_
      )
    )
  }
  
  model_df <- data.frame(
    nutrient = nutrient[complete],
    energy = energy[complete]
  )
  
  fit <- stats::lm(
    nutrient ~ energy,
    data = model_df
  )
  
  mean_energy <- mean(energy[complete])
  
  predicted_at_mean <- as.numeric(
    stats::predict(
      fit,
      newdata = data.frame(
        energy = mean_energy
      )
    )
  )
  
  adjusted[complete] <- stats::residuals(fit) +
    predicted_at_mean
  
  fit_summary <- summary(fit)
  
  list(
    adjusted = adjusted,
    n_complete = sum(complete),
    r_squared = fit_summary$r.squared,
    p_energy = coef(fit_summary)[2, 4],
    predicted_at_mean_energy = predicted_at_mean
  )
}

nutrient_adjusted <- nutrient_raw
nutrient_adjustment_diagnostics <- list()

for (nutrient in DEEP_DIVE_MICRONUTRIENTS) {
  fit <- willett_adjust_one(
    nutrient = nutrient_raw[[nutrient]],
    energy = nutrient_raw$total_energy_kcal
  )
  
  adjusted_name <- paste0(
    nutrient,
    "_energy_adj"
  )
  
  nutrient_adjusted[[adjusted_name]] <- fit$adjusted
  
  nutrient_adjustment_diagnostics[[nutrient]] <- tibble::tibble(
    nutrient = nutrient,
    adjusted_variable = adjusted_name,
    n_complete = fit$n_complete,
    r_squared = fit$r_squared,
    p_energy = fit$p_energy,
    predicted_at_mean_energy = fit$predicted_at_mean_energy
  )
}

nutrient_adjustment_diagnostics <- dplyr::bind_rows(
  nutrient_adjustment_diagnostics
)

adjusted_nutrient_cols <- paste0(
  DEEP_DIVE_MICRONUTRIENTS,
  "_energy_adj"
)

readr::write_csv(
  nutrient_adjustment_diagnostics,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D0_micronutrient_energy_adjustment_diagnostics.csv"
  )
)

# =============================================================================
# 10. D0: ENERGY-ADJUSTED MICRONUTRIENTS BY DIETARY PROFILE
# =============================================================================

D0_micronutrients <- purrr::map_dfr(
  seq_along(DEEP_DIVE_MICRONUTRIENTS),
  function(i) {
    nutrient <- DEEP_DIVE_MICRONUTRIENTS[i]
    adjusted_col <- adjusted_nutrient_cols[i]
    
    x_diverse <- nutrient_adjusted[[adjusted_col]][
      nutrient_adjusted$fcm_profile == DIVERSE_PROFILE
    ]
    
    x_limited <- nutrient_adjusted[[adjusted_col]][
      nutrient_adjusted$fcm_profile == LIMITED_PROFILE
    ]
    
    x_diverse <- x_diverse[is.finite(x_diverse)]
    x_limited <- x_limited[is.finite(x_limited)]
    
    wt <- if (
      length(x_diverse) >= 3L &&
      length(x_limited) >= 3L
    ) {
      tryCatch(
        stats::wilcox.test(
          x_diverse,
          x_limited,
          exact = FALSE,
          correct = TRUE
        ),
        error = function(e) NULL
      )
    } else {
      NULL
    }
    
    med_diverse <- if (length(x_diverse) > 0L) {
      stats::median(x_diverse)
    } else {
      NA_real_
    }
    
    med_limited <- if (length(x_limited) > 0L) {
      stats::median(x_limited)
    } else {
      NA_real_
    }
    
    tibble::tibble(
      nutrient = nutrient,
      adjusted_variable = adjusted_col,
      n_food_diverse = length(x_diverse),
      n_food_limited = length(x_limited),
      median_food_diverse = med_diverse,
      median_food_limited = med_limited,
      median_diff_diverse_minus_limited =
        med_diverse - med_limited,
      cliffs_delta = cliffs_delta(
        x_diverse,
        x_limited
      ),
      wilcoxon_p = if (!is.null(wt)) {
        wt$p.value
      } else {
        NA_real_
      }
    )
  }
) %>%
  dplyr::mutate(
    fdr = stats::p.adjust(
      wilcoxon_p,
      method = "BH"
    ),
    direction = dplyr::case_when(
      median_diff_diverse_minus_limited > 0 ~
        "Food-diverse > Food-limited",
      median_diff_diverse_minus_limited < 0 ~
        "Food-diverse < Food-limited",
      is.na(median_diff_diverse_minus_limited) ~
        NA_character_,
      TRUE ~ "No median difference"
    ),
    cliffs_magnitude = cliffs_delta_magnitude(
      cliffs_delta
    )
  ) %>%
  dplyr::arrange(wilcoxon_p)

readr::write_csv(
  D0_micronutrients,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D0_Nutrients_By_Diet_Profile.csv"
  )
)

# Pretty nutrient labels used by D0 and D1 figures.
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

# -----------------------------------------------------------------------------
# 10A. D0 PLOT — ENERGY-ADJUSTED MICRONUTRIENTS BY DIETARY PROFILE
# -----------------------------------------------------------------------------
# Style mirrors the prior thesis nutrient figure:
#   open circle   = FDR ≥ 0.10
#   filled circle = FDR < 0.10
# Positive Cliff's delta = higher energy-adjusted intake in Food-diverse.
# -----------------------------------------------------------------------------

if (nrow(D0_micronutrients) > 0L) {
  d0_plot_data <- D0_micronutrients %>%
    dplyr::filter(is.finite(cliffs_delta)) %>%
    dplyr::mutate(
      nutrient_label = pretty_nutrient_label(nutrient),
      fdr_status = dplyr::if_else(
        !is.na(fdr) & fdr < FDR_THRESHOLD,
        paste0("FDR < ", sprintf("%.2f", FDR_THRESHOLD)),
        paste0("FDR ≥ ", sprintf("%.2f", FDR_THRESHOLD))
      ),
      fdr_status = factor(
        fdr_status,
        levels = c(
          paste0("FDR ≥ ", sprintf("%.2f", FDR_THRESHOLD)),
          paste0("FDR < ", sprintf("%.2f", FDR_THRESHOLD))
        )
      )
    ) %>%
    dplyr::arrange(dplyr::desc(cliffs_delta))
  
  desired_order <- d0_plot_data$nutrient_label
  d0_plot_data <- d0_plot_data %>%
    dplyr::mutate(
      nutrient_plot = factor(nutrient_label, levels = rev(desired_order))
    )
  
  p_d0_nutrients <- ggplot2::ggplot(
    d0_plot_data,
    ggplot2::aes(x = cliffs_delta, y = nutrient_plot)
  ) +
    ggplot2::geom_vline(
      xintercept = 0,
      linewidth = 0.55,
      colour = "black"
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(
        x = 0,
        xend = cliffs_delta,
        y = nutrient_plot,
        yend = nutrient_plot
      ),
      linewidth = 0.55,
      colour = "grey25"
    ) +
    ggplot2::geom_point(
      ggplot2::aes(fill = fdr_status),
      shape = 21,
      size = 3.2,
      stroke = 0.8,
      colour = "black"
    ) +
    ggplot2::scale_fill_manual(
      values = stats::setNames(
        c("white", "black"),
        c(
          paste0("FDR ≥ ", sprintf("%.2f", FDR_THRESHOLD)),
          paste0("FDR < ", sprintf("%.2f", FDR_THRESHOLD))
        )
      ),
      drop = FALSE
    ) +
    ggplot2::labs(
      title = "Energy-adjusted micronutrients by dietary profile",
      subtitle = "Positive values indicate higher energy-adjusted intake in Food-diverse participants",
      x = "Cliff's delta: Food-diverse vs Food-limited",
      y = NULL,
      fill = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12.5) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
      plot.title = ggplot2::element_text(face = "bold", size = 15),
      plot.subtitle = ggplot2::element_text(size = 11),
      axis.title.x = ggplot2::element_text(size = 12),
      axis.text.y = ggplot2::element_text(size = 11),
      legend.position = "top",
      legend.direction = "horizontal"
    )
  
  d0_nutrient_paths <- save_plot_pair(
    p_d0_nutrients,
    stub = "D0_energy_adjusted_micronutrients_Cliffs_delta",
    width = 10.5,
    height = max(6.5, 0.34 * nrow(d0_plot_data) + 2.2)
  )
  
  deep_dive_plot_paths$D0_micronutrients_cliffs <- d0_nutrient_paths
}

# =============================================================================
# 11. D1: MICRONUTRIENTS VS CYTOKINE PC1/PC2
# =============================================================================
# Primary = Spearman correlation.
# Adjusted follow-up =
#   PC ~ nutrient_z + age_z + bmi_z + case_status + fcm_profile
# =============================================================================

D1_analysis_data <- nutrient_adjusted %>%
  dplyr::select(
    study_id,
    dplyr::all_of(adjusted_nutrient_cols)
  ) %>%
  dplyr::inner_join(
    cytokine_pc_scores %>%
      dplyr::select(
        study_id,
        PC1,
        PC2
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    assignments %>%
      dplyr::select(
        study_id,
        fcm_profile
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    clinical %>%
      dplyr::select(
        study_id,
        age,
        bmi
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    stage1$panel_results$cytokines$processed_cohort_a %>%
      dplyr::select(
        study_id,
        case_status
      ),
    by = "study_id"
  )

run_spearman <- function(x, y) {
  complete <- is.finite(x) &
    is.finite(y)
  
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

D1_Nutrient_PC <- purrr::map_dfr(
  seq_along(DEEP_DIVE_MICRONUTRIENTS),
  function(i) {
    nutrient <- DEEP_DIVE_MICRONUTRIENTS[i]
    adjusted_col <- adjusted_nutrient_cols[i]
    
    purrr::map_dfr(
      c("PC1", "PC2"),
      function(pc) {
        test <- run_spearman(
          D1_analysis_data[[adjusted_col]],
          D1_analysis_data[[pc]]
        )
        
        tibble::tibble(
          nutrient = nutrient,
          adjusted_variable = adjusted_col,
          PC = pc,
          n_complete = test$n,
          spearman_rho = test$rho,
          spearman_p = test$p
        )
      }
    )
  }
) %>%
  dplyr::group_by(PC) %>%
  dplyr::mutate(
    FDR_within_PC = stats::p.adjust(
      spearman_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_all_PC_tests = stats::p.adjust(
      spearman_p,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(spearman_p)

fit_D1_adjusted <- function(
    adjusted_col,
    nutrient,
    pc
) {
  dat <- D1_analysis_data %>%
    dplyr::transmute(
      outcome = .data[[pc]],
      nutrient_value = .data[[adjusted_col]],
      age,
      bmi,
      case_status,
      fcm_profile
    ) %>%
    dplyr::filter(
      is.finite(outcome),
      is.finite(nutrient_value),
      is.finite(age),
      is.finite(bmi),
      case_status %in% c("Cirrhosis", "HCC"),
      fcm_profile %in% c(LIMITED_PROFILE, DIVERSE_PROFILE)
    )
  
  empty <- function() {
    tibble::tibble(
      nutrient = nutrient,
      adjusted_variable = adjusted_col,
      PC = pc,
      n_complete = nrow(dat),
      beta_nutrient = NA_real_,
      beta_ci_low = NA_real_,
      beta_ci_high = NA_real_,
      nutrient_p = NA_real_,
      model_r_squared = NA_real_
    )
  }
  
  if (
    nrow(dat) < 15L ||
    stats::sd(dat$nutrient_value) <=
    sqrt(.Machine$double.eps)
  ) {
    return(empty())
  }
  
  dat <- dat %>%
    dplyr::mutate(
      nutrient_z = as.numeric(
        scale(nutrient_value)
      ),
      age_z = as.numeric(scale(age)),
      bmi_z = as.numeric(scale(bmi)),
      case_status = factor(
        case_status,
        levels = c("Cirrhosis", "HCC")
      ),
      fcm_profile = factor(
        fcm_profile,
        levels = c(
          LIMITED_PROFILE,
          DIVERSE_PROFILE
        )
      )
    )
  
  fit <- tryCatch(
    stats::lm(
      outcome ~ nutrient_z + age_z + bmi_z +
        case_status + fcm_profile,
      data = dat
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(empty())
  }
  
  coefs <- summary(fit)$coefficients
  
  if (!"nutrient_z" %in% rownames(coefs)) {
    return(empty())
  }
  
  estimate <- coefs["nutrient_z", "Estimate"]
  se <- coefs["nutrient_z", "Std. Error"]
  
  tibble::tibble(
    nutrient = nutrient,
    adjusted_variable = adjusted_col,
    PC = pc,
    n_complete = nrow(dat),
    beta_nutrient = estimate,
    beta_ci_low = estimate - 1.96 * se,
    beta_ci_high = estimate + 1.96 * se,
    nutrient_p = coefs[
      "nutrient_z",
      "Pr(>|t|)"
    ],
    model_r_squared = summary(fit)$r.squared
  )
}

D1_Adjusted <- purrr::map_dfr(
  seq_along(DEEP_DIVE_MICRONUTRIENTS),
  function(i) {
    nutrient <- DEEP_DIVE_MICRONUTRIENTS[i]
    adjusted_col <- adjusted_nutrient_cols[i]
    
    purrr::map_dfr(
      c("PC1", "PC2"),
      function(pc) {
        fit_D1_adjusted(
          adjusted_col,
          nutrient,
          pc
        )
      }
    )
  }
) %>%
  dplyr::group_by(PC) %>%
  dplyr::mutate(
    FDR_within_PC = stats::p.adjust(
      nutrient_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_all_PC_tests = stats::p.adjust(
      nutrient_p,
      method = "BH"
    )
  ) %>%
  dplyr::arrange(nutrient_p)

readr::write_csv(
  D1_Nutrient_PC,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D1_Nutrients_PC.csv"
  )
)

readr::write_csv(
  D1_Adjusted,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D1_Adjusted.csv"
  )
)

# -----------------------------------------------------------------------------
# 11A. D1 PLOTS — MICRONUTRIENTS VS CYTOKINE PC1/PC2
# -----------------------------------------------------------------------------
# Two complementary heatmaps are saved:
#   1. Primary Spearman rho (the prespecified D1 primary analysis).
#   2. Adjusted nutrient regression coefficient from PC ~ nutrient + age + BMI + disease + profile.
# Stars are based on BH-FDR across the complete 18 x 2 test family:
#   *  FDR < 0.10
#   ** FDR < 0.05
# -----------------------------------------------------------------------------

format_fdr_star <- function(fdr) {
  dplyr::case_when(
    !is.na(fdr) & fdr < 0.05 ~ "**",
    !is.na(fdr) & fdr < FDR_THRESHOLD ~ "*",
    TRUE ~ ""
  )
}

# ---- D1 primary Spearman heatmap --------------------------------------------
if (nrow(D1_Nutrient_PC) > 0L) {
  d1_spearman_plot_data <- D1_Nutrient_PC %>%
    dplyr::mutate(
      nutrient_label = pretty_nutrient_label(nutrient),
      PC = factor(PC, levels = c("PC1", "PC2")),
      significance = format_fdr_star(FDR_all_PC_tests),
      cell_label = dplyr::if_else(
        is.finite(spearman_rho),
        paste0(sprintf("%+.2f", spearman_rho), significance),
        "NA"
      )
    )
  
  nutrient_order <- d1_spearman_plot_data %>%
    dplyr::group_by(nutrient_label) %>%
    dplyr::summarise(
      max_abs_rho = max(abs(spearman_rho), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(max_abs_rho) %>%
    dplyr::pull(nutrient_label)
  
  d1_spearman_plot_data <- d1_spearman_plot_data %>%
    dplyr::mutate(
      nutrient_label = factor(
        nutrient_label,
        levels = nutrient_order
      )
    )
  
  rho_limit <- max(abs(d1_spearman_plot_data$spearman_rho), na.rm = TRUE)
  if (!is.finite(rho_limit) || rho_limit <= 0) rho_limit <- 1
  rho_limit <- min(1, max(0.25, rho_limit))
  
  p_d1_spearman <- ggplot2::ggplot(
    d1_spearman_plot_data,
    ggplot2::aes(x = PC, y = nutrient_label, fill = spearman_rho)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = cell_label),
      size = 3.4,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-rho_limit, rho_limit),
      oob = scales::squish,
      name = "Spearman rho"
    ) +
    ggplot2::labs(
      title = "D1: Energy-adjusted micronutrients vs cytokine PCA axes",
      subtitle = paste0(
        "Primary Spearman analysis | * BH-FDR < ", FDR_THRESHOLD,
        ", ** BH-FDR < 0.05 across all nutrient × PC tests"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  d1_spearman_paths <- save_plot_pair(
    p_d1_spearman,
    stub = "D1_micronutrients_vs_cytokine_PCs_Spearman_heatmap",
    width = 8.5,
    height = 8.5
  )
  
  deep_dive_plot_paths$D1_spearman <- d1_spearman_paths
}

# ---- D1 adjusted regression-coefficient heatmap -------------------------------
if (nrow(D1_Adjusted) > 0L) {
  d1_adjusted_plot_data <- D1_Adjusted %>%
    dplyr::mutate(
      nutrient_label = pretty_nutrient_label(nutrient),
      PC = factor(PC, levels = c("PC1", "PC2")),
      significance = format_fdr_star(FDR_all_PC_tests),
      cell_label = dplyr::if_else(
        is.finite(beta_nutrient),
        paste0(sprintf("%+.2f", beta_nutrient), significance),
        "NA"
      )
    )
  
  nutrient_order_beta <- d1_adjusted_plot_data %>%
    dplyr::group_by(nutrient_label) %>%
    dplyr::summarise(
      max_abs_beta = max(abs(beta_nutrient), na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::arrange(max_abs_beta) %>%
    dplyr::pull(nutrient_label)
  
  d1_adjusted_plot_data <- d1_adjusted_plot_data %>%
    dplyr::mutate(
      nutrient_label = factor(
        nutrient_label,
        levels = nutrient_order_beta
      )
    )
  
  beta_limit <- max(abs(d1_adjusted_plot_data$beta_nutrient), na.rm = TRUE)
  if (!is.finite(beta_limit) || beta_limit <= 0) beta_limit <- 1
  
  p_d1_adjusted <- ggplot2::ggplot(
    d1_adjusted_plot_data,
    ggplot2::aes(x = PC, y = nutrient_label, fill = beta_nutrient)
  ) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.8) +
    ggplot2::geom_text(
      ggplot2::aes(label = cell_label),
      size = 3.4,
      fontface = "bold"
    ) +
    ggplot2::scale_fill_gradient2(
      low = "#4575B4",
      mid = "white",
      high = "#D73027",
      midpoint = 0,
      limits = c(-beta_limit, beta_limit),
      oob = scales::squish,
      name = "Adjusted coefficient (beta)"
    ) +
    ggplot2::labs(
      title = "D1: Covariate-adjusted micronutrient associations with cytokine PCA axes",
      subtitle = paste0(
        "Coefficient = change in PC score per 1-SD higher nutrient, adjusted for age, BMI, disease and diet profile | * BH-FDR < ",
        FDR_THRESHOLD,
        ", ** BH-FDR < 0.05"
      ),
      x = NULL,
      y = NULL
    ) +
    ggplot2::theme_minimal(base_size = 11.5) +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold"),
      axis.text.x = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    )
  
  d1_adjusted_paths <- save_plot_pair(
    p_d1_adjusted,
    stub = "D1_micronutrients_vs_cytokine_PCs_adjusted_coefficient_heatmap",
    width = 8.5,
    height = 8.5
  )
  
  deep_dive_plot_paths$D1_adjusted <- d1_adjusted_paths
}

# =============================================================================
# 12. D2: MICRONUTRIENTS VS EXPLORATORY PRIORITY CYTOKINES
# =============================================================================
# Primary = Spearman correlation across all nutrient x priority-cytokine pairs.
# Adjusted follow-up =
#   cytokine_z ~ nutrient_z + age_z + bmi_z + case_status + fcm_profile
# BH-FDR is reported across the full tested pair family.
# =============================================================================

priority_feature_names <- intersect(
  priority_cytokines$feature,
  cytokine_features
)

D2_analysis_data <- nutrient_adjusted %>%
  dplyr::select(
    study_id,
    dplyr::all_of(adjusted_nutrient_cols)
  ) %>%
  dplyr::inner_join(
    cytokine_z_data %>%
      dplyr::select(
        study_id,
        case_status,
        dplyr::all_of(priority_feature_names)
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    assignments %>%
      dplyr::select(
        study_id,
        fcm_profile
      ),
    by = "study_id"
  ) %>%
  dplyr::left_join(
    clinical %>%
      dplyr::select(
        study_id,
        age,
        bmi
      ),
    by = "study_id"
  )

D2_Nutrients_Cyt <- purrr::map_dfr(
  seq_along(DEEP_DIVE_MICRONUTRIENTS),
  function(i) {
    nutrient <- DEEP_DIVE_MICRONUTRIENTS[i]
    adjusted_col <- adjusted_nutrient_cols[i]
    
    purrr::map_dfr(
      priority_feature_names,
      function(feature) {
        test <- run_spearman(
          D2_analysis_data[[adjusted_col]],
          D2_analysis_data[[feature]]
        )
        
        tibble::tibble(
          nutrient = nutrient,
          adjusted_variable = adjusted_col,
          cytokine = feature,
          n_complete = test$n,
          spearman_rho = test$rho,
          spearman_p = test$p
        )
      }
    )
  }
) %>%
  dplyr::group_by(cytokine) %>%
  dplyr::mutate(
    FDR_within_cytokine = stats::p.adjust(
      spearman_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_all_nutrient_cytokine_tests =
      stats::p.adjust(
        spearman_p,
        method = "BH"
      )
  ) %>%
  dplyr::arrange(spearman_p)

fit_D2_adjusted <- function(
    adjusted_col,
    nutrient,
    feature
) {
  dat <- D2_analysis_data %>%
    dplyr::transmute(
      cytokine_value = .data[[feature]],
      nutrient_value = .data[[adjusted_col]],
      age,
      bmi,
      case_status,
      fcm_profile
    ) %>%
    dplyr::filter(
      is.finite(cytokine_value),
      is.finite(nutrient_value),
      is.finite(age),
      is.finite(bmi),
      case_status %in% c("Cirrhosis", "HCC"),
      fcm_profile %in% c(LIMITED_PROFILE, DIVERSE_PROFILE)
    )
  
  empty <- function() {
    tibble::tibble(
      nutrient = nutrient,
      adjusted_variable = adjusted_col,
      cytokine = feature,
      n_complete = nrow(dat),
      beta_nutrient = NA_real_,
      beta_ci_low = NA_real_,
      beta_ci_high = NA_real_,
      nutrient_p = NA_real_,
      model_r_squared = NA_real_
    )
  }
  
  if (
    nrow(dat) < 15L ||
    stats::sd(dat$nutrient_value) <=
    sqrt(.Machine$double.eps) ||
    stats::sd(dat$cytokine_value) <=
    sqrt(.Machine$double.eps)
  ) {
    return(empty())
  }
  
  dat <- dat %>%
    dplyr::mutate(
      nutrient_z = as.numeric(
        scale(nutrient_value)
      ),
      age_z = as.numeric(scale(age)),
      bmi_z = as.numeric(scale(bmi)),
      case_status = factor(
        case_status,
        levels = c("Cirrhosis", "HCC")
      ),
      fcm_profile = factor(
        fcm_profile,
        levels = c(
          LIMITED_PROFILE,
          DIVERSE_PROFILE
        )
      )
    )
  
  fit <- tryCatch(
    stats::lm(
      cytokine_value ~ nutrient_z +
        age_z + bmi_z +
        case_status + fcm_profile,
      data = dat
    ),
    error = function(e) NULL
  )
  
  if (is.null(fit)) {
    return(empty())
  }
  
  coefs <- summary(fit)$coefficients
  
  if (!"nutrient_z" %in% rownames(coefs)) {
    return(empty())
  }
  
  estimate <- coefs["nutrient_z", "Estimate"]
  se <- coefs["nutrient_z", "Std. Error"]
  
  tibble::tibble(
    nutrient = nutrient,
    adjusted_variable = adjusted_col,
    cytokine = feature,
    n_complete = nrow(dat),
    beta_nutrient = estimate,
    beta_ci_low = estimate - 1.96 * se,
    beta_ci_high = estimate + 1.96 * se,
    nutrient_p = coefs[
      "nutrient_z",
      "Pr(>|t|)"
    ],
    model_r_squared = summary(fit)$r.squared
  )
}

D2_Adjusted <- purrr::map_dfr(
  seq_along(DEEP_DIVE_MICRONUTRIENTS),
  function(i) {
    nutrient <- DEEP_DIVE_MICRONUTRIENTS[i]
    adjusted_col <- adjusted_nutrient_cols[i]
    
    purrr::map_dfr(
      priority_feature_names,
      function(feature) {
        fit_D2_adjusted(
          adjusted_col,
          nutrient,
          feature
        )
      }
    )
  }
) %>%
  dplyr::group_by(cytokine) %>%
  dplyr::mutate(
    FDR_within_cytokine = stats::p.adjust(
      nutrient_p,
      method = "BH"
    )
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    FDR_all_nutrient_cytokine_tests =
      stats::p.adjust(
        nutrient_p,
        method = "BH"
      )
  ) %>%
  dplyr::arrange(nutrient_p)

readr::write_csv(
  D2_Nutrients_Cyt,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D2_Nutrients_Cyt.csv"
  )
)

readr::write_csv(
  D2_Adjusted,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "D2_Adjusted.csv"
  )
)

# =============================================================================
# SAVE
# =============================================================================

micronutrient_output <- list(
  nutrient_adjustment_diagnostics = nutrient_adjustment_diagnostics,
  D0_micronutrients_by_profile = D0_micronutrients,
  D1_nutrients_PC = D1_Nutrient_PC,
  D1_adjusted = D1_Adjusted,
  D2_nutrients_cytokines = D2_Nutrients_Cyt,
  D2_adjusted = D2_Adjusted,
  plots = deep_dive_plot_paths
)

save_rds_checked(
  micronutrient_output,
  file.path(DIRS$intermediate, "04d_micronutrient.rds")
)

cat("\n============================================================\n")
cat("04d MICRONUTRIENT DEEP DIVE COMPLETE\n")
cat("============================================================\n")
cat("\nD1 nutrient-PC top associations:\n")
print(
  D1_Nutrient_PC %>%
    dplyr::slice_min(spearman_p, n = 10L, with_ties = FALSE),
  n = Inf, width = Inf
)
cat("\nD2 nutrient-cytokine top associations:\n")
print(
  D2_Nutrients_Cyt %>%
    dplyr::slice_min(spearman_p, n = 10L, with_ties = FALSE),
  n = Inf, width = Inf
)
cat("============================================================\n")