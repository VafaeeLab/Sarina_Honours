# =============================================================================
# 04c_cytokine_dispersion.R
# Profile-level feature screen, cytokine variability, PCA, and priority
# cytokine construction
#
# PURPOSE
#   1. Aim 2 profile-level feature screen (Food-diverse vs Food-limited) across
#      ALL retained omics features: Wilcoxon + Cliff's delta + BH-FDR.
#   2. Cliff's delta lollipop plots per panel/condition.
#   3. Diet-profile PERMANOVA/PERMDISP across each omics panel (now in 04b).
#   4. Cytokine feature-wise Brown-Forsythe + MAD ratio (primary, outlier-
#      excluded, covariate-residual sensitivity).
#   5. MAD-ratio plots.
#   6. Cytokine PCA on frozen transformed/imputed then z-standardised cytokines.
#   7. Priority-cytokine construction (union of profile-level, dispersion, and
#      PCA loading evidence).
#
# INPUTS
#   RDS_PATHS$stage1
#   RDS_PATHS$biomarker
#
# OUTPUTS
#   intermediate/04c_cytokine_dispersion.rds
#   results/04_deep_dive/aim2_profile_feature_screen_all.csv
#   results/04_deep_dive/aim2_profile_feature_screen_top5_per_panel.csv
#   results/04_deep_dive/aim2_cytokine_brown_forsythe_mad_ratio.csv
#   results/04_deep_dive/aim2_cytokine_pca_scores.csv
#   results/04_deep_dive/aim2_cytokine_pca_loadings.csv
#   results/04_deep_dive/feature_priority_cytokines.csv
#   results/04_deep_dive/feature_priority_cytokine_sources.csv
#   figures/04_deep_dive/aim2_cliffs_delta_*.png/.pdf
#   figures/04_deep_dive/aim2_cytokine_MAD_ratio_*.png/.pdf
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")
source("pipeline/04_COMMON.R")

set.seed(MASTER_SEED)

deep_dive_plot_paths <- list()

run_profile_feature_screen <- function(
    panel_name,
    condition = c("B_full", "A_excluded")
) {
  condition <- match.arg(condition)
  
  obj <- get_panel_data(
    panel_name = panel_name,
    condition = condition
  )
  
  dat <- obj$data
  feature_cols <- obj$features
  
  dat <- dat %>%
    dplyr::filter(
      fcm_profile %in% c(LIMITED_PROFILE, DIVERSE_PROFILE)
    )
  
  n_diverse <- sum(dat$fcm_profile == DIVERSE_PROFILE)
  n_limited <- sum(dat$fcm_profile == LIMITED_PROFILE)
  
  if (n_diverse < 3L || n_limited < 3L) {
    warning(
      panel_name, " / ", condition,
      ": too few participants for profile-level feature screen.",
      call. = FALSE
    )
    return(tibble::tibble())
  }
  
  results <- purrr::map_dfr(
    feature_cols,
    function(feature) {
      x_diverse <- suppressWarnings(
        as.numeric(dat[[feature]][dat$fcm_profile == DIVERSE_PROFILE])
      )
      x_limited <- suppressWarnings(
        as.numeric(dat[[feature]][dat$fcm_profile == LIMITED_PROFILE])
      )
      
      x_diverse <- x_diverse[is.finite(x_diverse)]
      x_limited <- x_limited[is.finite(x_limited)]
      
      if (length(x_diverse) < 3L || length(x_limited) < 3L) {
        return(
          tibble::tibble(
            feature = feature,
            n_food_diverse = length(x_diverse),
            n_food_limited = length(x_limited),
            median_food_diverse = if (length(x_diverse) > 0L) {
              stats::median(x_diverse)
            } else {
              NA_real_
            },
            median_food_limited = if (length(x_limited) > 0L) {
              stats::median(x_limited)
            } else {
              NA_real_
            },
            median_diff_diverse_minus_limited = NA_real_,
            cliffs_delta = NA_real_,
            wilcoxon_p = NA_real_
          )
        )
      }
      
      wt <- tryCatch(
        stats::wilcox.test(
          x_diverse,
          x_limited,
          exact = FALSE,
          correct = TRUE
        ),
        error = function(e) NULL
      )
      
      med_diverse <- stats::median(x_diverse)
      med_limited <- stats::median(x_limited)
      
      tibble::tibble(
        feature = feature,
        n_food_diverse = length(x_diverse),
        n_food_limited = length(x_limited),
        median_food_diverse = med_diverse,
        median_food_limited = med_limited,
        median_diff_diverse_minus_limited = med_diverse - med_limited,
        cliffs_delta = cliffs_delta(x_diverse, x_limited),
        wilcoxon_p = if (!is.null(wt)) wt$p.value else NA_real_
      )
    }
  )
  
  results %>%
    dplyr::mutate(
      panel_name = panel_name,
      panel = obj$label,
      condition = condition,
      fdr_within_panel = stats::p.adjust(wilcoxon_p, method = "BH"),
      direction = dplyr::case_when(
        median_diff_diverse_minus_limited > 0 ~
          "Food-diverse > Food-limited",
        median_diff_diverse_minus_limited < 0 ~
          "Food-diverse < Food-limited",
        is.na(median_diff_diverse_minus_limited) ~ NA_character_,
        TRUE ~ "No median difference"
      ),
      cliffs_magnitude = cliffs_delta_magnitude(cliffs_delta)
    ) %>%
    dplyr::select(
      panel_name,
      panel,
      condition,
      dplyr::everything()
    ) %>%
    dplyr::arrange(wilcoxon_p)
}

profile_feature_results <- tidyr::crossing(
  panel_name = names(stage1$panel_results),
  condition = c("B_full", "A_excluded")
) %>%
  purrr::pmap_dfr(
    function(panel_name, condition) {
      run_profile_feature_screen(panel_name, condition)
    }
  )

if (nrow(profile_feature_results) > 0L) {
  profile_feature_results <- profile_feature_results %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(
      fdr_all_omics = stats::p.adjust(wilcoxon_p, method = "BH")
    ) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(condition, panel_name, wilcoxon_p)
}

DEEP_DIVE_TOP_FEATURES_PER_PANEL <- if (exists(
  "DEEP_DIVE_TOP_FEATURES_PER_PANEL",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_TOP_FEATURES_PER_PANEL)
} else {
  5L
}

profile_feature_top <- profile_feature_results %>%
  dplyr::group_by(panel_name, condition) %>%
  dplyr::slice_min(
    wilcoxon_p,
    n = DEEP_DIVE_TOP_FEATURES_PER_PANEL,
    with_ties = FALSE
  ) %>%
  dplyr::ungroup()

readr::write_csv(
  profile_feature_results,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_profile_feature_screen_all.csv"
  )
)

readr::write_csv(
  profile_feature_top,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_profile_feature_screen_top5_per_panel.csv"
  )
)

# -----------------------------------------------------------------------------
# 5A. CLIFF'S DELTA PLOTS — FOOD-DIVERSE VS FOOD-LIMITED
# -----------------------------------------------------------------------------
# Plot style mirrors the thesis figure convention:
#   - lollipop from zero to Cliff's delta
#   - open circle  = FDR ≥ 0.10
#   - filled circle = FDR < 0.10
#   - positive delta = higher in Food-diverse
# All available features are shown (not only significant features).
# -----------------------------------------------------------------------------

cliffs_plot_paths <- list()

if (nrow(profile_feature_results) > 0L) {
  cliffs_plot_data <- profile_feature_results %>%
    dplyr::filter(is.finite(cliffs_delta)) %>%
    dplyr::mutate(
      fdr_status = dplyr::if_else(
        !is.na(fdr_within_panel) & fdr_within_panel < FDR_THRESHOLD,
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
    )
  
  cliff_groups <- cliffs_plot_data %>%
    dplyr::distinct(panel_name, panel, condition)
  
  for (ii in seq_len(nrow(cliff_groups))) {
    panel_name <- cliff_groups$panel_name[ii]
    panel_label <- cliff_groups$panel[ii]
    condition_name <- cliff_groups$condition[ii]
    
    plot_df <- cliffs_plot_data %>%
      dplyr::filter(
        panel_name == !!panel_name,
        condition == !!condition_name
      ) %>%
      dplyr::arrange(dplyr::desc(cliffs_delta))
    
    if (nrow(plot_df) == 0L) next
    
    # ggplot draws the first factor level at the bottom. Reverse the desired
    # top-to-bottom ordering so the largest positive delta appears at the top.
    desired_order <- plot_df$feature
    plot_df <- plot_df %>%
      dplyr::mutate(
        feature_plot = factor(feature, levels = rev(desired_order))
      )
    
    condition_text <- dplyr::case_when(
      condition_name == "B_full" ~ "Primary: full Cohort A",
      condition_name == "A_excluded" ~ "Sensitivity: Stage 1 T2 flags excluded",
      TRUE ~ condition_name
    )
    
    p_cliffs <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = cliffs_delta, y = feature_plot)
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
          y = feature_plot,
          yend = feature_plot
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
        title = paste0(panel_label, ": Food-diverse vs Food-limited"),
        subtitle = paste0(
          condition_text,
          " | Positive values indicate higher levels in Food-diverse participants"
        ),
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
    
    stub <- paste0(
      "aim2_cliffs_delta_",
      safe_plot_stub(panel_name),
      "_",
      safe_plot_stub(condition_name)
    )
    
    paths <- save_plot_pair(
      p_cliffs,
      stub = stub,
      width = 10.5,
      height = max(6.0, 0.31 * nrow(plot_df) + 2.5)
    )
    
    cliffs_plot_paths[[paste(panel_name, condition_name, sep = "__")]] <- paths
  }
}

# =============================================================================
# 6. DIET-PROFILE PERMANOVA + PERMDISP ACROSS EACH OMICS PANEL
# =============================================================================
# This tests the Food-diverse vs Food-limited multivariate omics structure,
# distinct from the HCC-vs-Cirrhosis PERMANOVA/PERMDISP above.
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
# 7. CYTOKINE FEATURE-WISE DISPERSION: BROWN-FORSYTHE + MAD RATIO
# =============================================================================
# MAD ratio = MAD(Food-diverse) / MAD(Food-limited).
#   > 1 : more variable in Food-diverse
#   < 1 : more variable in Food-limited
#
# Three versions are retained:
#   primary_full       = Cohort A, all participants
#   outlier_excluded   = Stage 1 cytokine T2 flags excluded
#   covariate_residual = residuals after feature ~ profile + age + BMI + disease
# =============================================================================

raw_mad <- function(x) {
  x <- x[is.finite(x)]
  
  if (length(x) == 0L) {
    return(NA_real_)
  }
  
  med <- stats::median(x)
  stats::median(abs(x - med))
}

brown_forsythe_one_feature <- function(
    values,
    groups,
    diverse_label = DIVERSE_PROFILE,
    limited_label = LIMITED_PROFILE
) {
  valid <- is.finite(values) &
    !is.na(groups) &
    groups %in% c(limited_label, diverse_label)
  
  values <- values[valid]
  groups <- factor(
    groups[valid],
    levels = c(limited_label, diverse_label)
  )
  
  n_diverse <- sum(groups == diverse_label)
  n_limited <- sum(groups == limited_label)
  
  if (n_diverse < 3L || n_limited < 3L) {
    return(
      tibble::tibble(
        n_food_diverse = n_diverse,
        n_food_limited = n_limited,
        mad_food_diverse = NA_real_,
        mad_food_limited = NA_real_,
        mad_ratio = NA_real_,
        brown_forsythe_F = NA_real_,
        brown_forsythe_p = NA_real_
      )
    )
  }
  
  group_medians <- tapply(
    values,
    groups,
    stats::median,
    na.rm = TRUE
  )
  
  absolute_deviation <- abs(
    values - group_medians[as.character(groups)]
  )
  
  fit <- tryCatch(
    stats::lm(absolute_deviation ~ groups),
    error = function(e) NULL
  )
  
  anova_table <- if (!is.null(fit)) {
    tryCatch(stats::anova(fit), error = function(e) NULL)
  } else {
    NULL
  }
  
  mad_diverse <- raw_mad(values[groups == diverse_label])
  mad_limited <- raw_mad(values[groups == limited_label])
  
  tibble::tibble(
    n_food_diverse = n_diverse,
    n_food_limited = n_limited,
    mad_food_diverse = mad_diverse,
    mad_food_limited = mad_limited,
    mad_ratio = dplyr::if_else(
      is.finite(mad_limited) && mad_limited > 0,
      mad_diverse / mad_limited,
      NA_real_
    ),
    brown_forsythe_F = if (
      !is.null(anova_table) &&
      nrow(anova_table) >= 1L
    ) {
      anova_table[1, "F value"]
    } else {
      NA_real_
    },
    brown_forsythe_p = if (
      !is.null(anova_table) &&
      nrow(anova_table) >= 1L
    ) {
      anova_table[1, "Pr(>F)"]
    } else {
      NA_real_
    }
  )
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

run_cytokine_dispersion <- function(
    mode = c(
      "primary_full",
      "outlier_excluded",
      "covariate_residual"
    )
) {
  mode <- match.arg(mode)
  
  base_condition <- if (mode == "outlier_excluded") {
    "A_excluded"
  } else {
    "B_full"
  }
  
  obj <- get_cytokine_z_data(base_condition)
  dat <- obj$data
  features <- obj$features
  
  if (mode != "covariate_residual") {
    results <- purrr::map_dfr(
      features,
      function(feature) {
        brown_forsythe_one_feature(
          values = dat[[feature]],
          groups = dat$fcm_profile
        ) %>%
          dplyr::mutate(feature = feature, .before = 1)
      }
    )
  } else {
    residual_data <- dat %>%
      dplyr::left_join(
        clinical %>% dplyr::select(study_id, age, bmi),
        by = "study_id"
      ) %>%
      dplyr::filter(
        is.finite(age),
        is.finite(bmi),
        case_status %in% c("Cirrhosis", "HCC"),
        fcm_profile %in% c(LIMITED_PROFILE, DIVERSE_PROFILE)
      ) %>%
      dplyr::mutate(
        age_z = as.numeric(scale(age)),
        bmi_z = as.numeric(scale(bmi)),
        case_status = factor(
          case_status,
          levels = c("Cirrhosis", "HCC")
        ),
        fcm_profile = factor(
          fcm_profile,
          levels = c(LIMITED_PROFILE, DIVERSE_PROFILE)
        )
      )
    
    results <- purrr::map_dfr(
      features,
      function(feature) {
        model_df <- residual_data %>%
          dplyr::transmute(
            value = .data[[feature]],
            fcm_profile,
            age_z,
            bmi_z,
            case_status
          ) %>%
          dplyr::filter(
            is.finite(value),
            is.finite(age_z),
            is.finite(bmi_z)
          )
        
        if (
          nrow(model_df) < 10L ||
          length(unique(model_df$fcm_profile)) < 2L
        ) {
          return(
            brown_forsythe_one_feature(
              numeric(0),
              character(0)
            ) %>%
              dplyr::mutate(feature = feature, .before = 1)
          )
        }
        
        fit <- tryCatch(
          stats::lm(
            value ~ fcm_profile + age_z + bmi_z + case_status,
            data = model_df
          ),
          error = function(e) NULL
        )
        
        if (is.null(fit)) {
          return(
            brown_forsythe_one_feature(
              numeric(0),
              character(0)
            ) %>%
              dplyr::mutate(feature = feature, .before = 1)
          )
        }
        
        brown_forsythe_one_feature(
          values = stats::residuals(fit),
          groups = model_df$fcm_profile
        ) %>%
          dplyr::mutate(feature = feature, .before = 1)
      }
    )
  }
  
  results %>%
    dplyr::mutate(
      analysis = mode,
      brown_forsythe_fdr = stats::p.adjust(
        brown_forsythe_p,
        method = "BH"
      ),
      variability_direction = dplyr::case_when(
        mad_ratio > 1 ~ "Food-diverse more variable",
        mad_ratio < 1 ~ "Food-limited more variable",
        is.na(mad_ratio) ~ NA_character_,
        TRUE ~ "Equal MAD"
      ),
      .before = 1
    ) %>%
    dplyr::arrange(brown_forsythe_p)
}

cytokine_dispersion <- dplyr::bind_rows(
  run_cytokine_dispersion("primary_full"),
  run_cytokine_dispersion("outlier_excluded"),
  run_cytokine_dispersion("covariate_residual")
)

readr::write_csv(
  cytokine_dispersion,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_cytokine_brown_forsythe_mad_ratio.csv"
  )
)

# -----------------------------------------------------------------------------
# 7A. MAD-RATIO PLOTS — CYTOKINE VARIABILITY
# -----------------------------------------------------------------------------
# Plot style mirrors the prior thesis cytokine figure:
#   - raw MAD ratio on the x-axis
#   - dashed vertical reference line at MAD ratio = 1
#   - each point labelled with Brown-Forsythe raw p and BH-FDR
#   - ratio > 1 = greater variability in Food-diverse
#   - ratio < 1 = greater variability in Food-limited
# -----------------------------------------------------------------------------

mad_plot_paths <- list()

if (nrow(cytokine_dispersion) > 0L) {
  mad_plot_data <- cytokine_dispersion %>%
    dplyr::filter(
      is.finite(mad_ratio),
      mad_ratio > 0
    ) %>%
    dplyr::mutate(
      point_label = paste0(
        "p=", sprintf("%.3f", brown_forsythe_p),
        "\nFDR=", sprintf("%.3f", brown_forsythe_fdr)
      )
    )
  
  for (analysis_name in unique(mad_plot_data$analysis)) {
    plot_df <- mad_plot_data %>%
      dplyr::filter(analysis == !!analysis_name) %>%
      dplyr::arrange(dplyr::desc(mad_ratio))
    
    if (nrow(plot_df) == 0L) next
    
    desired_order <- plot_df$feature
    plot_df <- plot_df %>%
      dplyr::mutate(
        feature_plot = factor(feature, levels = rev(desired_order))
      )
    
    analysis_title <- dplyr::case_when(
      analysis_name == "primary_full" ~ "Before: All cytokines",
      analysis_name == "outlier_excluded" ~ "After: Stage 1 outliers excluded",
      analysis_name == "covariate_residual" ~ "Adjusted: Covariate-residual cytokines",
      TRUE ~ analysis_name
    )
    
    # Expand the right side enough for p/FDR labels while keeping all ratios.
    mad_range <- range(plot_df$mad_ratio, finite = TRUE)
    span <- diff(mad_range)
    if (!is.finite(span) || span <= 0) span <- max(0.5, mad_range[2] * 0.25)
    x_right <- mad_range[2] + 0.55 * span
    x_left <- min(mad_range[1], 1) - 0.08 * span
    
    p_mad <- ggplot2::ggplot(
      plot_df,
      ggplot2::aes(x = mad_ratio, y = feature_plot)
    ) +
      ggplot2::geom_vline(
        xintercept = 1,
        linetype = "dashed",
        linewidth = 0.6,
        colour = "black"
      ) +
      ggplot2::geom_point(
        shape = 21,
        size = 3.3,
        stroke = 0.8,
        fill = "black",
        colour = "black"
      ) +
      ggplot2::geom_text(
        ggplot2::aes(label = point_label),
        hjust = -0.10,
        vjust = 0.5,
        size = 3.6,
        lineheight = 0.9,
        check_overlap = FALSE
      ) +
      ggplot2::coord_cartesian(
        xlim = c(x_left, x_right),
        clip = "off"
      ) +
      ggplot2::labs(
        title = analysis_title,
        subtitle = paste0(
          "MAD ratio > 1 indicates greater variability in ", DIVERSE_PROFILE,
          "; Brown-Forsythe p-values with BH-FDR"
        ),
        x = "MAD ratio: Food-diverse / Food-limited",
        y = NULL
      ) +
      ggplot2::theme_minimal(base_size = 12.5) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_line(colour = "grey90"),
        plot.title = ggplot2::element_text(face = "bold", size = 15),
        plot.subtitle = ggplot2::element_text(size = 10.5),
        axis.title.x = ggplot2::element_text(size = 12),
        axis.text.y = ggplot2::element_text(size = 11),
        plot.margin = ggplot2::margin(10, 55, 10, 10)
      )
    
    stub <- paste0(
      "aim2_cytokine_MAD_ratio_",
      safe_plot_stub(analysis_name)
    )
    
    paths <- save_plot_pair(
      p_mad,
      stub = stub,
      width = 10.5,
      height = max(6.0, 0.34 * nrow(plot_df) + 2.4)
    )
    
    mad_plot_paths[[analysis_name]] <- paths
  }
  
  deep_dive_plot_paths$mad_ratio <- mad_plot_paths
}

# =============================================================================
# 8. CYTOKINE PCA FOR AIM 2 + PRIORITY-CYTOKINE CONSTRUCTION
# =============================================================================
# For the Aim 2 cytokine PCA, use the frozen transformed/imputed cytokines and
# z-standardise each cytokine before PCA. Diet, disease and nutrients do not
# enter the PCA itself.
# =============================================================================

cytokine_z <- get_cytokine_z_data("B_full")
cytokine_z_data <- cytokine_z$data
cytokine_features <- cytokine_z$features

cytokine_X_z <- as.matrix(
  cytokine_z_data[, cytokine_features, drop = FALSE]
)
storage.mode(cytokine_X_z) <- "double"

cytokine_pca <- stats::prcomp(
  cytokine_X_z,
  center = FALSE,
  scale. = FALSE
)

cytokine_pc_variance <- cytokine_pca$sdev^2 /
  sum(cytokine_pca$sdev^2)

cytokine_pc_scores <- tibble::tibble(
  study_id = cytokine_z_data$study_id,
  case_status = cytokine_z_data$case_status,
  fcm_profile = cytokine_z_data$fcm_profile,
  PC1 = cytokine_pca$x[, 1],
  PC2 = cytokine_pca$x[, 2]
)

cytokine_pc_loadings <- tibble::tibble(
  feature = rownames(cytokine_pca$rotation),
  PC1_loading = cytokine_pca$rotation[, 1],
  PC2_loading = cytokine_pca$rotation[, 2]
) %>%
  dplyr::mutate(
    abs_PC1_loading = abs(PC1_loading),
    abs_PC2_loading = abs(PC2_loading)
  )

readr::write_csv(
  cytokine_pc_scores,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_cytokine_pca_scores.csv"
  )
)

readr::write_csv(
  cytokine_pc_loadings,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "aim2_cytokine_pca_loadings.csv"
  )
)

DEEP_DIVE_PRIORITY_TOP_N <- if (exists(
  "DEEP_DIVE_PRIORITY_TOP_N",
  inherits = TRUE
)) {
  as.integer(DEEP_DIVE_PRIORITY_TOP_N)
} else {
  5L
}

cytokine_level_A <- profile_feature_results %>%
  dplyr::filter(
    panel_name == "cytokines",
    condition == "B_full"
  )

cytokine_dispersion_B <- cytokine_dispersion %>%
  dplyr::filter(analysis == "primary_full")

priority_A <- dplyr::bind_rows(
  cytokine_level_A %>%
    dplyr::filter(
      !is.na(fdr_within_panel),
      fdr_within_panel < FDR_THRESHOLD
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "A_level_FDR"
    ),
  cytokine_level_A %>%
    dplyr::slice_min(
      wilcoxon_p,
      n = min(
        DEEP_DIVE_PRIORITY_TOP_N,
        nrow(cytokine_level_A)
      ),
      with_ties = FALSE
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "A_level_top_raw_p"
    )
)

priority_B <- dplyr::bind_rows(
  cytokine_dispersion_B %>%
    dplyr::filter(
      !is.na(brown_forsythe_fdr),
      brown_forsythe_fdr < FDR_THRESHOLD
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "B_dispersion_FDR"
    ),
  cytokine_dispersion_B %>%
    dplyr::slice_min(
      brown_forsythe_p,
      n = min(
        DEEP_DIVE_PRIORITY_TOP_N,
        nrow(cytokine_dispersion_B)
      ),
      with_ties = FALSE
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "B_dispersion_top_raw_p"
    )
)

priority_C <- dplyr::bind_rows(
  cytokine_pc_loadings %>%
    dplyr::slice_max(
      abs_PC1_loading,
      n = min(
        DEEP_DIVE_PRIORITY_TOP_N,
        nrow(cytokine_pc_loadings)
      ),
      with_ties = FALSE
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "C_PC1_loading"
    ),
  cytokine_pc_loadings %>%
    dplyr::slice_max(
      abs_PC2_loading,
      n = min(
        DEEP_DIVE_PRIORITY_TOP_N,
        nrow(cytokine_pc_loadings)
      ),
      with_ties = FALSE
    ) %>%
    dplyr::transmute(
      feature,
      priority_source = "C_PC2_loading"
    )
)

priority_cytokine_sources <- dplyr::bind_rows(
  priority_A,
  priority_B,
  priority_C
) %>%
  dplyr::distinct()

priority_cytokines <- priority_cytokine_sources %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    priority_sources = paste(
      sort(unique(priority_source)),
      collapse = "; "
    ),
    n_priority_sources = dplyr::n_distinct(priority_source),
    .groups = "drop"
  ) %>%
  dplyr::arrange(
    dplyr::desc(n_priority_sources),
    feature
  )

# Add Stage 3 HCC evidence as CONTEXT ONLY. It does not determine D2 inclusion.
hcc_cytokine_context <- biomarker$univariate %>%
  dplyr::filter(
    panel_name == "cytokines",
    condition %in% c("A_excluded", "B_full")
  ) %>%
  dplyr::group_by(feature) %>%
  dplyr::summarise(
    min_hcc_fdr = min(fdr_p, na.rm = TRUE),
    hcc_fdr_any = any(
      !is.na(fdr_p) &
        fdr_p < FDR_THRESHOLD
    ),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    min_hcc_fdr = dplyr::if_else(
      is.infinite(min_hcc_fdr),
      NA_real_,
      min_hcc_fdr
    )
  )

priority_cytokines <- priority_cytokines %>%
  dplyr::left_join(
    hcc_cytokine_context,
    by = "feature"
  )

readr::write_csv(
  priority_cytokines,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "feature_priority_cytokines.csv"
  )
)

readr::write_csv(
  priority_cytokine_sources,
  file.path(
    DEEP_DIVE_RESULTS_DIR,
    "feature_priority_cytokine_sources.csv"
  )
)

# =============================================================================
# SAVE
# =============================================================================

cytokine_dispersion_output <- list(
  profile_feature_screen = profile_feature_results,
  profile_feature_top = profile_feature_top,
  cytokine_dispersion = cytokine_dispersion,
  cytokine_pca = list(
    model = cytokine_pca,
    variance = cytokine_pc_variance,
    scores = cytokine_pc_scores,
    loadings = cytokine_pc_loadings
  ),
  priority_cytokines = priority_cytokines,
  priority_cytokine_sources = priority_cytokine_sources,
  plots = deep_dive_plot_paths
)

save_rds_checked(
  cytokine_dispersion_output,
  file.path(DIRS$intermediate, "04c_cytokine_dispersion.rds")
)

cat("\n============================================================\n")
cat("04c CYTOKINE DISPERSION COMPLETE\n")
cat("============================================================\n")
cat("\nAim 2 profile-level top features:\n")
print(profile_feature_top, n = Inf, width = Inf)

cat("\nCytokine Brown-Forsythe / MAD-ratio top results:\n")
print(
  cytokine_dispersion %>%
    dplyr::group_by(analysis) %>%
    dplyr::slice_min(
      brown_forsythe_p,
      n = 5L,
      with_ties = FALSE
    ) %>%
    dplyr::ungroup(),
  n = Inf, width = Inf
)

cat("\nPriority cytokines for D2:\n")
print(priority_cytokines, n = Inf, width = Inf)
cat("============================================================\n")