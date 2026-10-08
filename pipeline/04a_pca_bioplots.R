# =============================================================================
# 04a_pca_biplots.R
# PCA biplots for each omics panel (primary Condition A)
#
# PURPOSE
#   Visualise participant-level omics structure in the first two PCs, coloured
#   by disease status and shaped by diet profile. Loading arrows show the top
#   contributing features. Matches Stage 1 settings (mean-centred, not
#   variance-scaled).
#
# INPUTS
#   RDS_PATHS$stage1
#   RDS_PATHS$biomarker
#
# OUTPUTS
#   intermediate/04a_pca_biplots.rds
#   results/04_deep_dive/pca_top_loadings.csv
#   figures/04_deep_dive/pca_biplot_*.png/.pdf
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")
source("pipeline/04_COMMON.R")

set.seed(MASTER_SEED)
ensure_dir(DEEP_DIVE_FIGURE_DIR)

# =============================================================================
# 1. PCA BIPLOTS
# =============================================================================

make_primary_pca_biplot <- function(panel_name) {
  obj <- get_panel_data(
    panel_name,
    condition = "A_excluded"
  )
  
  dat <- obj$data
  feature_cols <- obj$features
  
  if (nrow(dat) < 5L) {
    warning(panel_name, ": too few rows for PCA biplot.", call. = FALSE)
    return(NULL)
  }
  
  X <- finite_numeric_matrix(
    dat,
    feature_cols,
    paste0(panel_name, " primary PCA")
  )
  
  pca <- stats::prcomp(
    X,
    center = TRUE,
    scale. = FALSE
  )
  
  variance <- pca$sdev^2 / sum(pca$sdev^2)
  
  scores <- tibble::tibble(
    study_id = dat$study_id,
    case_status = factor(
      dat$case_status,
      levels = c("Cirrhosis", "HCC")
    ),
    fcm_profile = factor(
      dat$fcm_profile,
      levels = profile_levels
    ),
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2]
  )
  
  loadings <- tibble::tibble(
    feature = rownames(pca$rotation),
    PC1 = pca$rotation[, 1],
    PC2 = pca$rotation[, 2]
  ) %>%
    dplyr::mutate(
      loading_strength = sqrt(PC1^2 + PC2^2)
    ) %>%
    dplyr::slice_max(
      loading_strength,
      n = min(DEEP_DIVE_TOP_LOADINGS, nrow(.)),
      with_ties = FALSE
    )
  
  score_extent <- max(
    abs(c(scores$PC1, scores$PC2)),
    na.rm = TRUE
  )
  loading_extent <- max(
    abs(c(loadings$PC1, loadings$PC2)),
    na.rm = TRUE
  )
  
  arrow_multiplier <- if (
    is.finite(score_extent) &&
    is.finite(loading_extent) &&
    loading_extent > 0
  ) {
    0.75 * score_extent / loading_extent
  } else {
    1
  }
  
  loadings_plot <- loadings %>%
    dplyr::mutate(
      arrow_x = PC1 * arrow_multiplier,
      arrow_y = PC2 * arrow_multiplier
    )
  
  p <- ggplot2::ggplot(
    scores,
    ggplot2::aes(
      x = PC1,
      y = PC2,
      colour = case_status,
      shape = fcm_profile
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey80") +
    ggplot2::geom_point(size = 2.6, alpha = 0.82) +
    ggplot2::geom_segment(
      data = loadings_plot,
      ggplot2::aes(
        x = 0,
        y = 0,
        xend = arrow_x,
        yend = arrow_y
      ),
      inherit.aes = FALSE,
      arrow = grid::arrow(length = grid::unit(0.12, "inches")),
      linewidth = 0.45,
      colour = "grey30"
    ) +
    ggplot2::geom_text(
      data = loadings_plot,
      ggplot2::aes(
        x = arrow_x,
        y = arrow_y,
        label = feature
      ),
      inherit.aes = FALSE,
      size = 3,
      colour = "grey20",
      check_overlap = TRUE,
      vjust = -0.35
    ) +
    ggplot2::labs(
      title = paste0(obj$label, ": primary-condition PCA biplot"),
      subtitle = paste0(
        "Stage 1 T2 flags excluded; mean-centred, not variance-scaled | n = ",
        nrow(dat)
      ),
      x = paste0("PC1 (", round(100 * variance[1], 1), "% variance)"),
      y = paste0("PC2 (", round(100 * variance[2], 1), "% variance)"),
      colour = "Disease status",
      shape = "Diet profile"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      plot.title = ggplot2::element_text(face = "bold")
    )
  
  png_path <- file.path(
    DEEP_DIVE_FIGURE_DIR,
    paste0("pca_biplot_", panel_name, "_A_excluded.png")
  )
  
  pdf_path <- file.path(
    DEEP_DIVE_FIGURE_DIR,
    paste0("pca_biplot_", panel_name, "_A_excluded.pdf")
  )
  
  ggplot2::ggsave(
    png_path,
    p,
    width = 8.5,
    height = 6.5,
    dpi = 300,
    bg = "white"
  )
  
  ggplot2::ggsave(
    pdf_path,
    p,
    width = 8.5,
    height = 6.5,
    bg = "white"
  )
  
  list(
    panel_name = panel_name,
    panel = obj$label,
    n = nrow(dat),
    variance = tibble::tibble(
      component = seq_along(variance),
      variance_fraction = variance,
      cumulative_variance_fraction = cumsum(variance)
    ),
    top_loadings = loadings,
    scores = scores,
    png = png_path,
    pdf = pdf_path
  )
}

pca_outputs <- stats::setNames(
  lapply(
    names(stage1$panel_results),
    make_primary_pca_biplot
  ),
  names(stage1$panel_results)
)

pca_outputs <- pca_outputs[!vapply(pca_outputs, is.null, logical(1))]

pca_top_loadings <- purrr::imap_dfr(
  pca_outputs,
  function(x, panel_name) {
    x$top_loadings %>%
      dplyr::mutate(
        panel_name = panel_name,
        panel = x$panel,
        .before = 1
      )
  }
)

if (nrow(pca_top_loadings) > 0L) {
  readr::write_csv(
    pca_top_loadings,
    file.path(DEEP_DIVE_RESULTS_DIR, "pca_top_loadings.csv")
  )
}

# =============================================================================
# 2. SAVE
# =============================================================================

pca_biplot_output <- list(
  pca = pca_outputs,
  pca_top_loadings = pca_top_loadings
)

save_rds_checked(
  pca_biplot_output,
  file.path(DIRS$intermediate, "04a_pca_biplots.rds")
)

cat("\n============================================================\n")
cat("04a PCA BIPLOTS COMPLETE\n")
cat("Panels:", paste(names(pca_outputs), collapse = ", "), "\n")
cat("============================================================\n")