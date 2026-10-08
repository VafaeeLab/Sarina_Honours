# =============================================================================
# utils.R
# Shared helper functions for the thesis analysis pipeline.
#
# Source config.R before this file:
#   source("pipeline/config.R")
#   source("pipeline/utils.R")
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Small general helpers
# -----------------------------------------------------------------------------

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) y else x
}

assert_file_exists <- function(path, label = NULL) {
  if (!file.exists(path)) {
    if (is.null(label)) label <- path
    stop("Required file not found: ", label, " [", path, "]", call. = FALSE)
  }
  invisible(path)
}

assert_columns <- function(data, required, label = deparse(substitute(data))) {
  missing_cols <- setdiff(required, names(data))
  if (length(missing_cols) > 0L) {
    stop(
      label, " is missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

ensure_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  invisible(path)
}

ensure_parent_dir <- function(path) {
  ensure_dir(dirname(path))
  invisible(path)
}

read_rds_checked <- function(path, label = NULL) {
  assert_file_exists(path, label)
  readRDS(path)
}

save_rds_checked <- function(object, path) {
  ensure_parent_dir(path)
  saveRDS(object, path)
  invisible(path)
}

# -----------------------------------------------------------------------------
# 2. Participant-ID handling
# -----------------------------------------------------------------------------

# Matches the current ID normalisation: convert to character, replace NBSP,
# collapse whitespace, uppercase, and convert empty strings to NA.
normalise_id <- function(x) {
  out <- as.character(x)
  out <- stringr::str_replace_all(out, "\u00A0", " ")
  out <- stringr::str_squish(out)
  out <- stringr::str_to_upper(out)
  out[out == ""] <- NA_character_
  out
}

assert_unique_ids <- function(data, id_col = "study_id", label = "data") {
  assert_columns(data, id_col, label)
  
  ids <- data[[id_col]]
  duplicated_ids <- unique(ids[!is.na(ids) & duplicated(ids)])
  
  if (length(duplicated_ids) > 0L) {
    stop(
      label, " contains duplicated ", id_col, " values: ",
      paste(utils::head(duplicated_ids, 10L), collapse = ", "),
      if (length(duplicated_ids) > 10L) " ..." else "",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

# -----------------------------------------------------------------------------
# 3. FFQ helpers
# -----------------------------------------------------------------------------

normalise_response <- function(x) {
  out <- as.character(x)
  out <- stringr::str_squish(out)
  out <- stringr::str_to_lower(out)
  out[out %in% c("", "na", "n/a", "missing", "not answered")] <- NA_character_
  out
}

# User-selected food-group definition: mean daily frequency across the FFQ
# items belonging to the group. Only columns actually present are used.
safe_rowmean <- function(data, cols) {
  present <- intersect(cols, names(data))
  
  if (length(present) == 0L) {
    return(rep(NA_real_, nrow(data)))
  }
  
  rowMeans(
    as.data.frame(data[, present, drop = FALSE]),
    na.rm = TRUE
  )
}

# -----------------------------------------------------------------------------
# 4. Numeric-feature helpers
# -----------------------------------------------------------------------------

median_impute_numeric <- function(data, cols = names(data), label = "data") {
  assert_columns(data, cols, label)
  
  out <- data
  
  for (col in cols) {
    x <- suppressWarnings(as.numeric(out[[col]]))
    
    if (all(is.na(x))) {
      stop(
        label, ": cannot median-impute all-NA column '", col, "'.",
        call. = FALSE
      )
    }
    
    med <- stats::median(x, na.rm = TRUE)
    x[is.na(x)] <- med
    out[[col]] <- x
  }
  
  out
}

coerce_numeric_features <- function(data, cols, label = "data") {
  assert_columns(data, cols, label)
  
  out <- data
  introduced_na <- integer(length(cols))
  names(introduced_na) <- cols
  
  for (col in cols) {
    before_na <- is.na(out[[col]])
    converted <- suppressWarnings(as.numeric(out[[col]]))
    introduced_na[col] <- sum(is.na(converted) & !before_na)
    out[[col]] <- converted
  }
  
  bad <- introduced_na[introduced_na > 0L]
  if (length(bad) > 0L) {
    warning(
      label, ": numeric coercion introduced NA values in: ",
      paste(names(bad), bad, sep = "=", collapse = ", "),
      call. = FALSE
    )
  }
  
  out
}

# -----------------------------------------------------------------------------
# 5. Omics-panel helpers
# -----------------------------------------------------------------------------

get_panel_spec <- function(panel_name, panel_specs = OMICS_PANELS) {
  if (!panel_name %in% names(panel_specs)) {
    stop(
      "Unknown panel '", panel_name, "'. Valid panels: ",
      paste(names(panel_specs), collapse = ", "),
      call. = FALSE
    )
  }
  
  panel_specs[[panel_name]]
}

panel_uses_forced_log <- function(panel_name,
                                  force_log_panels = FORCE_LOG_PANELS) {
  panel_name %in% force_log_panels
}

# Load one omics panel in a common participant-level format.
# This function performs structural cleaning only:
#   * reads Cirrhosis + HCC files
#   * harmonises the ID column
#   * normalises IDs
#   * filters to the requested visit
#   * adds case_status
#   * converts feature columns to numeric

load_omics_panel <- function(panel_name,
                             panel_specs = OMICS_PANELS,
                             visit = OMICS_VISIT,
                             restrict_ids = NULL) {
  
  spec <- get_panel_spec(panel_name, panel_specs)
  
  assert_file_exists(spec$cirrhosis_file, paste0(spec$label, " Cirrhosis"))
  assert_file_exists(spec$hcc_file, paste0(spec$label, " HCC"))
  
  cirr <- readr::read_csv(
    spec$cirrhosis_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  hcc <- readr::read_csv(
    spec$hcc_file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  assert_columns(
    cirr,
    c(spec$cirrhosis_id_col, spec$visit_col),
    paste0(spec$label, " Cirrhosis")
  )
  assert_columns(
    hcc,
    c(spec$hcc_id_col, spec$visit_col),
    paste0(spec$label, " HCC")
  )
  
  names(cirr)[names(cirr) == spec$cirrhosis_id_col] <- "study_id"
  names(hcc)[names(hcc) == spec$hcc_id_col] <- "study_id"
  
  cirr <- cirr |>
    dplyr::mutate(
      study_id = normalise_id(study_id),
      case_status = "Cirrhosis"
    ) |>
    dplyr::filter(.data[[spec$visit_col]] == visit)
  
  hcc <- hcc |>
    dplyr::mutate(
      study_id = normalise_id(study_id),
      case_status = "HCC"
    ) |>
    dplyr::filter(.data[[spec$visit_col]] == visit)
  
  combined <- dplyr::bind_rows(cirr, hcc)
  
  if (!is.null(restrict_ids)) {
    restrict_ids <- normalise_id(restrict_ids)
    combined <- combined |>
      dplyr::filter(study_id %in% restrict_ids)
  }
  
  meta_cols <- c("study_id", "case_status", spec$visit_col)
  feature_cols <- setdiff(names(combined), meta_cols)
  
  if (length(feature_cols) == 0L) {
    stop(spec$label, " has no feature columns after metadata removal.", call. = FALSE)
  }
  
  combined <- coerce_numeric_features(
    combined,
    cols = feature_cols,
    label = spec$label
  )
  
  # Remove participants with no observed values across the entire panel.
  keep <- rowSums(!is.na(combined[, feature_cols, drop = FALSE])) > 0L
  combined <- combined[keep, , drop = FALSE]
  
  list(
    panel_name   = panel_name,
    label        = spec$label,
    data         = combined,
    meta_cols    = meta_cols,
    feature_cols = feature_cols,
    visit        = visit
  )
}

# Apply log1p only when the panel is explicitly listed in FORCE_LOG_PANELS.
# Negative values cause an error
apply_configured_log1p <- function(data,
                                   feature_cols,
                                   panel_name,
                                   force_log_panels = FORCE_LOG_PANELS) {
  
  assert_columns(data, feature_cols, panel_name)
  apply_log <- panel_uses_forced_log(panel_name, force_log_panels)
  
  if (!apply_log) {
    return(list(data = data, log_applied = FALSE))
  }
  
  x <- as.matrix(data[, feature_cols, drop = FALSE])
  storage.mode(x) <- "double"
  
  if (any(x < 0, na.rm = TRUE)) {
    negative_cols <- feature_cols[
      vapply(data[, feature_cols, drop = FALSE], function(z) {
        any(suppressWarnings(as.numeric(z)) < 0, na.rm = TRUE)
      }, logical(1))
    ]
    
    stop(
      "Cannot apply log1p to panel '", panel_name,
      "' because negative values were found in: ",
      paste(negative_cols, collapse = ", "),
      call. = FALSE
    )
  }
  
  out <- data
  out[, feature_cols] <- lapply(
    out[, feature_cols, drop = FALSE],
    function(z) log1p(as.numeric(z))
  )
  
  list(data = out, log_applied = TRUE)
}

# =============================================================================
# End utils.R
# =============================================================================
