# =============================================================================
# 0A_Goldberg.R
# Low-energy reporter (LER) classification using the Goldberg method.
#
# PURPOSE
#   Classify participants as low-energy reporters based on EI:BMR ratio using
#   the Mifflin-St Jeor equation for BMR estimation. This runs BEFORE 0_preprocess.R
#   and produces the Goldberg classification file consumed by the preprocessing step.
#
# INPUTS
#   DATA_DIR / Diet_HCC.csv
#   DATA_DIR / Diet_Cirr.csv
#   DATA_DIR / Clin_Path_MAFLD_OnlyHCC.csv          (via CLINICAL_SPEC)
#   DATA_DIR / Clin_Path_MAFLD_OnlyCirrhosis.csv     (via CLINICAL_SPEC)
#
# OUTPUTS
#   DATA_DIR / goldberg_reporter_classification.csv   (consumed by 0_preprocess.R)
#   DIRS$results / 0A_goldberg / goldberg_summary.csv
#
# This script is standalone and should be run before 0_preprocess.R.
# =============================================================================

source("pipeline/config.R")
source("pipeline/utils.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

set.seed(MASTER_SEED)

# =============================================================================
# 1. OUTPUT DIRECTORY
# =============================================================================

GOLDBERG_RESULTS_DIR <- file.path(DIRS$results, "0A_goldberg")
ensure_dir(GOLDBERG_RESULTS_DIR)

# =============================================================================
# 2. GOLDBERG-SPECIFIC HELPERS
# =============================================================================

#' Build a string listing data issues for a participant
make_issue_string <- function(has_energy, has_age, has_bmi, has_sex) {
  issues <- character(0)
  if (!has_energy) issues <- c(issues, "missing energy intake")
  if (!has_age)    issues <- c(issues, "missing age")
  if (!has_bmi)    issues <- c(issues, "missing BMI")
  if (!has_sex)    issues <- c(issues, "missing sex")
  if (length(issues) == 0L) return(NA_character_)
  paste(issues, collapse = "; ")
}

#' Check that the expected number of rows matches
check_expected_n <- function(data, expected_n, label) {
  actual_n <- nrow(data)
  if (!is.null(expected_n) && actual_n != expected_n) {
    warning(
      label, ": expected ", expected_n, " rows but got ", actual_n,
      call. = FALSE
    )
  }
  invisible(data)
}

# =============================================================================
# 3. MIFFLIN-ST JEOR BMR ESTIMATION
# =============================================================================
# BMR (kcal/day):
#   Male:   10 * weight_kg + 6.25 * height_cm - 5 * age - 161 + 166
#           => 10 * weight_kg + 6.25 * height_cm - 5 * age + 5
#   Female: 10 * weight_kg + 6.25 * height_cm - 5 * age - 161
#
# Since we have BMI but not raw weight and height, we use:
#   weight_kg ≈ BMI * height_m^2
# =============================================================================

DEFAULT_HEIGHT_M <- list(
  male   = 1.75,
  female = 1.62
)

estimate_bmr_mifflin <- function(age, bmi, sex, height_m = NULL) {
  sex_lower <- tolower(trimws(sex))

  if (is.null(height_m) || is.na(height_m)) {
    height_m <- if (sex_lower == "male") {
      DEFAULT_HEIGHT_M$male
    } else {
      DEFAULT_HEIGHT_M$female
    }
  }

  weight_kg <- bmi * height_m^2
  height_cm <- height_m * 100

  bmr <- 10 * weight_kg + 6.25 * height_cm - 5 * age

  if (sex_lower == "male") {
    bmr <- bmr + 5
  } else {
    bmr <- bmr - 161
  }

  bmr
}

# =============================================================================
# 4. GOLDBERG EI:BMR CLASSIFICATION
# =============================================================================
# Goldberg cutoff: EI:BMR < 1.1 => low-energy reporter (LER)
# =============================================================================

GOLDBERG_CUTOFF <- 1.1

classify_goldberg <- function(ei_kcal, bmr_kcal) {
  if (is.na(ei_kcal) || is.na(bmr_kcal) || bmr_kcal <= 0) {
    return(NA_character_)
  }
  ratio <- ei_kcal / bmr_kcal
  if (ratio < GOLDBERG_CUTOFF) {
    return("low energy reporter (LER)")
  } else {
    return("plausible reporter")
  }
}

# =============================================================================
# 5. LOAD RAW DATA
# =============================================================================

cat("\n============================================================\n")
cat("0A GOLDBERG — Loading raw data\n")
cat("============================================================\n")

# Diet data — energy intake
diet_hcc <- readr::read_csv(
  file.path(DATA_DIR, "Diet_HCC.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    study_id    = normalise_id(.data[[DIET_SPEC$id_col]]),
    case_status = "HCC"
  )

diet_cirr <- readr::read_csv(
  file.path(DATA_DIR, "Diet_Cirr.csv"),
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    study_id    = normalise_id(.data[[DIET_SPEC$id_col]]),
    case_status = "Cirrhosis"
  )

# Clinical data — age, BMI, sex
clin_hcc <- readr::read_csv(
  CLINICAL_SPEC$hcc_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    study_id = normalise_id(.data[[CLINICAL_SPEC$id_col]])
  )

clin_cirr <- readr::read_csv(
  CLINICAL_SPEC$cirrhosis_file,
  show_col_types = FALSE
) %>%
  dplyr::mutate(
    study_id = normalise_id(.data[[CLINICAL_SPEC$id_col]])
  )

cat("Diet HCC rows:", nrow(diet_hcc), "\n")
cat("Diet Cirrhosis rows:", nrow(diet_cirr), "\n")
cat("Clinical HCC rows:", nrow(clin_hcc), "\n")
cat("Clinical Cirrhosis rows:", nrow(clin_cirr), "\n")

# =============================================================================
# 6. MERGE AND COMPUTE EI:BMR
# =============================================================================

# Extract energy intake column (kJ, convert to kcal)
KJ_TO_KCAL <- 1 / 4.184

diet_all <- dplyr::bind_rows(diet_hcc, diet_cirr) %>%
  dplyr::mutate(
    energy_kj = as_numeric_clean(.data[[DIET_SPEC$energy_col]]),
    energy_kcal = energy_kj * KJ_TO_KCAL
  ) %>%
  dplyr::select(study_id, case_status, energy_kcal)

# Detect sex column name in clinical data (could be "Sex", "Gender", etc.)
detect_sex_col <- function(df) {
  candidates <- c("Sex", "sex", "Gender", "gender", "SEX", "GENDER")
  found <- intersect(candidates, colnames(df))
  if (length(found) == 0L) {
    warning(
      "No sex/gender column found in clinical data. ",
      "Available columns: ", paste(head(colnames(df), 20), collapse = ", "),
      call. = FALSE
    )
    return(NULL)
  }
  found[1]
}

sex_col_hcc  <- detect_sex_col(clin_hcc)
sex_col_cirr <- detect_sex_col(clin_cirr)

extract_clinical <- function(df, id_col, age_col, bmi_col, sex_col) {
  out <- df %>%
    dplyr::transmute(
      study_id = normalise_id(.data[[id_col]]),
      age      = as_numeric_clean(.data[[age_col]]),
      bmi      = as_numeric_clean(.data[[bmi_col]]),
      sex      = if (!is.null(sex_col)) {
        as.character(.data[[sex_col]])
      } else {
        NA_character_
      }
    )
  out
}

clin_all <- dplyr::bind_rows(
  extract_clinical(
    clin_hcc,
    CLINICAL_SPEC$id_col,
    CLINICAL_SPEC$age_col,
    CLINICAL_SPEC$bmi_col,
    sex_col_hcc
  ),
  extract_clinical(
    clin_cirr,
    CLINICAL_SPEC$id_col,
    CLINICAL_SPEC$age_col,
    CLINICAL_SPEC$bmi_col,
    sex_col_cirr
  )
)

# Merge diet + clinical
goldberg_data <- diet_all %>%
  dplyr::left_join(clin_all, by = "study_id") %>%
  dplyr::distinct(study_id, .keep_all = TRUE)

cat("\nMerged dataset rows:", nrow(goldberg_data), "\n")

# =============================================================================
# 7. COMPUTE BMR AND CLASSIFY
# =============================================================================

goldberg_data <- goldberg_data %>%
  dplyr::rowwise() %>%
  dplyr::mutate(
    has_energy = !is.na(energy_kcal) & is.finite(energy_kcal),
    has_age    = !is.na(age) & is.finite(age),
    has_bmi    = !is.na(bmi) & is.finite(bmi),
    has_sex    = !is.na(sex) & nchar(trimws(sex)) > 0,
    data_issues = make_issue_string(has_energy, has_age, has_bmi, has_sex),
    bmr_kcal = if (has_age & has_bmi & has_sex) {
      estimate_bmr_mifflin(age, bmi, sex)
    } else {
      NA_real_
    },
    ei_bmr_ratio = if (has_energy & !is.na(bmr_kcal) & bmr_kcal > 0) {
      energy_kcal / bmr_kcal
    } else {
      NA_real_
    },
    reporter_class = classify_goldberg(energy_kcal, bmr_kcal)
  ) %>%
  dplyr::ungroup()

# =============================================================================
# 8. SUMMARY
# =============================================================================

cat("\n============================================================\n")
cat("GOLDBERG CLASSIFICATION SUMMARY\n")
cat("============================================================\n")

class_summary <- goldberg_data %>%
  dplyr::count(reporter_class, name = "n") %>%
  dplyr::mutate(pct = round(100 * n / sum(n), 1))

print(as.data.frame(class_summary))

cat("\nBy case status:\n")
cross_summary <- goldberg_data %>%
  dplyr::count(case_status, reporter_class, name = "n") %>%
  tidyr::pivot_wider(names_from = reporter_class, values_from = n, values_fill = 0L)

print(as.data.frame(cross_summary))

# Participants with data issues
n_issues <- sum(!is.na(goldberg_data$data_issues))
if (n_issues > 0L) {
  cat("\nParticipants with data issues:", n_issues, "\n")
}

# =============================================================================
# 9. SAVE OUTPUTS
# =============================================================================

# Primary output: classification file consumed by 0_preprocess.R
# This goes in DATA_DIR so the preprocessing script can find it via
# GOLDBERG_SPEC$file in config.R.
output_classification <- goldberg_data %>%
  dplyr::select(study_id, reporter_class, ei_bmr_ratio, data_issues)

readr::write_csv(
  output_classification,
  GOLDBERG_SPEC$file
)

cat(
  "\nClassification written to:",
  GOLDBERG_SPEC$file,
  "\n"
)

# Detailed summary for results
readr::write_csv(
  goldberg_data %>%
    dplyr::select(
      study_id, case_status, energy_kcal, age, bmi, sex,
      bmr_kcal, ei_bmr_ratio, reporter_class, data_issues
    ),
  file.path(GOLDBERG_RESULTS_DIR, "goldberg_summary.csv")
)

cat(
  "Detailed summary written to:",
  file.path(GOLDBERG_RESULTS_DIR, "goldberg_summary.csv"),
  "\n"
)

cat("\n============================================================\n")
cat("0A GOLDBERG CLASSIFICATION COMPLETE\n")
cat("============================================================\n")