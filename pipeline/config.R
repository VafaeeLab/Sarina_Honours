# =============================================================================
# config.R
# Shared configuration for the thesis analysis pipeline.
#
# SETTINGS ONLY: paths, cohort definitions and analysis parameters.
# =============================================================================

# -----------------------------------------------------------------------------
# 1. Reproducibility
# -----------------------------------------------------------------------------

MASTER_SEED <- 42L

# -----------------------------------------------------------------------------
# 2. Project directories
# -----------------------------------------------------------------------------

PROJECT_ROOT <- here::here()

DATA_DIR <- file.path(PROJECT_ROOT, "data")

DIRS <- list(
  intermediate = file.path(PROJECT_ROOT, "intermediate"),
  results      = file.path(PROJECT_ROOT, "results"),
  figures      = file.path(PROJECT_ROOT, "figures")
)

RDS_PATHS <- list(
  preprocess = file.path(DIRS$intermediate, "00_preprocessed_diet.rds"),
  clustering = file.path(DIRS$intermediate, "01_clustering.rds"),
  stage1     = file.path(DIRS$intermediate, "02_stage1_screen.rds"),
  biomarker  = file.path(DIRS$intermediate, "03_biomarker.rds"),
  deep_dive  = file.path(DIRS$intermediate, "04_deep_dive.rds")
)

# -----------------------------------------------------------------------------
# 3. Raw-data paths and schemas
# ----------------------------------------------------------------------------


DIET_SPEC <- list(
  file = file.path(DATA_DIR, "Diet_Data_MAFLD_Only.xlsx"),
  sheets = c(
    Cirrhosis = "Cirrhosis",
    HCC       = "HCC"
  ),
  id_col     = "Study ID",
  energy_col = "kjwithdf"
)

GOLDBERG_SPEC <- list(
  file               = file.path(DATA_DIR, "goldberg_reporter_classification.csv"),
  id_col             = "study_id",
  reporter_class_col = "reporter_class",
  ler_pattern        = "low energy reporter|\\(ler\\)"
)

CLINICAL_SPEC <- list(
  cirrhosis_file = file.path(DATA_DIR, "Clin_Path_MAFLD_Only(Cirrhosis).csv"),
  hcc_file       = file.path(DATA_DIR, "Clin_Path_MAFLD_Only(HCC).csv"),
  id_col         = "Study ID",
  age_col        = "Age at Baseline Visit",
  bmi_col        = "BMI"
)

OMICS_PANELS <- list(
  cytokines = list(
    label             = "Cytokines",
    cirrhosis_file    = file.path(DATA_DIR, "Cytokines_Cirrhosis.csv"),
    hcc_file          = file.path(DATA_DIR, "Cytokines_HCC.csv"),
    cirrhosis_id_col  = "SampleID",
    hcc_id_col        = "SampleID",
    visit_col         = "Visit #"
  ),
  small_molecules = list(
    label             = "Small Molecules",
    cirrhosis_file    = file.path(DATA_DIR, "SmallMolecules_Cirrhosis.csv"),
    hcc_file          = file.path(DATA_DIR, "SmallMolecules_HCC.csv"),
    cirrhosis_id_col  = "Study ID",
    hcc_id_col        = "Study ID",
    visit_col         = "Visit #"
  ),
  lipoproteins = list(
    label             = "Lipoproteins",
    cirrhosis_file    = file.path(DATA_DIR, "Lipoproteins(Cirrhosis).csv"),
    hcc_file          = file.path(DATA_DIR, "Lipoproteins(HCC).csv"),
    cirrhosis_id_col  = "sourceID",
    hcc_id_col        = "Study ID",
    visit_col         = "Visit #"
  )
)

# -----------------------------------------------------------------------------
# 4. Cohort definitions
# -----------------------------------------------------------------------------

# Primary diet analysis excludes low-energy reporters using the Goldberg file.
EXCLUDE_LER <- TRUE

# Omics analyses use baseline / Visit 1 measurements.
OMICS_VISIT <- 1L

# Biomarker-analysis condition labels retained from the current analysis.
BIOMARKER_CONDITIONS <- c(
  primary     = "A_excluded",
  sensitivity = "B_full",
  reference   = "C_full_cohort"
)

# -----------------------------------------------------------------------------
# 5. Diet preprocessing parameters
# -----------------------------------------------------------------------------

# Keep the current definition: each food-group score is the MEAN daily
# frequency across its constituent FFQ items.
FOOD_GROUP_AGGREGATION <- "mean"

# Remaining missing FFQ item frequencies are median-imputed before food-group
# aggregation, matching the current preprocessing pipeline.
FFQ_IMPUTATION <- "median"

# -----------------------------------------------------------------------------
# 6. Omics preprocessing / Stage 1 parameters
# -----------------------------------------------------------------------------

# Force log1p for cytokines ONLY. Small molecules and lipoproteins remain on
# their original scale unless this vector is deliberately changed later.
FORCE_LOG_PANELS <- c("cytokines")

# Current missing-value handling for retained omics features.
OMICS_IMPUTATION <- "median"

# Stage 1 feature-screening / PCA / Hotelling T-squared settings.
STAGE1_MAX_FEATURE_MISSING  <- 0.30
STAGE1_PA_ITER              <- 1000L
STAGE1_HOTELLING_CONFIDENCE <- 0.99

# -----------------------------------------------------------------------------
# 7. Clustering parameters
# -----------------------------------------------------------------------------

# Outcome-independent clustering audit. Each method is tuned with its own
# internal criterion; the frozen candidates are then compared using a common
# repeated-subsample Jaccard stability framework.
CLUSTER_K_MAX                  <- 5L
CLUSTER_PC_MAX_CAP             <- 8L
CLUSTER_PA_ITER                <- 1000L
CLUSTER_GAP_B                  <- 200L
CLUSTER_TUNE_RESAMPLES         <- 40L
CLUSTER_FINAL_RESAMPLES        <- 200L
CLUSTER_SUBSAMPLE_FRACTION     <- 0.80
CLUSTER_FCM_M_VALUES           <- c(1.5, 2.0, 2.5)
CLUSTER_UMAP_NEIGHBORS         <- c(10L, 15L, 25L)
CLUSTER_UMAP_MIN_DIST          <- c(0.05, 0.15, 0.30)

# Primary-method rule. Recommended thesis setting is "configured": complete
# the audit, preserve it, then freeze the pre-reviewed method below. Setting
# this to "stability" instead chooses the candidate with the highest minimum
# cluster Jaccard, then mean Jaccard, then resampling success rate.
CLUSTER_PRIMARY_SELECTION      <- "configured"
PRIMARY_CLUSTER_METHOD         <- "FCM"

# Expected frozen FCM result from the existing clustering audit. When FCM is
# the primary method, the pipeline stops if a rerun no longer selects these
# values, preventing an unnoticed change in the thesis clustering solution.
FCM_K                          <- 2L
FCM_M                          <- 1.5
FCM_N_STARTS                   <- 25L
CLUSTER_ASSERT_EXPECTED_FCM    <- TRUE

# -----------------------------------------------------------------------------
# 8. Biomarker-analysis parameters
# -----------------------------------------------------------------------------

FDR_THRESHOLD <- 0.10

# -----------------------------------------------------------------------------
# 9. Deep-dive parameters
# -----------------------------------------------------------------------------

DEEP_DIVE_PERMUTATIONS         <- 9999L
DEEP_DIVE_STANDARDIZE_DISTANCE <- TRUE
DEEP_DIVE_TOP_LOADINGS         <- 10L

# =============================================================================
# End config.R
# =============================================================================
