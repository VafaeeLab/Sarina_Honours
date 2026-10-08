# =============================================================================
# 00_preprocess.R
# Diet preprocessing for the thesis pipeline
#
# PURPOSE
#   Raw FFQ data -> 18 food-group scores -> Willett energy residualisation
#   -> z-standardised clustering matrix.
#
# INPUTS
#   Defined in config.R:
#     DIET_SPEC$file
#     DIET_SPEC$sheets
#     GOLDBERG_SPEC$file
#
# OUTPUT
#   RDS_PATHS$preprocess
#
# IMPORTANT METHOD DECISIONS
#   * Primary analysis excludes Goldberg low-energy reporters when EXCLUDE_LER=TRUE.
#   * Food-group scores are ROW MEANS of daily FFQ frequencies.
#   * Remaining FFQ item missingness is median-imputed before aggregation.
#   * Willett residual adjustment is applied before z-standardisation.
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

source("pipeline/config.R")
source("pipeline/utils.R")

set.seed(MASTER_SEED)

if (!identical(FOOD_GROUP_AGGREGATION, "mean")) {
  stop(
    "This pipeline is defined to use rowMeans() for food-group aggregation. ",
    "Set FOOD_GROUP_AGGREGATION <- 'mean' in config.R.",
    call. = FALSE
  )
}

if (!identical(FFQ_IMPUTATION, "median")) {
  stop(
    "This script currently implements median FFQ imputation only.",
    call. = FALSE
  )
}

assert_file_exists(DIET_SPEC$file, "diet workbook")
assert_file_exists(GOLDBERG_SPEC$file, "Goldberg classification file")
ensure_parent_dir(RDS_PATHS$preprocess)

# =============================================================================
# 1. LOAD DIET DATA
# =============================================================================

available_sheets <- readxl::excel_sheets(DIET_SPEC$file)
required_sheets  <- unname(DIET_SPEC$sheets)
missing_sheets   <- setdiff(required_sheets, available_sheets)

if (length(missing_sheets) > 0L) {
  stop(
    "Diet workbook is missing required sheet(s): ",
    paste(missing_sheets, collapse = ", "),
    call. = FALSE
  )
}

read_diet_sheet <- function(sheet_name, case_status) {
  dat <- readxl::read_excel(
    DIET_SPEC$file,
    sheet = sheet_name,
    .name_repair = "minimal"
  )
  
  assert_columns(
    dat,
    c(DIET_SPEC$id_col, DIET_SPEC$energy_col),
    paste0("diet sheet '", sheet_name, "'")
  )
  
  names(dat)[names(dat) == DIET_SPEC$id_col] <- "study_id"
  
  dat %>%
    dplyr::mutate(
      study_id = normalise_id(study_id),
      case_status = case_status
    )
}

diet_cirr <- read_diet_sheet(
  DIET_SPEC$sheets[["Cirrhosis"]],
  "Cirrhosis"
)

diet_hcc <- read_diet_sheet(
  DIET_SPEC$sheets[["HCC"]],
  "HCC"
)

assert_unique_ids(diet_cirr, "study_id", "diet Cirrhosis sheet")
assert_unique_ids(diet_hcc, "study_id", "diet HCC sheet")

cross_duplicates <- intersect(
  stats::na.omit(diet_cirr$study_id),
  stats::na.omit(diet_hcc$study_id)
)

if (length(cross_duplicates) > 0L) {
  stop(
    "Study IDs appear in both dietary sheets: ",
    paste(utils::head(cross_duplicates, 10L), collapse = ", "),
    call. = FALSE
  )
}

# =============================================================================
# 2. GOLDBERG LOW-ENERGY-REPORTER EXCLUSION
# =============================================================================

ler_ids <- character(0)

if (isTRUE(EXCLUDE_LER)) {
  ler_data <- readr::read_csv(
    GOLDBERG_SPEC$file,
    show_col_types = FALSE,
    name_repair = "minimal"
  )
  
  assert_columns(
    ler_data,
    c(
      GOLDBERG_SPEC$id_col,
      GOLDBERG_SPEC$reporter_class_col
    ),
    "Goldberg classification"
  )
  
  ler_ids <- ler_data %>%
    dplyr::transmute(
      study_id = normalise_id(.data[[GOLDBERG_SPEC$id_col]]),
      reporter_class = stringr::str_squish(
        as.character(.data[[GOLDBERG_SPEC$reporter_class_col]])
      )
    ) %>%
    dplyr::filter(
      stringr::str_detect(
        reporter_class,
        stringr::regex(
          GOLDBERG_SPEC$ler_pattern,
          ignore_case = TRUE
        )
      )
    ) %>%
    dplyr::pull(study_id) %>%
    unique()
  
  cat("Goldberg LER IDs:", length(ler_ids), "\n")
  cat(
    "  matched in Cirrhosis:",
    sum(diet_cirr$study_id %in% ler_ids),
    "\n"
  )
  cat(
    "  matched in HCC:",
    sum(diet_hcc$study_id %in% ler_ids),
    "\n"
  )
  
  diet_cirr <- diet_cirr %>%
    dplyr::filter(!study_id %in% ler_ids)
  
  diet_hcc <- diet_hcc %>%
    dplyr::filter(!study_id %in% ler_ids)
} else {
  cat("Goldberg LER exclusion: OFF\n")
}

# Align the two sheets before row-binding.
all_cols <- union(names(diet_cirr), names(diet_hcc))

for (col in setdiff(all_cols, names(diet_cirr))) {
  diet_cirr[[col]] <- NA
}

for (col in setdiff(all_cols, names(diet_hcc))) {
  diet_hcc[[col]] <- NA
}

diet_raw <- dplyr::bind_rows(
  diet_cirr %>% dplyr::select(dplyr::all_of(all_cols)),
  diet_hcc  %>% dplyr::select(dplyr::all_of(all_cols))
)

if (isTRUE(EXCLUDE_LER) && any(diet_raw$study_id %in% ler_ids)) {
  stop("LER exclusion failed: one or more LER participants remain.", call. = FALSE)
}

if (anyNA(diet_raw$study_id)) {
  stop("Missing study_id values remain in the analytic diet cohort.", call. = FALSE)
}

assert_unique_ids(diet_raw, "study_id", "combined diet cohort")

cat("\n=== ANALYTIC DIET COHORT ===\n")
cat("Participants:", nrow(diet_raw), "\n")
print(table(diet_raw$case_status, useNA = "ifany"))

# =============================================================================
# 3. FOOD-GROUP DICTIONARY
# =============================================================================

food_group_items <- list(
  veg = c(
    "Potato", "Pumpkin", "Sweet potato", "Cauliflower",
    "Green beans", "Spinach", "Cabbage", "Peas", "Broccoli",
    "Carrots", "Zucchini", "Capsicum", "Corn", "Mushrooms",
    "Tomatoes", "Lettuce", "Celery, cucumber", "Avocado",
    "Onion, spring onion etc."
  ),
  legumes = c(
    "Lentils", "Kidney/black beans", "Chickpeas", "Butter beans",
    "3/4 bean mix", "Cannelli beans", "Baked beans",
    "Non-soy veg meat replacements", "Non-soy veg products"
  ),
  fruit = c(
    "Canned fruit", "Fruit salad", "Dried fruit", "Apple or pear",
    "Orange", "Banana", "Peach, nectarine etc", "mango or paw-paw",
    "Pineapple", "Grapes, berries", "Melon"
  ),
  grains_cereals = c(
    "Muesli", "Rolled oats", "Breakfast cereal",
    "bread/roll/pita/toast", "English muffin/bagel/crumpet",
    "Rice", "Other grains", "Noodles", "Pasta"
  ),
  dairy_and_alternatives = c(
    "Flavoured milk", "Plain milk", "Cream or sour cream",
    "Frozen yoghurt", "Yoghurt (not frozen)",
    "Cottage cheese or ricotta", "Cheese", "Cheese spread",
    "Soy yogurt", "Soy cheese"
  ),
  red_meat = c(
    "Mince dish", "Beef/lamb + sauce, no vege",
    "Beef/lamb + sauce, w vege", "Beef/lamp, no vege",
    "Beef/lamb, w vege", "Pork + sauce, no vege",
    "Pork + sauce, w vege", "Pork, no vege", "Pork, w vege",
    "Liver"
  ),
  processed_meat = c(
    "Devon, salami", "Bacon, ham", "Sausages", "Hot dog",
    "Pie, sausage roll", "Savoury pastries", "Hamburgers"
  ),
  poultry = c(
    "Chicken + sauce, no vege", "Chicken + sauce, w vege",
    "Crumbed chicken", "Chicken, no vege", "Chicken, w vege"
  ),
  seafood = c(
    "Canned fish", "Oily fish", "White fish",
    "Crustaceans", "Molluscs", "Squid, calamari"
  ),
  eggs = c("Eggs"),
  nuts = c(
    "Almonds", "Walnuts", "Pecans", "Hazelnuts", "Cashews",
    "Pistachios", "Macadamias", "Pine nuts", "Mixed nuts",
    "Peanut butter, nutella"
  ),
  discretionary_sweet = c(
    "Cakes, muffins etc", "Sweet pies", "Other puddings",
    "Plain sweet biscuits", "Cream or chocolate biscuits",
    "Fruit bars", "Snack bars", "Muesli bars", "Chocolate",
    "Lollies, no chocolate", "Iceblock - creamy",
    "Iceblock - water", "soft drinks", "Cordial",
    "fruit-based juice", "Ice cream", "Sweet combination snacks",
    "Diet soft drink", "Jam, honey, syrup", "Jelly"
  ),
  discretionary_savoury = c(
    "Dry or savoury biscuits", "Savoury combination snacks",
    "Snack noodles", "Pizza",
    "Tacos, burritos, etc", "Hot chips", "Hot chips at home",
    "Hash browns", "Chips (not potato)", "Crisps", "Vegemite"
  ),
  alcohol = c("Beer", "Wine or cooler", "Spirits"),
  saturated = c("Butter","Margarine", "Shortening", "Mixed or blended oil"),
  unsaturated = c(
    "Sesame oil", "Oil sprays", "Olive oil, etc", "Sunflower oil",
    "Low fat salad dressing", "Salad dressing, mayonnaise"
  ),
  healthy_drinks = c("water", "Tea or coffee"),
  soy_plant = c(
    "Tempeh", "Silken/firm tofu", "Flavoured tofu",
    "Soy sausage", "Soy deli slices", "Soy frozen foods"
  )
)

food_group_vars <- names(food_group_items)

if (length(food_group_vars) != 18L) {
  stop(
    "Expected 18 food groups but found ",
    length(food_group_vars),
    ".",
    call. = FALSE
  )
}

all_items <- unlist(food_group_items, use.names = FALSE)
duplicated_items <- unique(all_items[duplicated(all_items)])

if (length(duplicated_items) > 0L) {
  stop(
    "FFQ items are assigned to more than one food group: ",
    paste(duplicated_items, collapse = ", "),
    call. = FALSE
  )
}

expected_items <- unique(all_items)
present_items  <- intersect(expected_items, names(diet_raw))
missing_items  <- setdiff(expected_items, names(diet_raw))

cat("\n=== FOOD-GROUP DICTIONARY COVERAGE ===\n")
cat("Expected FFQ items:", length(expected_items), "\n")
cat("Present:", length(present_items), "\n")
cat("Missing:", length(missing_items), "\n")

if (length(missing_items) > 0L) {
  warning(
    "Missing FFQ item column(s): ",
    paste(missing_items, collapse = ", "),
    call. = FALSE
  )
}

for (group_name in food_group_vars) {
  found <- intersect(food_group_items[[group_name]], names(diet_raw))
  
  if (length(found) == 0L) {
    stop(
      "No usable FFQ items were found for food group '",
      group_name,
      "'.",
      call. = FALSE
    )
  }
}

# =============================================================================
# 4. FFQ RESPONSE -> DAILY FREQUENCY
# =============================================================================

frequency_to_daily <- c(
  "never"                 = 0,
  "less than 1 per month" = 0.5 / 30.4375,
  "1-3 per month"         = 2 / 30.4375,
  "less than 1 per week"  = 0.5 / 7,
  "1 per week or less"    = 0.75 / 7,
  "1 per week"            = 1 / 7,
  "once per week"         = 1 / 7,
  "1-2 per week"          = 1.5 / 7,
  "2 or more per week"    = 2.5 / 7,
  "2-4 per week"          = 3 / 7,
  "3-4 per week"          = 3.5 / 7,
  "2-6 per week"          = 4 / 7,
  "5-6 per week"          = 5.5 / 7,
  "5 or more per week"    = 5.5 / 7,
  "5-7 per week"          = 6 / 7,
  "once a day"            = 1,
  "every day"             = 1,
  "1 per day"             = 1,
  "1 or more per day"     = 1,
  "1-2 per day"           = 1.5,
  "1-3 per day"           = 2,
  "2 or more per day"     = 2,
  "2-3 per day"           = 2.5,
  "4 or more per day"     = 4,
  "4-6 per day"           = 5,
  "less than 1 per day"   = 0.5,
  "4-6 glasses/day"       = 5,
  "7+ glasses/day"        = 7
)

map_to_daily <- function(x) {
  cleaned <- normalise_response(x)
  unname(frequency_to_daily[cleaned])
}

unmapped_found <- FALSE

for (item in present_items) {
  cleaned <- normalise_response(diet_raw[[item]])
  bad <- !is.na(cleaned) & !(cleaned %in% names(frequency_to_daily))
  
  if (any(bad)) {
    unmapped_values <- sort(unique(cleaned[bad]))
    
    warning(
      "Unmapped FFQ response(s) in '",
      item,
      "': ",
      paste(unmapped_values, collapse = ", "),
      call. = FALSE
    )
    
    unmapped_found <- TRUE
  }
}

if (unmapped_found) {
  warning(
    "Unmapped FFQ responses will become NA before median imputation. ",
    "Review the preceding warnings.",
    call. = FALSE
  )
}

ffq_numeric <- diet_raw %>%
  dplyr::transmute(
    dplyr::across(
      dplyr::all_of(present_items),
      map_to_daily
    )
  )

item_missingness <- tibble::tibble(
  item = names(ffq_numeric),
  missing_fraction = vapply(
    ffq_numeric,
    function(x) mean(is.na(x)),
    numeric(1)
  )
)

high_missing_items <- item_missingness %>%
  dplyr::filter(missing_fraction > 0.30)

if (nrow(high_missing_items) > 0L) {
  warning(
    nrow(high_missing_items),
    " FFQ item(s) have >30% missingness after conversion: ",
    paste(high_missing_items$item, collapse = ", "),
    call. = FALSE
  )
}

cat("\n=== FFQ CONVERSION ===\n")
cat("Items converted:", ncol(ffq_numeric), "\n")
cat("NAs before imputation:", sum(is.na(ffq_numeric)), "\n")

ffq_numeric <- median_impute_numeric(
  ffq_numeric,
  cols = names(ffq_numeric),
  label = "FFQ daily-frequency data"
)

if (anyNA(ffq_numeric)) {
  stop("NAs remain after FFQ median imputation.", call. = FALSE)
}

cat("NAs after imputation:", sum(is.na(ffq_numeric)), "\n")

# =============================================================================
# 5. AGGREGATE INTO 18 FOOD-GROUP SCORES USING ROW MEANS
# =============================================================================

food_group_columns <- lapply(
  food_group_vars,
  function(group_name) {
    safe_rowmean(
      ffq_numeric,
      food_group_items[[group_name]]
    )
  }
)

names(food_group_columns) <- food_group_vars

diet_groups <- tibble::as_tibble(food_group_columns) %>%
  dplyr::mutate(
    id = diet_raw$study_id,
    case_status = diet_raw$case_status,
    .before = 1
  )

if (nrow(diet_groups) != nrow(diet_raw)) {
  stop("Row mismatch after food-group aggregation.", call. = FALSE)
}

if (any(!is.finite(as.matrix(diet_groups[, food_group_vars, drop = FALSE])))) {
  stop("Non-finite food-group scores were created.", call. = FALSE)
}

cat("\n=== FOOD GROUPS ===\n")
cat("Food groups:", length(food_group_vars), "\n")
cat("Aggregation: rowMeans\n")

# =============================================================================
# 6. TOTAL ENERGY
# =============================================================================

energy_raw <- suppressWarnings(
  as.numeric(diet_raw[[DIET_SPEC$energy_col]])
)

diet_groups$total_energy_kj   <- energy_raw
diet_groups$total_energy_kcal <- energy_raw / 4.184

valid_energy <- (
  !is.na(diet_groups$total_energy_kcal) &
    is.finite(diet_groups$total_energy_kcal) &
    diet_groups$total_energy_kcal > 0
)

invalid_energy_ids <- diet_groups$id[!valid_energy]

if (length(invalid_energy_ids) > 0L) {
  warning(
    length(invalid_energy_ids),
    " participant(s) have invalid total energy and will be excluded before ",
    "Willett adjustment: ",
    paste(utils::head(invalid_energy_ids, 10L), collapse = ", "),
    call. = FALSE
  )
}

diet_for_willett <- diet_groups %>%
  dplyr::filter(valid_energy)

if (nrow(diet_for_willett) < 10L) {
  stop(
    "Fewer than 10 participants have valid energy values.",
    call. = FALSE
  )
}

cat("\n=== ENERGY ===\n")
cat("Valid energy records:", nrow(diet_for_willett), "\n")
print(
  diet_for_willett %>%
    dplyr::group_by(case_status) %>%
    dplyr::summarise(
      n = dplyr::n(),
      mean_kcal = round(mean(total_energy_kcal)),
      sd_kcal = round(stats::sd(total_energy_kcal)),
      median_kcal = round(stats::median(total_energy_kcal)),
      .groups = "drop"
    )
)

# =============================================================================
# 7. WILLETT ENERGY RESIDUALISATION
# =============================================================================

mean_energy <- mean(diet_for_willett$total_energy_kcal)
adj_cols <- paste0(food_group_vars, "_adj")

residual_diagnostics <- tibble::tibble(
  food_group = food_group_vars,
  n_complete = NA_integer_,
  r_squared = NA_real_,
  slope_per_1000kcal = NA_real_,
  p_energy = NA_real_,
  predicted_at_mean_energy = NA_real_
)

cat("\n=== WILLETT ENERGY RESIDUALISATION ===\n")

for (i in seq_along(food_group_vars)) {
  fg  <- food_group_vars[[i]]
  adj <- adj_cols[[i]]
  
  y      <- diet_for_willett[[fg]]
  energy <- diet_for_willett$total_energy_kcal
  
  complete <- is.finite(y) & is.finite(energy)
  n_ok <- sum(complete)
  
  if (n_ok < 10L) {
    stop(
      "Too few complete observations for Willett model: ",
      fg,
      " (n=",
      n_ok,
      ").",
      call. = FALSE
    )
  }
  
  model_data <- data.frame(
    outcome = y[complete],
    total_energy_kcal = energy[complete]
  )
  
  fit <- stats::lm(
    outcome ~ total_energy_kcal,
    data = model_data
  )
  
  fit_summary <- summary(fit)
  
  pred_at_mean <- as.numeric(
    stats::predict(
      fit,
      newdata = data.frame(total_energy_kcal = mean_energy)
    )
  )
  
  adjusted_complete <- stats::residuals(fit) + pred_at_mean
  
  adj_vec <- rep(NA_real_, nrow(diet_for_willett))
  adj_vec[complete] <- adjusted_complete
  diet_for_willett[[adj]] <- adj_vec
  
  residual_diagnostics$n_complete[[i]] <- n_ok
  residual_diagnostics$r_squared[[i]] <- fit_summary$r.squared
  residual_diagnostics$slope_per_1000kcal[[i]] <- stats::coef(fit)[[2]] * 1000
  residual_diagnostics$p_energy[[i]] <- stats::coef(fit_summary)[2, 4]
  residual_diagnostics$predicted_at_mean_energy[[i]] <- pred_at_mean
}

adj_cols_present <- adj_cols[
  vapply(
    diet_for_willett[adj_cols],
    function(x) any(is.finite(x)),
    logical(1)
  )
]

if (length(adj_cols_present) != length(food_group_vars)) {
  stop(
    "Not all 18 adjusted food groups were successfully created.",
    call. = FALSE
  )
}

diet_final <- diet_for_willett

if (anyNA(diet_final[, adj_cols_present, drop = FALSE])) {
  stop(
    "Adjusted food-group values contain NA after Willett residualisation. ",
    "Review upstream FFQ/energy processing rather than imputing at this stage.",
    call. = FALSE
  )
}

# =============================================================================
# 8. Z-STANDARDISE THE FROZEN CLUSTERING MATRIX
# =============================================================================

adj_matrix <- diet_final %>%
  dplyr::select(dplyr::all_of(adj_cols_present))

adj_sds <- vapply(adj_matrix, stats::sd, numeric(1))
zero_var <- names(adj_sds)[
  !is.finite(adj_sds) |
    adj_sds <= sqrt(.Machine$double.eps)
]

if (length(zero_var) > 0L) {
  stop(
    "Zero-variance adjusted food-group column(s): ",
    paste(zero_var, collapse = ", "),
    call. = FALSE
  )
}

adj_scaled <- scale(adj_matrix)
colnames(adj_scaled) <- adj_cols_present
rownames(adj_scaled) <- diet_final$id

scaled_means <- colMeans(adj_scaled)
scaled_sds   <- apply(adj_scaled, 2, stats::sd)

if (max(abs(scaled_means)) > 1e-10) {
  stop("Scaled food-group means are not approximately zero.", call. = FALSE)
}

if (max(abs(scaled_sds - 1)) > 1e-10) {
  stop("Scaled food-group SDs are not approximately one.", call. = FALSE)
}

cat("\n=== FINAL PREPROCESSED DATA ===\n")
cat("Participants:", nrow(diet_final), "\n")
cat("Food groups:", ncol(adj_scaled), "\n")
print(table(diet_final$case_status, useNA = "ifany"))

# =============================================================================
# 9. SAVE STAGE OUTPUT
# =============================================================================

preprocess_output <- list(
  diet_final = diet_final,
  adj_scaled = adj_scaled,
  adj_cols = adj_cols_present,
  food_group_vars = food_group_vars,
  food_group_items = food_group_items,
  item_missingness = item_missingness,
  residual_diagnostics = residual_diagnostics,
  cohort = list(
    participant_ids = as.character(diet_final$id),
    excluded_ler_ids = as.character(ler_ids),
    invalid_energy_ids = as.character(invalid_energy_ids)
  ),
  settings = list(
    seed = MASTER_SEED,
    exclude_ler = EXCLUDE_LER,
    food_group_aggregation = FOOD_GROUP_AGGREGATION,
    ffq_imputation = FFQ_IMPUTATION,
    diet_file = DIET_SPEC$file,
    diet_sheets = DIET_SPEC$sheets,
    goldberg_file = GOLDBERG_SPEC$file
  )
)

save_rds_checked(
  preprocess_output,
  RDS_PATHS$preprocess
)

cat("\n============================================================\n")
cat("PREPROCESSING COMPLETE\n")
cat("Saved:", RDS_PATHS$preprocess, "\n")
cat("============================================================\n")