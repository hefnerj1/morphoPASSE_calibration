# ============================================================================
# MorphoPASSE posterior-probability calibration: revised analysis
#
# This script:
#   1. imports the revised simulated and known-sex MorphoPASSE data;
#   2. removes erroneous duplicate known-sex records after verifying that
#      duplicate rows are identical;
#   3. merges the balanced random labels assigned after trait simulation;
#   4. estimates random-label classification benchmarks using out-of-bag
#      Random Forest predictions;
#   5. produces the manuscript tables and figures needed for revision;
#   6. saves cleaned data, permutation results, and session information.
#
# Default inputs (edit these paths if needed):
#   raw_Data_morpho.xlsx
#   simulated_combined_analysis_50.csv
#
# Optional command-line use:
#   Rscript morphopasse_revision_analysis.R \
#     raw_Data_morpho.xlsx simulated_combined_analysis_50.csv output_directory
# ============================================================================

# ------------------------------
# 0. Packages and settings
# ------------------------------
required_packages <- c(
  "readxl", "dplyr", "tidyr", "purrr", "stringr",
  "randomForest", "ggplot2", "patchwork", "scales"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]

if (length(missing_packages) > 0) {
  stop(
    "Install the following packages before running this script: ",
    paste(missing_packages, collapse = ", ")
  )
}

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(randomForest)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
input_workbook <- if (length(args) >= 1) args[1] else "raw_Data_morpho.xlsx"
random_label_file <- if (length(args) >= 2) args[2] else "simulated_combined_analysis_50.csv"
output_dir <- if (length(args) >= 3) args[3] else "morphopasse_revision_outputs"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 1234
n_perm <- 500
n_trees_perm <- 500
n_trees_independent <- 2000

set.seed(seed)

cranial_traits <- c("G", "NC", "SOL", "SOR", "MPL", "MPR")
pelvic_traits <- c("MAL", "MAR", "VAL", "VAR", "SPCL", "SPCR")
all_traits <- c(cranial_traits, pelvic_traits)

# ------------------------------
# 1. Utility functions
# ------------------------------
assert_columns <- function(df, required, object_name = "data frame") {
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop(
      object_name, " is missing required columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

standardize_sex <- function(x) {
  x <- toupper(trimws(as.character(x)))
  x[x %in% c("FEMALE", "F")] <- "F"
  x[x %in% c("MALE", "M")] <- "M"
  x
}

balanced_accuracy <- function(truth, prediction) {
  truth <- factor(truth, levels = c("F", "M"))
  prediction <- factor(prediction, levels = c("F", "M"))
  recalls <- vapply(levels(truth), function(level) {
    idx <- truth == level
    if (!any(idx, na.rm = TRUE)) return(NA_real_)
    mean(prediction[idx] == truth[idx], na.rm = TRUE)
  }, numeric(1))
  mean(recalls, na.rm = TRUE)
}

format_median_range <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) return(NA_character_)
  med <- median(x)
  med_txt <- if (abs(med - round(med)) < .Machine$double.eps^0.5) {
    sprintf("%.0f", med)
  } else {
    sprintf("%.1f", med)
  }
  sprintf("%s (%d-%d)", med_txt, min(x), max(x))
}

verify_and_deduplicate_known <- function(df, module_name) {
  assert_columns(df, "Ind_ID", module_name)

  duplicated_ids <- df %>%
    count(Ind_ID, name = "row_count") %>%
    filter(row_count > 1)

  if (nrow(duplicated_ids) == 0) return(df)

  comparison_cols <- setdiff(names(df), "Ind_ID")

  conflicts <- df %>%
    group_by(Ind_ID) %>%
    summarise(
      across(all_of(comparison_cols), ~ n_distinct(.x, na.rm = FALSE)),
      .groups = "drop"
    ) %>%
    filter(if_any(all_of(comparison_cols), ~ .x > 1))

  if (nrow(conflicts) > 0) {
    stop(
      "Duplicate IDs in ", module_name,
      " contain conflicting values. Inspect IDs: ",
      paste(conflicts$Ind_ID, collapse = ", ")
    )
  }

  message(
    module_name, ": removed ", nrow(df) - n_distinct(df$Ind_ID),
    " erroneous duplicate rows; retained ", n_distinct(df$Ind_ID),
    " unique individuals."
  )

  df %>% distinct(Ind_ID, .keep_all = TRUE)
}

read_morphopasse_sheet <- function(sheet_name, dataset_name, region_name, traits) {
  df <- read_excel(input_workbook, sheet = sheet_name)

  if ("ID" %in% names(df) && !"Ind_ID" %in% names(df)) {
    df <- df %>% rename(Ind_ID = ID)
  }

  required <- c("Ind_ID", traits, "CasePrediction", "postprob")
  if (dataset_name == "known") required <- c(required, "Sex")
  assert_columns(df, required, paste0("sheet '", sheet_name, "'"))

  df <- df %>%
    mutate(
      Ind_ID = as.character(Ind_ID),
      CasePrediction = standardize_sex(CasePrediction),
      postprob = as.numeric(postprob),
      dataset = dataset_name,
      region = region_name
    )

  for (trait in traits) df[[trait]] <- as.integer(df[[trait]])

  if (dataset_name == "known") {
    df <- df %>% mutate(Sex = standardize_sex(Sex))
    df <- verify_and_deduplicate_known(df, sheet_name)
  } else {
    df$Sex <- NA_character_
  }

  if (any(!df$CasePrediction %in% c("F", "M"))) {
    stop("Unexpected CasePrediction labels in sheet '", sheet_name, "'.")
  }
  if (dataset_name == "known" && any(!df$Sex %in% c("F", "M"))) {
    stop("Unexpected documented Sex labels in sheet '", sheet_name, "'.")
  }
  if (any(df$postprob < 0 | df$postprob > 1, na.rm = TRUE)) {
    stop("Posterior probabilities outside [0,1] in sheet '", sheet_name, "'.")
  }

  df
}

rf_random_label_benchmark <- function(
    df, traits, dataset_name, region_name,
    n_perm = 500, ntree = 500, seed = 1234) {

  X <- df %>% dplyr::select(all_of(traits))
  if (anyNA(X)) stop("Missing trait data in ", dataset_name, " ", region_name, ".")

  n <- nrow(X)
  if (n %% 2 != 0) {
    stop("Balanced random-label assignment requires an even sample size; n = ", n)
  }

  base_labels <- rep(c("F", "M"), each = n / 2)
  out <- vector("list", n_perm)

  set.seed(seed)
  for (i in seq_len(n_perm)) {
    y <- factor(sample(base_labels, size = n, replace = FALSE), levels = c("F", "M"))

    fit <- randomForest(
      x = X,
      y = y,
      ntree = ntree,
      mtry = max(1, floor(sqrt(length(traits)))),
      importance = FALSE,
      keep.forest = FALSE
    )

    pred <- fit$predicted
    valid <- !is.na(pred)

    out[[i]] <- tibble(
      dataset = dataset_name,
      region = region_name,
      permutation = i,
      oob_accuracy = mean(pred[valid] == y[valid]),
      oob_balanced_accuracy = balanced_accuracy(y[valid], pred[valid])
    )
  }

  bind_rows(out)
}

extract_vote_for_class <- function(votes, classes) {
  class_index <- match(classes, colnames(votes))
  if (anyNA(class_index)) {
    stop("A requested class is absent from the Random Forest vote matrix.")
  }
  votes[cbind(seq_len(nrow(votes)), class_index)]
}

fit_known_independent_rf <- function(df, traits, region_name, ntree, seed) {
  X <- df %>% dplyr::select(all_of(traits))
  y <- factor(df$Sex, levels = c("F", "M"))

  # Equal per-class sampling mirrors the equal-prior logic used in MorphoPASSE.
  class_sample_n <- min(table(y))
  class_sampsize <- setNames(rep(as.integer(class_sample_n), 2), levels(y))

  set.seed(seed)
  fit <- randomForest(
    x = X,
    y = y,
    ntree = ntree,
    mtry = max(1, floor(sqrt(length(traits)))),
    strata = y,
    sampsize = class_sampsize,
    importance = FALSE,
    keep.forest = TRUE
  )

  rf_pred <- as.character(fit$predicted)
  rf_votes_for_mp_class <- extract_vote_for_class(fit$votes, df$CasePrediction)
  rf_max_vote <- apply(fit$votes, 1, max, na.rm = TRUE)

  tibble(
    dataset = "known",
    region = region_name,
    n = nrow(df),
    morphopasse_accuracy = mean(df$CasePrediction == df$Sex),
    morphopasse_balanced_accuracy = balanced_accuracy(df$Sex, df$CasePrediction),
    independent_rf_oob_accuracy = mean(rf_pred == df$Sex),
    independent_rf_oob_balanced_accuracy = balanced_accuracy(df$Sex, rf_pred),
    prediction_agreement = mean(rf_pred == df$CasePrediction),
    spearman_rho_mp_postprob_vs_rf_vote_for_mp_class = suppressWarnings(
      cor(df$postprob, rf_votes_for_mp_class, method = "spearman", use = "complete.obs")
    ),
    spearman_rho_mp_postprob_vs_rf_max_vote = suppressWarnings(
      cor(df$postprob, rf_max_vote, method = "spearman", use = "complete.obs")
    )
  )
}

fit_simulated_internal_consistency <- function(df, traits, region_name, ntree, seed) {
  X <- df %>% dplyr::select(all_of(traits))
  y <- factor(df$CasePrediction, levels = c("F", "M"))

  class_sample_n <- min(table(y))
  class_sampsize <- setNames(rep(as.integer(class_sample_n), 2), levels(y))

  set.seed(seed)
  fit <- randomForest(
    x = X,
    y = y,
    ntree = ntree,
    mtry = max(1, floor(sqrt(length(traits)))),
    strata = y,
    sampsize = class_sampsize,
    importance = FALSE,
    keep.forest = TRUE
  )

  rf_pred <- as.character(fit$predicted)
  rf_vote_for_mp_class <- extract_vote_for_class(fit$votes, df$CasePrediction)

  tibble(
    dataset = "simulated",
    region = region_name,
    n = nrow(df),
    oob_accuracy_predicting_morphopasse_class = mean(rf_pred == df$CasePrediction),
    oob_balanced_accuracy_predicting_morphopasse_class = balanced_accuracy(
      df$CasePrediction, rf_pred
    ),
    spearman_rho_mp_postprob_vs_rf_vote_for_mp_class = suppressWarnings(
      cor(df$postprob, rf_vote_for_mp_class, method = "spearman", use = "complete.obs")
    )
  )
}

matrix_to_long <- function(mat, dataset_name, region_name, trait_order) {
  out <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(out) <- c("trait_1", "trait_2", "rho")
  out %>%
    mutate(
      dataset = dataset_name,
      region = region_name,
      trait_1 = factor(trait_1, levels = trait_order),
      trait_2 = factor(trait_2, levels = rev(trait_order)),
      i = match(as.character(trait_1), trait_order),
      j = match(as.character(trait_2), trait_order)
    ) %>%
    filter(i >= match(as.character(trait_2), trait_order))
}

# ------------------------------
# 2. Import and validate data
# ------------------------------
if (!file.exists(input_workbook)) stop("Input workbook not found: ", input_workbook)
if (!file.exists(random_label_file)) stop("Random-label file not found: ", random_label_file)

sim_cranial <- read_morphopasse_sheet(
  "simulated cranial data", "simulated", "cranial", cranial_traits
)
sim_pelvic <- read_morphopasse_sheet(
  "simulated pelvic data", "simulated", "pelvic", pelvic_traits
)
known_cranial <- read_morphopasse_sheet(
  "known cranial data", "known", "cranial", cranial_traits
)
known_pelvic <- read_morphopasse_sheet(
  "known pelvic data", "known", "pelvic", pelvic_traits
)

random_labels <- read.csv(
  random_label_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
) %>%
  mutate(
    Ind_ID = as.character(Ind_ID),
    RandomSex = standardize_sex(RandomSex)
  )

assert_columns(random_labels, c("Ind_ID", "RandomSex", all_traits), "random-label file")

if (anyDuplicated(random_labels$Ind_ID)) stop("Duplicate Ind_ID values in random-label file.")
if (any(!random_labels$RandomSex %in% c("F", "M"))) {
  stop("Unexpected RandomSex labels in random-label file.")
}
label_counts <- table(factor(random_labels$RandomSex, levels = c("F", "M")))
if (any(label_counts != c(F = 25L, M = 25L))) {
  warning(
    "RandomSex is not an exact 25/25 split; observed: ",
    paste(names(label_counts), as.integer(label_counts), collapse = ", ")
  )
}

# Confirm that the revised simulated traits in the MorphoPASSE workbook match
# the simulation/label file exactly.
sim_traits_from_workbook <- sim_cranial %>%
  dplyr::select(Ind_ID, all_of(cranial_traits)) %>%
  inner_join(
    sim_pelvic %>% dplyr::select(Ind_ID, all_of(pelvic_traits)),
    by = "Ind_ID"
  )

if (!setequal(sim_traits_from_workbook$Ind_ID, random_labels$Ind_ID)) {
  stop("Simulated IDs differ between the workbook and random-label file.")
}

validation <- sim_traits_from_workbook %>%
  inner_join(random_labels, by = "Ind_ID", suffix = c(".workbook", ".labels"))

for (trait in all_traits) {
  workbook_col <- paste0(trait, ".workbook")
  labels_col <- paste0(trait, ".labels")
  if (any(validation[[workbook_col]] != validation[[labels_col]], na.rm = TRUE)) {
    stop("Trait mismatch between workbook and random-label file for ", trait, ".")
  }
}

sim_cranial <- sim_cranial %>%
  left_join(random_labels %>% dplyr::select(Ind_ID, RandomSex), by = "Ind_ID")
sim_pelvic <- sim_pelvic %>%
  left_join(random_labels %>% dplyr::select(Ind_ID, RandomSex), by = "Ind_ID")
known_cranial$RandomSex <- NA_character_
known_pelvic$RandomSex <- NA_character_

# Bind the four modules. Traits not used in a given module remain NA.
dat_clean <- bind_rows(sim_cranial, sim_pelvic, known_cranial, known_pelvic) %>%
  arrange(dataset, region, Ind_ID)

write.csv(
  dat_clean,
  file.path(output_dir, "cleaned_all_modules.csv"),
  row.names = FALSE
)

write.csv(
  bind_rows(
    sim_cranial %>% mutate(module = "simulated_cranial"),
    sim_pelvic %>% mutate(module = "simulated_pelvic"),
    known_cranial %>% mutate(module = "known_cranial"),
    known_pelvic %>% mutate(module = "known_pelvic")
  ),
  file.path(output_dir, "cleaned_module_records.csv"),
  row.names = FALSE
)

# ------------------------------
# 3. Table 1: group summaries and requested trait ranges
# ------------------------------
make_group_summary <- function(df, traits, grouping_column, dataset_label, region_label) {
  df %>%
    group_by(Group = .data[[grouping_column]]) %>%
    summarise(
      n = n(),
      mean_postprob = mean(postprob, na.rm = TRUE),
      across(all_of(traits), format_median_range),
      .groups = "drop"
    ) %>%
    mutate(
      Dataset = dataset_label,
      Region = region_label,
      Grouping = ifelse(dataset_label == "Known-sex", "Documented sex", "Predicted sex"),
      Group = recode(Group, F = "Female", M = "Male"),
      mean_postprob = round(mean_postprob, 4)
    ) %>%
    dplyr::select(Dataset, Region, Grouping, Group, n, mean_postprob, all_of(traits))
}

table1_cranial <- bind_rows(
  make_group_summary(known_cranial, cranial_traits, "Sex", "Known-sex", "Cranial"),
  make_group_summary(sim_cranial, cranial_traits, "CasePrediction", "Simulated", "Cranial")
)

table1_pelvic <- bind_rows(
  make_group_summary(known_pelvic, pelvic_traits, "Sex", "Known-sex", "Pelvic"),
  make_group_summary(sim_pelvic, pelvic_traits, "CasePrediction", "Simulated", "Pelvic")
)

write.csv(table1_cranial, file.path(output_dir, "Table1A_Cranial_Group_Summaries.csv"), row.names = FALSE)
write.csv(table1_pelvic, file.path(output_dir, "Table1B_Pelvic_Group_Summaries.csv"), row.names = FALSE)

# ------------------------------
# 4. Random-label OOB Random Forest benchmarks
# ------------------------------
message("Running ", n_perm, " random-label permutations per dataset-region combination...")

perm_sim_cranial <- rf_random_label_benchmark(
  sim_cranial, cranial_traits, "simulated", "cranial",
  n_perm, n_trees_perm, seed + 1
)
perm_sim_pelvic <- rf_random_label_benchmark(
  sim_pelvic, pelvic_traits, "simulated", "pelvic",
  n_perm, n_trees_perm, seed + 2
)
perm_known_cranial <- rf_random_label_benchmark(
  known_cranial, cranial_traits, "known", "cranial",
  n_perm, n_trees_perm, seed + 3
)
perm_known_pelvic <- rf_random_label_benchmark(
  known_pelvic, pelvic_traits, "known", "pelvic",
  n_perm, n_trees_perm, seed + 4
)

permutation_results <- bind_rows(
  perm_sim_cranial,
  perm_sim_pelvic,
  perm_known_cranial,
  perm_known_pelvic
)

write.csv(
  permutation_results,
  file.path(output_dir, "Permutation_OOB_Accuracy_All.csv"),
  row.names = FALSE
)

permutation_summary <- permutation_results %>%
  group_by(dataset, region) %>%
  summarise(
    n_permutations = n(),
    theoretical_chance = 0.50,
    mean_oob_accuracy = mean(oob_accuracy),
    sd_oob_accuracy = sd(oob_accuracy),
    p95_oob_accuracy = as.numeric(quantile(oob_accuracy, 0.95)),
    p99_oob_accuracy = as.numeric(quantile(oob_accuracy, 0.99)),
    mean_oob_balanced_accuracy = mean(oob_balanced_accuracy),
    .groups = "drop"
  ) %>%
  mutate(
    benchmark_type = ifelse(
      dataset == "simulated",
      "Trait-signal ceiling under no assigned-label signal",
      "Random-label benchmark after disrupting documented sex"
    )
  )

table2_main <- permutation_summary %>%
  filter(dataset == "simulated") %>%
  mutate(
    Dataset = "Simulated",
    Region = str_to_title(region),
    n = ifelse(region == "cranial", nrow(sim_cranial), nrow(sim_pelvic)),
    n_traits = 6
  ) %>%
  dplyr::select(
    Dataset, Region, n, n_traits, theoretical_chance,
    mean_oob_accuracy, sd_oob_accuracy,
    p95_oob_accuracy, p99_oob_accuracy
  )

write.csv(
  table2_main,
  file.path(output_dir, "Table2_Trait_Signal_Ceiling_Simulated.csv"),
  row.names = FALSE
)
write.csv(
  permutation_summary,
  file.path(output_dir, "TableS2_Random_Label_Benchmarks_All.csv"),
  row.names = FALSE
)

# ------------------------------
# 5. Table 3: MorphoPASSE posterior summaries
# ------------------------------
table3 <- dat_clean %>%
  group_by(dataset, region) %>%
  summarise(
    n = n(),
    predicted_female = sum(CasePrediction == "F", na.rm = TRUE),
    predicted_male = sum(CasePrediction == "M", na.rm = TRUE),
    mean_postprob = mean(postprob, na.rm = TRUE),
    median_postprob = median(postprob, na.rm = TRUE),
    proportion_ge_0_90 = mean(postprob >= 0.90, na.rm = TRUE),
    proportion_ge_0_95 = mean(postprob >= 0.95, na.rm = TRUE),
    min_postprob = min(postprob, na.rm = TRUE),
    max_postprob = max(postprob, na.rm = TRUE),
    observed_accuracy = ifelse(
      first(dataset) == "known",
      mean(CasePrediction == Sex, na.rm = TRUE),
      NA_real_
    ),
    agreement_with_assigned_random_labels = ifelse(
      first(dataset) == "simulated",
      mean(CasePrediction == RandomSex, na.rm = TRUE),
      NA_real_
    ),
    .groups = "drop"
  ) %>%
  left_join(
    permutation_summary %>%
      dplyr::select(dataset, region, mean_oob_accuracy, p95_oob_accuracy, p99_oob_accuracy),
    by = c("dataset", "region")
  ) %>%
  mutate(
    Dataset = recode(dataset, simulated = "Simulated", known = "Known-sex"),
    Region = str_to_title(region),
    benchmark_95_99 = sprintf("%.3f-%.3f", p95_oob_accuracy, p99_oob_accuracy)
  ) %>%
  dplyr::select(
    Dataset, Region, n, predicted_female, predicted_male,
    mean_postprob, median_postprob,
    proportion_ge_0_90, proportion_ge_0_95,
    min_postprob, max_postprob,
    observed_accuracy, agreement_with_assigned_random_labels,
    mean_oob_accuracy, benchmark_95_99
  )

write.csv(
  table3,
  file.path(output_dir, "Table3_MorphoPASSE_Posterior_Summaries.csv"),
  row.names = FALSE
)

# ------------------------------
# 6. Table 4: independent Random Forest comparison
# ------------------------------
table4 <- bind_rows(
  fit_known_independent_rf(
    known_cranial, cranial_traits, "cranial", n_trees_independent, seed + 101
  ),
  fit_known_independent_rf(
    known_pelvic, pelvic_traits, "pelvic", n_trees_independent, seed + 102
  )
) %>%
  mutate(
    Dataset = "Known-sex",
    Region = str_to_title(region)
  ) %>%
  dplyr::select(
    Dataset, Region, n,
    morphopasse_accuracy,
    morphopasse_balanced_accuracy,
    independent_rf_oob_accuracy,
    independent_rf_oob_balanced_accuracy,
    prediction_agreement,
    spearman_rho_mp_postprob_vs_rf_vote_for_mp_class,
    spearman_rho_mp_postprob_vs_rf_max_vote
  )

write.csv(
  table4,
  file.path(output_dir, "Table4_Independent_RF_Known_Sex.csv"),
  row.names = FALSE
)

# Optional supplement: internal consistency within the simulated datasets.
table_s1 <- bind_rows(
  fit_simulated_internal_consistency(
    sim_cranial, cranial_traits, "cranial", n_trees_independent, seed + 201
  ),
  fit_simulated_internal_consistency(
    sim_pelvic, pelvic_traits, "pelvic", n_trees_independent, seed + 202
  )
) %>%
  mutate(
    Dataset = "Simulated",
    Region = str_to_title(region)
  ) %>%
  dplyr::select(
    Dataset, Region, n,
    oob_accuracy_predicting_morphopasse_class,
    oob_balanced_accuracy_predicting_morphopasse_class,
    spearman_rho_mp_postprob_vs_rf_vote_for_mp_class
  )

write.csv(
  table_s1,
  file.path(output_dir, "TableS1_Simulated_Internal_Consistency.csv"),
  row.names = FALSE
)

# Write all tables to one workbook when writexl is available.
if (requireNamespace("writexl", quietly = TRUE)) {
  writexl::write_xlsx(
    list(
      Table1_Cranial = table1_cranial,
      Table1_Pelvic = table1_pelvic,
      Table2_Simulated_Ceiling = table2_main,
      Table3_Posteriors = table3,
      Table4_Independent_RF = table4,
      TableS1_Sim_Internal = table_s1,
      TableS2_All_Benchmarks = permutation_summary
    ),
    path = file.path(output_dir, "Manuscript_Tables_Revised.xlsx")
  )
}

# ------------------------------
# 7. Correlation summaries for the revised simulation
# ------------------------------
cor_known_cranial <- cor(known_cranial[cranial_traits], method = "spearman", use = "pairwise.complete.obs")
cor_sim_cranial <- cor(sim_cranial[cranial_traits], method = "spearman", use = "pairwise.complete.obs")
cor_known_pelvic <- cor(known_pelvic[pelvic_traits], method = "spearman", use = "pairwise.complete.obs")
cor_sim_pelvic <- cor(sim_pelvic[pelvic_traits], method = "spearman", use = "pairwise.complete.obs")

bilateral_pairs <- tibble(
  region = c("cranial", "cranial", "pelvic", "pelvic", "pelvic"),
  left = c("SOL", "MPL", "MAL", "VAL", "SPCL"),
  right = c("SOR", "MPR", "MAR", "VAR", "SPCR")
)

table_s3_bilateral <- bilateral_pairs %>%
  rowwise() %>%
  mutate(
    known_rho = if (region == "cranial") {
      cor_known_cranial[left, right]
    } else {
      cor_known_pelvic[left, right]
    },
    simulated_rho = if (region == "cranial") {
      cor_sim_cranial[left, right]
    } else {
      cor_sim_pelvic[left, right]
    },
    absolute_difference = abs(known_rho - simulated_rho)
  ) %>%
  ungroup() %>%
  mutate(
    Region = str_to_title(region),
    Trait_pair = paste(left, right, sep = "-")
  ) %>%
  dplyr::select(Region, Trait_pair, known_rho, simulated_rho, absolute_difference)

write.csv(
  table_s3_bilateral,
  file.path(output_dir, "TableS3_Bilateral_Correlations_Known_vs_Simulated.csv"),
  row.names = FALSE
)

# ------------------------------
# 8. Figure 1: conceptual trait structure figure
# ------------------------------
set.seed(seed)
n_concept <- 120
latent_trait <- rnorm(n_concept)

concept_dat <- tibble(
  individual = seq_len(n_concept),
  trait_1 = latent_trait + rnorm(n_concept, 0, 0.55),
  trait_2 = latent_trait + rnorm(n_concept, 0, 0.55),
  random_sex = sample(rep(c("Female", "Male"), each = n_concept / 2))
)

p1a <- ggplot(concept_dat, aes(trait_1, trait_2)) +
  geom_point(shape = 21, size = 2.8, stroke = 0.7, fill = "grey75") +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8, linetype = "dashed") +
  labs(
    title = "A. Correlated trait structure",
    x = "Trait 1 expression",
    y = "Trait 2 expression"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

p1b <- ggplot(concept_dat, aes(trait_1, trait_2, shape = random_sex)) +
  geom_point(size = 2.8, stroke = 0.9) +
  labs(
    title = "B. Balanced random-label assignment",
    x = "Trait 1 expression",
    y = "Trait 2 expression",
    shape = "Assigned label"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

flow_nodes <- tibble(
  x = c(1.0, 2.5, 4.0, 5.5),
  y = 1,
  label = c(
    "Generate correlated\nordinal traits",
    "Assign balanced labels\nafter trait generation",
    "Estimate OOB accuracy\nunder repeated random labels",
    "Compare empirical benchmark\nwith MorphoPASSE confidence"
  )
)

p1c <- ggplot(flow_nodes, aes(x, y)) +
  geom_label(
    aes(label = label),
    size = 3.7,
    label.size = 0.6,
    fill = "white",
    label.padding = grid::unit(0.25, "lines")
  ) +
  geom_segment(
    data = tibble(
      x = c(1.55, 3.05, 4.55),
      xend = c(1.95, 3.45, 4.95),
      y = 1,
      yend = 1
    ),
    aes(x = x, xend = xend, y = y, yend = yend),
    inherit.aes = FALSE,
    arrow = arrow(length = grid::unit(0.15, "inches")),
    linewidth = 0.8
  ) +
  annotate(
    "text", x = 3.25, y = 0.62,
    label = "Trait structure is retained; assigned sex-label signal remains absent",
    fontface = "italic", size = 4
  ) +
  coord_cartesian(xlim = c(0.35, 6.15), ylim = c(0.45, 1.25), clip = "off") +
  labs(title = "C. Analysis logic") +
  theme_void(base_size = 13) +
  theme(plot.title = element_text(face = "bold", hjust = 0))

figure1 <- (p1a | p1b) / p1c +
  plot_layout(heights = c(1, 0.72)) +
  plot_annotation(
    title = "Trait structures and sex signals",
    theme = theme(
      plot.title = element_text(face = "bold", size = 18),
      plot.subtitle = element_text(size = 13)
    )
  )

print(figure1)

ggsave(
  file.path(output_dir, "Figure1_Trait_Structure_vs_Sex_Signal.png"),
  figure1,
  width = 11,
  height = 8.5,
  dpi = 300,
  bg = "white"
)

# ------------------------------
# 9. Figure 2: empirical multi-panel results
# ------------------------------
perm_plot_df <- permutation_results %>%
  filter(dataset == "simulated") %>%
  mutate(region = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic")))

perm_plot_summary <- permutation_summary %>%
  filter(dataset == "simulated") %>%
  mutate(region = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic")))

p2a <- ggplot(perm_plot_df, aes(oob_accuracy)) +
  geom_density(linewidth = 0.9) +
  geom_vline(data = perm_plot_summary, aes(xintercept = mean_oob_accuracy), linewidth = 0.9) +
  geom_vline(data = perm_plot_summary, aes(xintercept = p95_oob_accuracy), linetype = "dashed", linewidth = 0.8) +
  geom_vline(data = perm_plot_summary, aes(xintercept = p99_oob_accuracy), linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = 0.50, linetype = "longdash", linewidth = 0.8, color = "firebrick") +
  facet_wrap(~region, ncol = 2) +
  labs(
    title = "A. Random-label trait-signal benchmark",
      x = "Out-of-bag classification accuracy",
    y = "Density"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

post_plot_df <- dat_clean %>%
  mutate(
    dataset_label = factor(dataset, levels = c("simulated", "known"), labels = c("Simulated", "Known-sex")),
    region_label = factor(region, levels = c("cranial", "pelvic"), labels = c("Cranial", "Pelvic"))
  )

p2b <- ggplot(post_plot_df, aes(postprob)) +
  geom_histogram(bins = 15, boundary = 0.5, closed = "left") +
  geom_vline(xintercept = 0.95, linetype = "dashed", linewidth = 0.6) +
  facet_grid(dataset_label ~ region_label) +
  coord_cartesian(xlim = c(0.50, 1.00)) +
  labs(
    title = "B. MorphoPASSE posterior probabilities",
    x = "Posterior probability for the predicted class",
    y = "Number of individuals"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank()
  )

summary_plot_df <- table3 %>%
  transmute(
    dataset = factor(Dataset, levels = c("Simulated", "Known-sex")),
    region = factor(Region, levels = c("Cranial", "Pelvic")),
    post_mean = mean_postprob,
    post_median = median_postprob,
    accuracy = observed_accuracy,
    benchmark_mean = mean_oob_accuracy
  ) %>%
  left_join(
    permutation_summary %>%
      transmute(
        dataset = factor(recode(dataset, simulated = "Simulated", known = "Known-sex"), levels = c("Simulated", "Known-sex")),
        region = factor(str_to_title(region), levels = c("Cranial", "Pelvic")),
        benchmark_p95 = p95_oob_accuracy,
        benchmark_p99 = p99_oob_accuracy
      ),
    by = c("dataset", "region")
  ) %>%
  mutate(x_pos = as.numeric(region))

summary_points <- summary_plot_df %>%
  dplyr::select(dataset, region, x_pos, post_mean, post_median, benchmark_mean, accuracy) %>%
  pivot_longer(
    cols = c(benchmark_mean, post_mean, post_median, accuracy),
    names_to = "stat",
    values_to = "value"
  ) %>%
  filter(!is.na(value)) %>%
  mutate(
    stat = recode(
      stat,
      benchmark_mean = "Random-label benchmark mean",
      post_mean = "MorphoPASSE mean postprob",
      post_median = "MorphoPASSE median postprob",
      accuracy = "Observed accuracy (known-sex only)"
    ),
    label = sprintf("%.2f", value),
    point_size = ifelse(stat == "Random-label benchmark mean", 2.7, 3.5)
  )

p2c <- ggplot(summary_plot_df, aes(x = x_pos)) +
  geom_hline(yintercept = 0.50, linewidth = 0.9, linetype = "dashed", color = "firebrick") +
  annotate(
    "curve",
    x = 1.35, y = 0.60, xend = 1.03, yend = 0.505,
    curvature = -0.35,
    arrow = arrow(length = grid::unit(0.16, "inches")),
    linewidth = 0.7,
    color = "firebrick"
  ) +
  annotate(
    "text",
    x = 1.38, y = 0.62,
    label = "Random expectation = 0.50",
    hjust = 0,
    size = 3.3,
    color = "firebrick"
  ) +
  geom_linerange(
    aes(ymin = benchmark_p95, ymax = benchmark_p99),
    linewidth = 10,
    alpha = 0.18
  ) +
  geom_point(aes(y = benchmark_p95), shape = 95, size = 10, alpha = 0.55) +
  geom_point(
    data = summary_points,
    aes(x = x_pos, y = value, shape = stat, size = point_size),
    inherit.aes = FALSE,
    stroke = 1.1
  ) +
  scale_size_identity() +
  geom_text(
    data = summary_points %>% filter(stat != "Random-label benchmark mean"),
    aes(x = x_pos, y = value, label = label),
    inherit.aes = FALSE,
    nudge_x = 0.18,
    vjust = -0.5,
    size = 3.2,
    fontface = "bold"
  ) +
  facet_wrap(~dataset) +
  coord_cartesian(ylim = c(0.45, 1.02), clip = "off") +
  scale_x_continuous(breaks = c(1, 2), labels = c("Cranial", "Pelvic")) +
  scale_shape_manual(values = c(
    "Random-label benchmark mean" = 16,
    "MorphoPASSE mean postprob" = 17,
    "MorphoPASSE median postprob" = 15,
    "Observed accuracy (known-sex only)" = 1
  )) +
  labs(
    title = "C. Random-label benchmark versus posterior confidence",
    x = NULL,
    y = "Accuracy or probability",
    shape = NULL
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    panel.grid.minor = element_blank(),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 25, 10, 10)
  )

figure2 <- (p2a / p2b / p2c) + plot_layout(heights = c(1.0, 1.15, 1.1))

print(figure2)

ggsave(
  file.path(output_dir, "Figure2_Empirical_MultiPanel.png"),
  figure2,
  width = 11,
  height = 13,
  dpi = 300,
  bg = "white"
)

# ------------------------------
# 10. Optional Figure S1: source versus simulated correlations
# ------------------------------
make_corr_plot <- function(known_matrix, simulated_matrix, region_label, trait_order) {
  corr_long <- bind_rows(
    matrix_to_long(known_matrix, "Known-sex", region_label, trait_order),
    matrix_to_long(simulated_matrix, "Simulated", region_label, trait_order)
  ) %>%
    mutate(dataset = factor(dataset, levels = c("Known-sex", "Simulated")))

  ggplot(corr_long, aes(trait_1, trait_2, fill = rho)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 3.0) +
    scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Spearman rho") +
    facet_wrap(~dataset, ncol = 2) +
    labs(
      title = paste0(region_label, " trait-correlation structure"),
      x = NULL,
      y = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1),
      strip.text = element_text(face = "bold"),
      legend.position = "bottom"
    )
}

p_corr_cranial <- make_corr_plot(
  cor_known_cranial, cor_sim_cranial, "Cranial", cranial_traits
)
p_corr_pelvic <- make_corr_plot(
  cor_known_pelvic, cor_sim_pelvic, "Pelvic", pelvic_traits
)

figure_s1 <- p_corr_cranial / p_corr_pelvic +
  plot_annotation(
    title = "Known-sex source and simulated trait correlations",
    theme = theme(plot.title = element_text(face = "bold", size = 17))
  )

print(figure_s1)

ggsave(
  file.path(output_dir, "FigureS1_Correlation_Structure.png"),
  figure_s1,
  width = 11,
  height = 10,
  dpi = 300,
  bg = "white"
)

# ------------------------------
# 11. Captions, summary, and reproducibility files
# ------------------------------
caption_text <- c(
  "FIGURE 1. Trait structure is not equivalent to sex signal. (A) Morphological traits may covary even when no sex differences are introduced. (B) Balanced random assignment of sex labels after trait generation does not create a true relationship between the labels and the simulated morphology. (C) The analysis retains the observed trait structure while using out-of-bag Random Forest performance under repeated random labels as an empirical benchmark for comparison with MorphoPASSE posterior confidence.",
  "",
  "FIGURE 2. Random-label classification benchmarks and MorphoPASSE posterior confidence. (A) Out-of-bag Random Forest accuracy across 500 balanced random-label permutations of the revised simulated cranial and pelvic datasets. Solid, dashed, and dotted vertical lines denote the mean, 95th percentile, and 99th percentile, respectively; the firebrick line denotes the theoretical 0.50 expectation. (B) MorphoPASSE posterior probabilities for simulated and known-sex data. (C) Mean and median MorphoPASSE posterior probabilities compared with random-label benchmark distributions and observed classification accuracy for the known-sex data.",
  "",
  "FIGURE S1. Spearman trait-correlation matrices for the deduplicated known-sex source data and revised simulated datasets. The revised simulations were generated without sex labels but preserve bilateral and within-module correlation patterns estimated from the pooled known-sex data.",
  "",
  "TABLE 1. Group-specific mean MorphoPASSE posterior probabilities and median morphoscopic trait scores with observed ranges for known-sex and revised simulated cranial and pelvic samples.",
  "",
  "TABLE 2. Out-of-bag Random Forest classification performance for revised simulated cranial and pelvic trait sets under 500 balanced random sex-label assignments. The mean and upper percentiles provide the empirical no-signal benchmark used to evaluate MorphoPASSE posterior confidence.",
  "",
  "TABLE 3. Distributional properties of MorphoPASSE posterior probabilities by dataset and skeletal region, including observed accuracy for known-sex data and agreement with the independently assigned random labels for simulated data.",
  "",
  "TABLE 4. MorphoPASSE performance and an independent, class-balanced Random Forest analysis of the known-sex datasets. Spearman correlations compare MorphoPASSE posterior probabilities with independent out-of-bag vote proportions supporting the MorphoPASSE-predicted class."
)

writeLines(caption_text, file.path(output_dir, "Table_and_Figure_Captions.txt"))

summary_lines <- c(
  "MORPHOPASSE REVISED ANALYSIS SUMMARY",
  paste0("Input workbook: ", normalizePath(input_workbook)),
  paste0("Random-label file: ", normalizePath(random_label_file)),
  paste0("Known cranial unique N: ", nrow(known_cranial)),
  paste0("Known pelvic unique N: ", nrow(known_pelvic)),
  paste0("Simulated cranial N: ", nrow(sim_cranial)),
  paste0("Simulated pelvic N: ", nrow(sim_pelvic)),
  paste0("Permutations per dataset-region: ", n_perm),
  paste0("Trees per permutation forest: ", n_trees_perm),
  paste0("Trees per independent forest: ", n_trees_independent),
  "",
  "Permutation summaries:",
  capture.output(print(permutation_summary)),
  "",
  "Posterior summaries:",
  capture.output(print(table3)),
  "",
  "Independent RF comparison:",
  capture.output(print(table4))
)

writeLines(summary_lines, file.path(output_dir, "Analysis_Summary.txt"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_Session_Info.txt"))

message("Analysis complete. Outputs saved to: ", normalizePath(output_dir))
