# ============================================================
# A10 FILTER AND CLASSIFY ChEMBL ACTIVITY RECORDS
#
# Input:
#   results/A9_chembl_activities_complete.csv
#
# Outputs:
#   results/A10_activity_quality_controlled.csv
#   results/A10_quantitative_activities.csv
#   results/A10_priority_inhibitory_activities.csv
#   results/A10_compound_protein_summary.csv
#   results/A10_activity_filter_summary.csv
#   results/A10_candidate_molecule_ids.csv
# ============================================================

# ============================================================
# 1. CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them using:\n",
      "install.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

# ============================================================
# 2. DEFINE FILE PATHS
# ============================================================

input_file <-
  "results/A9_chembl_activities_complete.csv"

quality_controlled_file <-
  "results/A10_activity_quality_controlled.csv"

quantitative_file <-
  "results/A10_quantitative_activities.csv"

priority_file <-
  "results/A10_priority_inhibitory_activities.csv"

compound_protein_summary_file <-
  "results/A10_compound_protein_summary.csv"

filter_summary_file <-
  "results/A10_activity_filter_summary.csv"

candidate_ids_file <-
  "results/A10_candidate_molecule_ids.csv"

output_files <- c(
  quality_controlled_file,
  quantitative_file,
  priority_file,
  compound_protein_summary_file,
  filter_summary_file,
  candidate_ids_file
)

if (!dir.exists("results")) {
  dir.create(
    "results",
    recursive = TRUE
  )
}

# ============================================================
# 3. LOAD INPUT DATA
# ============================================================

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ",
      input_file,
      "\n\nCurrent working directory:\n",
      getwd(),
      "\n\nRun this script from the RESKO project root."
    )
  )
}

activities <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(activities) == 0) {
  stop(
    "The ChEMBL activity file contains no records."
  )
}

# ============================================================
# 4. CHECK ESSENTIAL COLUMNS
# ============================================================

essential_columns <- c(
  "queried_protein",
  "queried_target_chembl_id",
  "molecule_chembl_id",
  "standard_type",
  "standard_value",
  "standard_units"
)

missing_columns <- setdiff(
  essential_columns,
  names(activities)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing essential columns: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# ============================================================
# 5. ADD OPTIONAL COLUMNS WHEN ABSENT
# ============================================================

optional_columns <- c(
  "standard_relation",
  "pchembl_value",
  "activity_comment",
  "data_validity_comment",
  "potential_duplicate",
  "assay_chembl_id",
  "assay_type",
  "document_chembl_id",
  "bao_label",
  "target_organism",
  "target_pref_name",
  "canonical_smiles"
)

for (column_name in optional_columns) {
  if (!column_name %in% names(activities)) {
    activities[[column_name]] <- NA_character_
  }
}

# ============================================================
# 6. HELPER FUNCTIONS
# ============================================================

normalise_text <- function(value) {
  value <- as.character(value)
  value <- stringr::str_trim(value)
  value[value %in% c("", "NA", "NaN", "NULL")] <- NA_character_
  value
}

collapse_unique <- function(value) {
  value <- value[!is.na(value)]
  value <- unique(as.character(value))
  value <- sort(value)

  if (length(value) == 0) {
    return(NA_character_)
  }

  paste(value, collapse = "; ")
}

safe_numeric_min <- function(value) {
  value <- value[is.finite(value)]

  if (length(value) == 0) {
    return(NA_real_)
  }

  min(value)
}

safe_numeric_median <- function(value) {
  value <- value[is.finite(value)]

  if (length(value) == 0) {
    return(NA_real_)
  }

  stats::median(value)
}

safe_numeric_max <- function(value) {
  value <- value[is.finite(value)]

  if (length(value) == 0) {
    return(NA_real_)
  }

  max(value)
}

make_csv_safe <- function(data) {
  data <- dplyr::ungroup(data)
  data <- as.data.frame(data, stringsAsFactors = FALSE)

  for (column_name in names(data)) {
    column <- data[[column_name]]

    if (is.list(column)) {
      data[[column_name]] <- vapply(
        column,
        function(value) {
          value <- unlist(value, recursive = TRUE, use.names = FALSE)
          value <- value[!is.na(value)]

          if (length(value) == 0) {
            return(NA_character_)
          }

          paste(as.character(value), collapse = "; ")
        },
        character(1)
      )
    } else if (is.factor(column)) {
      data[[column_name]] <- as.character(column)
    }
  }

  data
}

write_csv_verified <- function(data, path) {
  safe_data <- make_csv_safe(data)

  readr::write_csv(
    safe_data,
    path,
    na = ""
  )

  if (!file.exists(path)) {
    stop(
      paste0(
        "Output file was not created: ",
        path
      )
    )
  }

  file_size <- file.info(path)$size

  if (is.na(file_size) || file_size <= 0) {
    stop(
      paste0(
        "Output file is empty or invalid: ",
        path
      )
    )
  }

  invisible(path)
}

# ============================================================
# 7. CLEAN AND STANDARDISE CORE FIELDS
# ============================================================

activities_qc <- activities |>
  dplyr::mutate(
    queried_protein = normalise_text(queried_protein),
    queried_target_chembl_id = normalise_text(queried_target_chembl_id),
    molecule_chembl_id = normalise_text(molecule_chembl_id),
    standard_type = normalise_text(standard_type),
    standard_units = normalise_text(standard_units),
    standard_relation = normalise_text(standard_relation),
    data_validity_comment = normalise_text(data_validity_comment),
    assay_chembl_id = normalise_text(assay_chembl_id),
    standard_value_numeric = suppressWarnings(as.numeric(standard_value)),
    pchembl_value_numeric = suppressWarnings(as.numeric(pchembl_value)),
    has_numeric_value = is.finite(standard_value_numeric),
    has_molecule_id = !is.na(molecule_chembl_id),
    is_nanomolar = standard_units == "nM",
    is_exact_or_upper_bound =
      is.na(standard_relation) |
      standard_relation %in% c("=", "<", "<="),
    has_validity_warning = !is.na(data_validity_comment),
    is_potential_duplicate =
      normalise_text(potential_duplicate) %in%
      c("1", "TRUE", "True", "true")
  )

# ============================================================
# 8. CLASSIFY ACTIVITY TYPES
# ============================================================

activities_qc <- activities_qc |>
  dplyr::mutate(
    evidence_class = dplyr::case_when(
      standard_type == "IC50" ~ "Inhibitory potency",
      standard_type == "Kd" ~ "Binding affinity",
      standard_type == "ED50" ~ "Functional dose response",
      standard_type == "Inhibition" ~ "Percentage inhibition",
      standard_type == "Ka" ~ "Association constant",
      standard_type %in% c("Activity", "FC", "Ratio") ~ "Contextual activity",
      TRUE ~ "Other activity"
    ),
    is_quantitative_concentration =
      standard_type %in% c("IC50", "Kd", "ED50") &
      is_nanomolar &
      has_numeric_value &
      has_molecule_id,
    is_inhibitory_measurement =
      standard_type == "IC50" &
      is_nanomolar &
      has_numeric_value &
      has_molecule_id,
    is_binding_measurement =
      standard_type == "Kd" &
      is_nanomolar &
      has_numeric_value &
      has_molecule_id,
    is_functional_measurement =
      standard_type == "ED50" &
      is_nanomolar &
      has_numeric_value &
      has_molecule_id
  )

# ============================================================
# 9. ASSIGN EVIDENCE TIERS
# ============================================================

activities_qc <- activities_qc |>
  dplyr::mutate(
    evidence_tier = dplyr::case_when(
      is_inhibitory_measurement &
        is_exact_or_upper_bound &
        !has_validity_warning ~
        "Tier 1: quantitative inhibition",

      is_binding_measurement &
        is_exact_or_upper_bound &
        !has_validity_warning ~
        "Tier 2: quantitative binding",

      is_functional_measurement &
        is_exact_or_upper_bound &
        !has_validity_warning ~
        "Tier 3: quantitative functional response",

      is_quantitative_concentration ~
        "Quantitative record requiring review",

      TRUE ~
        "Contextual record"
    ),
    potency_band = dplyr::case_when(
      !is_quantitative_concentration ~
        "Not concentration-ranked",
      standard_value_numeric <= 10 ~
        "Very high apparent potency (<=10 nM)",
      standard_value_numeric <= 100 ~
        "High apparent potency (>10 to 100 nM)",
      standard_value_numeric <= 1000 ~
        "Moderate apparent potency (>100 to 1000 nM)",
      standard_value_numeric <= 10000 ~
        "Lower apparent potency (>1 to 10 uM)",
      standard_value_numeric > 10000 ~
        "Weak apparent potency (>10 uM)",
      TRUE ~
        "Unclassified"
    )
  )

# ============================================================
# 10. CREATE QUANTITATIVE ACTIVITY SUBSET
# ============================================================

quantitative_activities <- activities_qc |>
  dplyr::filter(
    is_quantitative_concentration
  ) |>
  dplyr::arrange(
    queried_protein,
    standard_type,
    standard_value_numeric
  )

# ============================================================
# 11. CREATE PRIORITY INHIBITORY SUBSET
#
# The 10,000 nM threshold is deliberately permissive for the
# first-pass network and can be tightened later.
# ============================================================

priority_inhibitory_activities <- activities_qc |>
  dplyr::filter(
    is_inhibitory_measurement,
    is_exact_or_upper_bound,
    !has_validity_warning,
    standard_value_numeric <= 10000
  ) |>
  dplyr::arrange(
    queried_protein,
    standard_value_numeric
  )

# ============================================================
# 12. AGGREGATE REPEATED MOLECULE-PROTEIN MEASUREMENTS
#
# IC50, Kd and ED50 remain separate activity types.
# ============================================================

if (nrow(quantitative_activities) == 0) {
  compound_protein_summary <- tibble::tibble(
    molecule_chembl_id = character(),
    queried_protein = character(),
    queried_target_chembl_id = character(),
    standard_type = character(),
    standard_units = character(),
    evidence_class = character(),
    evidence_tier = character(),
    measurement_count = integer(),
    minimum_standard_value = double(),
    median_standard_value = double(),
    maximum_standard_value = double(),
    maximum_pchembl_value = double(),
    unique_assay_count = integer(),
    validity_warning_count = integer(),
    potential_duplicate_count = integer()
  )
} else {
  compound_protein_summary <- quantitative_activities |>
    dplyr::group_by(
      molecule_chembl_id,
      queried_protein,
      queried_target_chembl_id,
      standard_type,
      standard_units
    ) |>
    dplyr::summarise(
      evidence_class = dplyr::first(evidence_class),
      evidence_tier = collapse_unique(evidence_tier),
      measurement_count = dplyr::n(),
      minimum_standard_value = safe_numeric_min(standard_value_numeric),
      median_standard_value = safe_numeric_median(standard_value_numeric),
      maximum_standard_value = safe_numeric_max(standard_value_numeric),
      maximum_pchembl_value = safe_numeric_max(pchembl_value_numeric),
      unique_assay_count = dplyr::n_distinct(
        assay_chembl_id[
          !is.na(assay_chembl_id)
        ]
      ),
      validity_warning_count = sum(
        has_validity_warning,
        na.rm = TRUE
      ),
      potential_duplicate_count = sum(
        is_potential_duplicate,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      standard_type,
      minimum_standard_value
    )
}

# ============================================================
# 13. CREATE CANDIDATE MOLECULE LIST
# ============================================================

if (nrow(compound_protein_summary) == 0) {
  candidate_molecule_ids <- tibble::tibble(
    molecule_chembl_id = character(),
    proteins = character(),
    activity_types = character(),
    strongest_record_nM = double(),
    total_measurements = integer()
  )
} else {
  candidate_molecule_ids <- compound_protein_summary |>
    dplyr::group_by(
      molecule_chembl_id
    ) |>
    dplyr::summarise(
      proteins = collapse_unique(queried_protein),
      activity_types = collapse_unique(standard_type),
      strongest_record_nM = safe_numeric_min(minimum_standard_value),
      total_measurements = sum(
        measurement_count,
        na.rm = TRUE
      ),
      .groups = "drop"
    ) |>
    dplyr::arrange(
      strongest_record_nM
    )
}

# ============================================================
# 14. CREATE FILTER SUMMARY
# ============================================================

filter_summary <- tibble::tibble(
  category = c(
    "All ChEMBL records",
    "Unique ChEMBL molecules",
    "Quantitative concentration records",
    "IC50 records in nM",
    "Kd records in nM",
    "ED50 records in nM",
    "Priority inhibitory records",
    "Compound-protein summaries",
    "Candidate molecules with quantitative evidence"
  ),
  count = c(
    nrow(activities_qc),
    dplyr::n_distinct(
      activities_qc$molecule_chembl_id[
        activities_qc$has_molecule_id
      ]
    ),
    nrow(quantitative_activities),
    sum(activities_qc$is_inhibitory_measurement, na.rm = TRUE),
    sum(activities_qc$is_binding_measurement, na.rm = TRUE),
    sum(activities_qc$is_functional_measurement, na.rm = TRUE),
    nrow(priority_inhibitory_activities),
    nrow(compound_protein_summary),
    nrow(candidate_molecule_ids)
  )
)

# ============================================================
# 15. WRITE AND VERIFY ALL OUTPUT FILES
# ============================================================

write_csv_verified(
  activities_qc,
  quality_controlled_file
)

write_csv_verified(
  quantitative_activities,
  quantitative_file
)

write_csv_verified(
  priority_inhibitory_activities,
  priority_file
)

write_csv_verified(
  compound_protein_summary,
  compound_protein_summary_file
)

write_csv_verified(
  filter_summary,
  filter_summary_file
)

write_csv_verified(
  candidate_molecule_ids,
  candidate_ids_file
)

file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(
    file.info(output_files)$size
  )
)

if (!all(file_check$exists)) {
  stop(
    "At least one A10 output file was not created."
  )
}

if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) {
  stop(
    "At least one A10 output file is empty or invalid."
  )
}

# ============================================================
# 16. PRINT RESULTS
# ============================================================

cat("\n")
cat("A10 activity filtering completed successfully.\n")
cat("----------------------------------------------\n")

print(
  filter_summary,
  n = Inf
)

cat("\nRecords by evidence tier:\n")

print(
  activities_qc |>
    dplyr::count(
      evidence_tier,
      sort = TRUE
    ),
  n = Inf
)

cat("\nPriority inhibitory records by protein:\n")

print(
  priority_inhibitory_activities |>
    dplyr::count(
      queried_protein,
      sort = TRUE
    ),
  n = Inf
)

cat("\nOutput verification:\n")

print(
  file_check,
  n = Inf
)

cat("\nA10 files created:\n")

for (file_name in output_files) {
  cat(
    "- ",
    normalizePath(
      file_name,
      winslash = "/",
      mustWork = TRUE
    ),
    "\n",
    sep = ""
  )
}
