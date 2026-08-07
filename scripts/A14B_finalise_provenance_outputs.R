# ============================================================
# A14B FINALISE PROVENANCE REVIEW OUTPUTS
#
# Purpose:
# - Rebuild the duplicate-value review from the saved enriched file
# - Create a provenance summary
# - Write and verify the two missing A14 outputs
#
# Input:
#   results/A14_activity_provenance_enriched.csv
#   results/A14_assay_annotations.csv
#   results/A14_document_annotations.csv
#
# Outputs:
#   results/A14_duplicate_value_review.csv
#   results/A14_provenance_summary.csv
# ============================================================

required_packages <- c("readr", "dplyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them with:\ninstall.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

activity_file <- "results/A14_activity_provenance_enriched.csv"
assay_file <- "results/A14_assay_annotations.csv"
document_file <- "results/A14_document_annotations.csv"
duplicate_output_file <- "results/A14_duplicate_value_review.csv"
summary_output_file <- "results/A14_provenance_summary.csv"

input_files <- c(activity_file, assay_file, document_file)
missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n")
    )
  )
}

activities <- readr::read_csv(
  activity_file,
  show_col_types = FALSE,
  progress = FALSE
)

assays <- readr::read_csv(
  assay_file,
  show_col_types = FALSE,
  progress = FALSE
)

documents <- readr::read_csv(
  document_file,
  show_col_types = FALSE,
  progress = FALSE
)

# Add optional columns if an API response did not provide them.
optional_activity_columns <- c(
  "activity_id",
  "molecule_chembl_id",
  "queried_protein",
  "standard_type",
  "standard_value_numeric",
  "standard_units",
  "assay_chembl_id",
  "final_document_chembl_id",
  "final_assay_type",
  "assay_description",
  "confidence_score",
  "relationship_type",
  "doi",
  "pubmed_id",
  "direct_binding_assay",
  "functional_assay",
  "high_confidence_target_assignment",
  "single_protein_high_confidence"
)

for (column_name in optional_activity_columns) {
  if (!column_name %in% names(activities)) {
    activities[[column_name]] <- NA
  }
}

optional_assay_columns <- c(
  "assay_chembl_id",
  "assay_retrieved",
  "confidence_score"
)

for (column_name in optional_assay_columns) {
  if (!column_name %in% names(assays)) {
    assays[[column_name]] <- NA
  }
}

optional_document_columns <- c(
  "document_chembl_id",
  "document_retrieved"
)

for (column_name in optional_document_columns) {
  if (!column_name %in% names(documents)) {
    documents[[column_name]] <- NA
  }
}

collapse_unique <- function(value) {
  value <- value[!is.na(value)]
  value <- sort(unique(as.character(value)))

  if (length(value) == 0) {
    return(NA_character_)
  }

  paste(value, collapse = "; ")
}

count_nonmissing_distinct <- function(value) {
  value <- value[!is.na(value) & as.character(value) != ""]
  dplyr::n_distinct(value)
}

# Quantitative records only.
quantitative_records <- activities |>
  dplyr::mutate(
    standard_value_numeric = suppressWarnings(
      as.numeric(standard_value_numeric)
    ),
    standard_type = as.character(standard_type),
    standard_units = as.character(standard_units),
    assay_chembl_id = as.character(assay_chembl_id),
    final_document_chembl_id = as.character(final_document_chembl_id)
  ) |>
  dplyr::filter(
    is.finite(standard_value_numeric),
    standard_type %in% c("IC50", "Kd", "ED50")
  )

# Build the duplicate-value review. The .drop = FALSE behaviour is not
# required because an explicit empty schema is created when no groups match.
if (nrow(quantitative_records) == 0) {
  duplicate_value_review <- tibble::tibble(
    molecule_chembl_id = character(),
    queried_protein = character(),
    standard_value_numeric = double(),
    standard_units = character(),
    record_count = integer(),
    activity_types = character(),
    activity_type_count = integer(),
    activity_ids = character(),
    assay_ids = character(),
    unique_assay_count = integer(),
    document_ids = character(),
    unique_document_count = integer(),
    assay_types = character(),
    assay_descriptions = character(),
    confidence_scores = character(),
    relationship_types = character(),
    dois = character(),
    pubmed_ids = character(),
    has_same_value_across_activity_types = logical(),
    likely_same_assay_reannotation = logical(),
    likely_same_document_evidence = logical()
  )
} else {
  duplicate_value_review <- quantitative_records |>
    dplyr::group_by(
      molecule_chembl_id,
      queried_protein,
      standard_value_numeric,
      standard_units
    ) |>
    dplyr::summarise(
      record_count = dplyr::n(),
      activity_types = collapse_unique(standard_type),
      activity_type_count = dplyr::n_distinct(standard_type),
      activity_ids = collapse_unique(activity_id),
      assay_ids = collapse_unique(assay_chembl_id),
      unique_assay_count = count_nonmissing_distinct(assay_chembl_id),
      document_ids = collapse_unique(final_document_chembl_id),
      unique_document_count = count_nonmissing_distinct(final_document_chembl_id),
      assay_types = collapse_unique(final_assay_type),
      assay_descriptions = collapse_unique(assay_description),
      confidence_scores = collapse_unique(confidence_score),
      relationship_types = collapse_unique(relationship_type),
      dois = collapse_unique(doi),
      pubmed_ids = collapse_unique(pubmed_id),
      has_same_value_across_activity_types =
        dplyr::n_distinct(standard_type) > 1,
      likely_same_assay_reannotation =
        dplyr::n_distinct(standard_type) > 1 &&
        count_nonmissing_distinct(assay_chembl_id) <= 1,
      likely_same_document_evidence =
        dplyr::n_distinct(standard_type) > 1 &&
        count_nonmissing_distinct(final_document_chembl_id) <= 1,
      .groups = "drop"
    ) |>
    dplyr::filter(
      record_count > 1 |
        has_same_value_across_activity_types
    ) |>
    dplyr::arrange(
      dplyr::desc(has_same_value_across_activity_types),
      molecule_chembl_id,
      queried_protein,
      standard_value_numeric
    )
}

# Robust logical counting even when columns were imported as text.
as_logical_flag <- function(value) {
  if (is.logical(value)) {
    return(value)
  }

  as.character(value) %in% c("TRUE", "True", "true", "1")
}

assays_retrieved <- sum(
  as_logical_flag(assays$assay_retrieved),
  na.rm = TRUE
)

documents_retrieved <- sum(
  as_logical_flag(documents$document_retrieved),
  na.rm = TRUE
)

binding_records <- sum(
  as_logical_flag(activities$direct_binding_assay),
  na.rm = TRUE
)

functional_records <- sum(
  as_logical_flag(activities$functional_assay),
  na.rm = TRUE
)

high_confidence_records <- sum(
  as_logical_flag(activities$high_confidence_target_assignment),
  na.rm = TRUE
)

single_protein_records <- sum(
  as_logical_flag(activities$single_protein_high_confidence),
  na.rm = TRUE
)

provenance_summary <- tibble::tibble(
  metric = c(
    "Candidate activity records curated",
    "Unique assays",
    "Assays successfully retrieved",
    "Unique source documents",
    "Documents successfully retrieved",
    "Binding-assay records",
    "Functional-assay records",
    "High-confidence target assignments (score >=8)",
    "Single-protein target assignments (score 9)",
    "Duplicate-value groups requiring review",
    "Identical values across multiple activity types",
    "Likely same-assay reannotations",
    "Likely same-document evidence groups"
  ),
  count = c(
    nrow(activities),
    count_nonmissing_distinct(activities$assay_chembl_id),
    assays_retrieved,
    count_nonmissing_distinct(activities$final_document_chembl_id),
    documents_retrieved,
    binding_records,
    functional_records,
    high_confidence_records,
    single_protein_records,
    nrow(duplicate_value_review),
    sum(
      duplicate_value_review$has_same_value_across_activity_types,
      na.rm = TRUE
    ),
    sum(
      duplicate_value_review$likely_same_assay_reannotation,
      na.rm = TRUE
    ),
    sum(
      duplicate_value_review$likely_same_document_evidence,
      na.rm = TRUE
    )
  )
)

# Use base write.csv for maximum compatibility with zero-row outputs.
utils::write.csv(
  as.data.frame(duplicate_value_review),
  duplicate_output_file,
  row.names = FALSE,
  na = ""
)

utils::write.csv(
  as.data.frame(provenance_summary),
  summary_output_file,
  row.names = FALSE,
  na = ""
)

output_files <- c(
  duplicate_output_file,
  summary_output_file
)

file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(file.info(output_files)$size)
)

if (!all(file_check$exists)) {
  stop("One or more A14B output files were not created.")
}

if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) {
  stop("One or more A14B output files are empty or invalid.")
}

cat("\n")
cat("A14 provenance outputs finalised successfully.\n")
cat("---------------------------------------------\n")

print(provenance_summary, n = Inf)

cat("\nDuplicate-value groups requiring review:\n")
print(
  duplicate_value_review |>
    dplyr::select(
      molecule_chembl_id,
      queried_protein,
      standard_value_numeric,
      standard_units,
      activity_types,
      assay_ids,
      document_ids,
      has_same_value_across_activity_types,
      likely_same_assay_reannotation,
      likely_same_document_evidence
    ),
  n = Inf,
  width = Inf
)

cat("\nOutput verification:\n")
print(file_check, n = Inf)
