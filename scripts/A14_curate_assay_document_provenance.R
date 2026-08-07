# ============================================================
# A14 CURATE ChEMBL ASSAY AND DOCUMENT PROVENANCE
#
# Purpose:
# - Retrieve assay metadata for quantitative candidate evidence
# - Retrieve source-document metadata
# - Add assay confidence, descriptions, formats and provenance
# - Identify records that may reflect duplicated evidence
# - Test whether identical Kd and ED50 values share provenance
#
# Input:
#   results/A10_activity_quality_controlled.csv
#
# Outputs:
#   results/A14_assay_annotations.csv
#   results/A14_document_annotations.csv
#   results/A14_activity_provenance_enriched.csv
#   results/A14_duplicate_value_review.csv
#   results/A14_provenance_summary.csv
# ============================================================

# ============================================================
# 1. CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "stringr",
  "httr",
  "jsonlite",
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
      "\n\nInstall them with:\n",
      "install.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

# ============================================================
# 2. FILE PATHS
# ============================================================

input_file <-
  "results/A10_activity_quality_controlled.csv"

assay_output_file <-
  "results/A14_assay_annotations.csv"

document_output_file <-
  "results/A14_document_annotations.csv"

enriched_output_file <-
  "results/A14_activity_provenance_enriched.csv"

duplicate_output_file <-
  "results/A14_duplicate_value_review.csv"

summary_output_file <-
  "results/A14_provenance_summary.csv"

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ", input_file,
      "\n\nCurrent working directory:\n", getwd()
    )
  )
}

# ============================================================
# 3. LOAD ACTIVITY DATA
# ============================================================

activities <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

essential_columns <- c(
  "molecule_chembl_id",
  "queried_protein",
  "standard_type",
  "standard_value_numeric",
  "standard_units",
  "assay_chembl_id",
  "document_chembl_id"
)

missing_columns <- setdiff(
  essential_columns,
  names(activities)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing essential columns in A10 data: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

# Add optional fields if absent.
optional_activity_columns <- c(
  "activity_id",
  "assay_type",
  "bao_label",
  "target_chembl_id",
  "queried_target_chembl_id",
  "target_pref_name",
  "target_organism",
  "data_validity_comment",
  "activity_comment",
  "potential_duplicate",
  "standard_relation",
  "pchembl_value_numeric",
  "evidence_class",
  "evidence_tier",
  "is_quantitative_concentration"
)

for (column_name in optional_activity_columns) {
  if (!column_name %in% names(activities)) {
    activities[[column_name]] <- NA
  }
}

# Restrict provenance curation to the six quantitative candidate molecules.
candidate_file <- "results/A10_candidate_molecule_ids.csv"

if (file.exists(candidate_file)) {
  candidate_ids <- readr::read_csv(
    candidate_file,
    show_col_types = FALSE,
    progress = FALSE
  ) |>
    dplyr::pull(molecule_chembl_id) |>
    unique()

  activities <- activities |>
    dplyr::filter(
      molecule_chembl_id %in% candidate_ids
    )
}

activities <- activities |>
  dplyr::mutate(
    assay_chembl_id = dplyr::na_if(
      stringr::str_trim(as.character(assay_chembl_id)),
      ""
    ),
    document_chembl_id = dplyr::na_if(
      stringr::str_trim(as.character(document_chembl_id)),
      ""
    ),
    standard_type = stringr::str_trim(as.character(standard_type)),
    standard_units = stringr::str_trim(as.character(standard_units)),
    standard_value_numeric = suppressWarnings(
      as.numeric(standard_value_numeric)
    )
  )

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

safe_character <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_character_)
  }

  value <- value[[1]]

  if (is.null(value) || length(value) == 0 || is.na(value)) {
    return(NA_character_)
  }

  as.character(value)
}

safe_numeric <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA_real_)
  }

  suppressWarnings(as.numeric(value[[1]]))
}

safe_logical <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA)
  }

  as.logical(value[[1]])
}

request_resource <- function(
  resource,
  identifier,
  maximum_retries = 3
) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/",
    resource,
    "/",
    identifier,
    ".json"
  )

  for (attempt in seq_len(maximum_retries)) {
    response <- tryCatch(
      httr::GET(
        url,
        httr::timeout(90),
        httr::user_agent("RESKO-eEF1A academic research")
      ),
      error = function(error) NULL
    )

    if (!is.null(response) && httr::status_code(response) == 200) {
      response_text <- httr::content(
        response,
        as = "text",
        encoding = "UTF-8"
      )

      return(
        jsonlite::fromJSON(
          response_text,
          simplifyVector = FALSE
        )
      )
    }

    Sys.sleep(attempt * 2)
  }

  NULL
}

collapse_unique <- function(value) {
  value <- sort(unique(as.character(value[!is.na(value)])))

  if (length(value) == 0) {
    return(NA_character_)
  }

  paste(value, collapse = "; ")
}

# ============================================================
# 5. RETRIEVE UNIQUE ASSAY RECORDS
# ============================================================

assay_ids <- sort(unique(activities$assay_chembl_id))
assay_ids <- assay_ids[!is.na(assay_ids)]

assay_rows <- vector("list", length(assay_ids))

for (index in seq_along(assay_ids)) {
  assay_id <- assay_ids[index]

  message(
    "Retrieving assay ", index, " of ", length(assay_ids), ": ", assay_id
  )

  assay <- request_resource("assay", assay_id)

  if (is.null(assay)) {
    assay_rows[[index]] <- tibble::tibble(
      assay_chembl_id = assay_id,
      assay_retrieved = FALSE,
      assay_type_api = NA_character_,
      assay_description = NA_character_,
      assay_organism = NA_character_,
      assay_tissue = NA_character_,
      assay_cell_type = NA_character_,
      assay_subcellular_fraction = NA_character_,
      bao_format = NA_character_,
      bao_label_api = NA_character_,
      confidence_score = NA_real_,
      target_chembl_id_api = NA_character_,
      relationship_type = NA_character_,
      document_chembl_id_api = NA_character_,
      src_assay_id = NA_character_,
      assay_category = NA_character_,
      assay_tax_id = NA_real_
    )

    next
  }

  assay_rows[[index]] <- tibble::tibble(
    assay_chembl_id = assay_id,
    assay_retrieved = TRUE,
    assay_type_api = safe_character(assay$assay_type),
    assay_description = safe_character(assay$description),
    assay_organism = safe_character(assay$assay_organism),
    assay_tissue = safe_character(assay$assay_tissue),
    assay_cell_type = safe_character(assay$assay_cell_type),
    assay_subcellular_fraction = safe_character(assay$assay_subcellular_fraction),
    bao_format = safe_character(assay$bao_format),
    bao_label_api = safe_character(assay$bao_label),
    confidence_score = safe_numeric(assay$confidence_score),
    target_chembl_id_api = safe_character(assay$target_chembl_id),
    relationship_type = safe_character(assay$relationship_type),
    document_chembl_id_api = safe_character(assay$document_chembl_id),
    src_assay_id = safe_character(assay$src_assay_id),
    assay_category = safe_character(assay$assay_category),
    assay_tax_id = safe_numeric(assay$assay_tax_id)
  )

  Sys.sleep(0.2)
}

assay_annotations <- dplyr::bind_rows(assay_rows)

# ============================================================
# 6. IDENTIFY AND RETRIEVE DOCUMENT RECORDS
# ============================================================

# Prefer document IDs from activities, while also including any
# document IDs returned by the assay endpoint.
document_ids <- sort(
  unique(
    c(
      activities$document_chembl_id,
      assay_annotations$document_chembl_id_api
    )
  )
)

document_ids <- document_ids[!is.na(document_ids)]

document_rows <- vector("list", length(document_ids))

for (index in seq_along(document_ids)) {
  document_id <- document_ids[index]

  message(
    "Retrieving document ", index, " of ", length(document_ids), ": ", document_id
  )

  document <- request_resource("document", document_id)

  if (is.null(document)) {
    document_rows[[index]] <- tibble::tibble(
      document_chembl_id = document_id,
      document_retrieved = FALSE,
      document_type = NA_character_,
      title = NA_character_,
      abstract = NA_character_,
      authors = NA_character_,
      journal = NA_character_,
      journal_full_title = NA_character_,
      year = NA_real_,
      volume = NA_character_,
      issue = NA_character_,
      first_page = NA_character_,
      last_page = NA_character_,
      doi = NA_character_,
      pubmed_id = NA_character_,
      patent_id = NA_character_
    )

    next
  }

  document_rows[[index]] <- tibble::tibble(
    document_chembl_id = document_id,
    document_retrieved = TRUE,
    document_type = safe_character(document$doc_type),
    title = safe_character(document$title),
    abstract = safe_character(document$abstract),
    authors = safe_character(document$authors),
    journal = safe_character(document$journal),
    journal_full_title = safe_character(document$journal_full_title),
    year = safe_numeric(document$year),
    volume = safe_character(document$volume),
    issue = safe_character(document$issue),
    first_page = safe_character(document$first_page),
    last_page = safe_character(document$last_page),
    doi = safe_character(document$doi),
    pubmed_id = safe_character(document$pubmed_id),
    patent_id = safe_character(document$patent_id)
  )

  Sys.sleep(0.2)
}

document_annotations <- dplyr::bind_rows(document_rows)

# ============================================================
# 7. JOIN PROVENANCE TO ACTIVITY RECORDS
# ============================================================

activities_enriched <- activities |>
  dplyr::left_join(
    assay_annotations,
    by = "assay_chembl_id"
  ) |>
  dplyr::mutate(
    final_document_chembl_id = dplyr::coalesce(
      document_chembl_id,
      document_chembl_id_api
    ),
    final_assay_type = dplyr::coalesce(
      assay_type_api,
      as.character(assay_type)
    ),
    final_bao_label = dplyr::coalesce(
      bao_label_api,
      as.character(bao_label)
    ),
    final_target_chembl_id = dplyr::coalesce(
      target_chembl_id_api,
      as.character(target_chembl_id),
      as.character(queried_target_chembl_id)
    ),
    high_confidence_target_assignment =
      !is.na(confidence_score) & confidence_score >= 8,
    single_protein_high_confidence =
      !is.na(confidence_score) & confidence_score == 9,
    direct_binding_assay =
      final_assay_type == "B",
    functional_assay =
      final_assay_type == "F"
  ) |>
  dplyr::left_join(
    document_annotations,
    by = c(
      "final_document_chembl_id" = "document_chembl_id"
    )
  )

# ============================================================
# 8. REVIEW IDENTICAL VALUES AND POSSIBLE DUPLICATION
#
# Groups together records with identical molecule, protein,
# units and numeric value, then checks whether multiple activity
# types share assays or documents.
# ============================================================

duplicate_value_review <- activities_enriched |>
  dplyr::filter(
    !is.na(standard_value_numeric),
    standard_type %in% c("IC50", "Kd", "ED50")
  ) |>
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
    unique_assay_count = dplyr::n_distinct(
      assay_chembl_id[!is.na(assay_chembl_id)]
    ),
    document_ids = collapse_unique(final_document_chembl_id),
    unique_document_count = dplyr::n_distinct(
      final_document_chembl_id[!is.na(final_document_chembl_id)]
    ),
    assay_types = collapse_unique(final_assay_type),
    assay_descriptions = collapse_unique(assay_description),
    confidence_scores = collapse_unique(confidence_score),
    relationship_types = collapse_unique(relationship_type),
    dois = collapse_unique(doi),
    pubmed_ids = collapse_unique(pubmed_id),
    has_same_value_across_activity_types =
      dplyr::n_distinct(standard_type) > 1,
    likely_same_assay_reannotation =
      dplyr::n_distinct(standard_type) > 1 &
      dplyr::n_distinct(assay_chembl_id[!is.na(assay_chembl_id)]) <= 1,
    likely_same_document_evidence =
      dplyr::n_distinct(standard_type) > 1 &
      dplyr::n_distinct(final_document_chembl_id[!is.na(final_document_chembl_id)]) <= 1,
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

# ============================================================
# 9. CREATE PROVENANCE SUMMARY
# ============================================================

provenance_summary <- dplyr::bind_rows(
  tibble::tibble(
    metric = "Candidate activity records curated",
    count = nrow(activities_enriched)
  ),
  tibble::tibble(
    metric = "Unique assays",
    count = dplyr::n_distinct(
      activities_enriched$assay_chembl_id[
        !is.na(activities_enriched$assay_chembl_id)
      ]
    )
  ),
  tibble::tibble(
    metric = "Assays successfully retrieved",
    count = sum(assay_annotations$assay_retrieved, na.rm = TRUE)
  ),
  tibble::tibble(
    metric = "Unique source documents",
    count = dplyr::n_distinct(
      activities_enriched$final_document_chembl_id[
        !is.na(activities_enriched$final_document_chembl_id)
      ]
    )
  ),
  tibble::tibble(
    metric = "Documents successfully retrieved",
    count = sum(document_annotations$document_retrieved, na.rm = TRUE)
  ),
  tibble::tibble(
    metric = "Binding-assay records",
    count = sum(activities_enriched$direct_binding_assay, na.rm = TRUE)
  ),
  tibble::tibble(
    metric = "Functional-assay records",
    count = sum(activities_enriched$functional_assay, na.rm = TRUE)
  ),
  tibble::tibble(
    metric = "High-confidence target assignments (score >=8)",
    count = sum(
      activities_enriched$high_confidence_target_assignment,
      na.rm = TRUE
    )
  ),
  tibble::tibble(
    metric = "Single-protein target assignments (score 9)",
    count = sum(
      activities_enriched$single_protein_high_confidence,
      na.rm = TRUE
    )
  ),
  tibble::tibble(
    metric = "Identical values across multiple activity types",
    count = sum(
      duplicate_value_review$has_same_value_across_activity_types,
      na.rm = TRUE
    )
  ),
  tibble::tibble(
    metric = "Likely same-assay reannotations",
    count = sum(
      duplicate_value_review$likely_same_assay_reannotation,
      na.rm = TRUE
    )
  )
)

# ============================================================
# 10. WRITE AND VERIFY OUTPUTS
# ============================================================

readr::write_csv(
  assay_annotations,
  assay_output_file,
  na = ""
)

readr::write_csv(
  document_annotations,
  document_output_file,
  na = ""
)

readr::write_csv(
  activities_enriched,
  enriched_output_file,
  na = ""
)

readr::write_csv(
  duplicate_value_review,
  duplicate_output_file,
  na = ""
)

readr::write_csv(
  provenance_summary,
  summary_output_file,
  na = ""
)

output_files <- c(
  assay_output_file,
  document_output_file,
  enriched_output_file,
  duplicate_output_file,
  summary_output_file
)

file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(file.info(output_files)$size)
)

if (!all(file_check$exists)) {
  stop("At least one A14 output file was not created.")
}

if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) {
  stop("At least one A14 output file is empty or invalid.")
}

# ============================================================
# 11. PRINT KEY RESULTS
# ============================================================

cat("\n")
cat("A14 assay and document provenance curation completed.\n")
cat("----------------------------------------------------\n")

print(provenance_summary, n = Inf)

cat("\nAssay confidence scores:\n")
print(
  assay_annotations |>
    dplyr::count(confidence_score, sort = TRUE),
  n = Inf
)

cat("\nIdentical-value groups requiring review:\n")
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
