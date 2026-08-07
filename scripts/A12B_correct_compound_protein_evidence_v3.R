# ============================================================
# A12B V3: CORRECT COMPOUND-PROTEIN EVIDENCE
# ============================================================

required_packages <- c("readr", "dplyr", "tibble")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(paste0(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))"
  ))
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

activity_file <- file.path(results_dir, "A14_activity_provenance_enriched.csv")
duplicate_file <- file.path(results_dir, "A14_duplicate_value_review.csv")
annotation_file <- file.path(results_dir, "A11_candidate_molecule_annotations.csv")
protein_file <- file.path(results_dir, "nodes_proteins.csv")

corrected_evidence_file <- file.path(results_dir, "A12B_compound_protein_evidence_corrected.csv")
corrected_edge_file <- file.path(results_dir, "edges_compound_protein_corrected.csv")
summary_file <- file.path(results_dir, "A12B_correction_summary.csv")

input_files <- c(activity_file, duplicate_file, annotation_file, protein_file)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(paste0(
    "Missing input files:\n", paste(missing_files, collapse = "\n"),
    "\n\nCurrent project root: ", project_root
  ))
}

activities <- readr::read_csv(activity_file, show_col_types = FALSE, progress = FALSE)
duplicate_review <- readr::read_csv(duplicate_file, show_col_types = FALSE, progress = FALSE)
annotations <- readr::read_csv(annotation_file, show_col_types = FALSE, progress = FALSE)
proteins <- readr::read_csv(protein_file, show_col_types = FALSE, progress = FALSE)

required_activity_columns <- c(
  "molecule_chembl_id", "queried_protein", "standard_type",
  "standard_value_numeric", "standard_units", "assay_chembl_id",
  "final_document_chembl_id", "confidence_score", "final_assay_type",
  "evidence_class", "evidence_tier"
)
missing_columns <- setdiff(required_activity_columns, names(activities))
if (length(missing_columns) > 0) {
  stop(paste0("Missing A14 columns: ", paste(missing_columns, collapse = ", ")))
}

optional_columns <- c(
  "activity_id", "standard_relation", "pchembl_value_numeric",
  "data_validity_comment", "potential_duplicate", "doi", "pubmed_id",
  "relationship_type", "bao_format", "final_bao_label",
  "assay_description", "assay_description.x", "assay_description.y",
  "description", "description.x", "description.y"
)
for (column_name in optional_columns) {
  if (!column_name %in% names(activities)) activities[[column_name]] <- NA
}

# Resolve assay description across possible column names created by joins.
assay_description_candidates <- c(
  "assay_description",
  "assay_description.x",
  "assay_description.y",
  "description",
  "description.x",
  "description.y"
)

activity_assay_description <- rep(NA_character_, nrow(activities))

for (candidate_column in assay_description_candidates) {
  if (candidate_column %in% names(activities)) {
    candidate_values <- as.character(activities[[candidate_column]])
    use_values <- (
      is.na(activity_assay_description) |
      activity_assay_description == ""
    ) & !is.na(candidate_values) & candidate_values != ""
    activity_assay_description[use_values] <- candidate_values[use_values]
  }
}

activities$assay_description_resolved <- activity_assay_description

required_annotation_columns <- c(
  "molecule_chembl_id", "display_name", "molecule_type",
  "max_phase", "development_status"
)
missing_columns <- setdiff(required_annotation_columns, names(annotations))
if (length(missing_columns) > 0) {
  stop(paste0("Missing A11 columns: ", paste(missing_columns, collapse = ", ")))
}
if (!"protein" %in% names(proteins)) stop("nodes_proteins.csv lacks protein column")

collapse_unique <- function(x) {
  x <- sort(unique(as.character(x[!is.na(x) & as.character(x) != ""])))
  if (length(x) == 0) NA_character_ else paste(x, collapse = "; ")
}
safe_min <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else min(x)
}
safe_median <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else stats::median(x)
}
safe_max <- function(x) {
  x <- suppressWarnings(as.numeric(x)); x <- x[is.finite(x)]
  if (length(x) == 0) NA_real_ else max(x)
}
count_distinct_nonmissing <- function(x) {
  x <- x[!is.na(x) & as.character(x) != ""]
  dplyr::n_distinct(x)
}

# Keys identifying the A14-confirmed Kd/ED50 duplicate groups.
duplicate_keys <- paste(
  duplicate_review$molecule_chembl_id,
  duplicate_review$queried_protein,
  format(as.numeric(duplicate_review$standard_value_numeric), scientific = FALSE, trim = TRUE),
  duplicate_review$standard_units,
  sep = "||"
)

candidate_ids <- unique(annotations$molecule_chembl_id)

activities_filtered <- activities |>
  dplyr::mutate(
    standard_value_numeric = suppressWarnings(as.numeric(standard_value_numeric)),
    confidence_score = suppressWarnings(as.numeric(confidence_score)),
    pchembl_value_numeric = suppressWarnings(as.numeric(pchembl_value_numeric)),
    duplicate_key = paste(
      molecule_chembl_id,
      queried_protein,
      format(standard_value_numeric, scientific = FALSE, trim = TRUE),
      standard_units,
      sep = "||"
    ),
    duplicate_kd_ed50 = duplicate_key %in% duplicate_keys,
    corrected_standard_type = dplyr::case_when(
      duplicate_kd_ed50 & standard_type %in% c("Kd", "ED50") ~ "Kd",
      TRUE ~ as.character(standard_type)
    ),
    corrected_relationship = dplyr::case_when(
      corrected_standard_type == "IC50" ~ "HAS_INHIBITORY_ACTIVITY_AGAINST",
      corrected_standard_type == "Kd" ~ "BINDS_TO",
      TRUE ~ "HAS_ACTIVITY_AGAINST"
    ),
    corrected_evidence_class = dplyr::case_when(
      corrected_standard_type == "IC50" ~ "Inhibitory activity from binding assay",
      corrected_standard_type == "Kd" ~ "Binding affinity",
      TRUE ~ as.character(evidence_class)
    )
  ) |>
  dplyr::filter(
    molecule_chembl_id %in% candidate_ids,
    standard_type %in% c("IC50", "Kd", "ED50"),
    is.finite(standard_value_numeric)
  )

corrected_evidence <- activities_filtered |>
  dplyr::group_by(
    molecule_chembl_id,
    queried_protein,
    corrected_standard_type,
    corrected_relationship,
    assay_chembl_id,
    final_document_chembl_id,
    standard_value_numeric,
    standard_units
  ) |>
  dplyr::summarise(
    original_activity_types = collapse_unique(standard_type),
    original_activity_ids = collapse_unique(activity_id),
    corrected_evidence_class = dplyr::first(corrected_evidence_class),
    evidence_tiers = collapse_unique(evidence_tier),
    standard_relations = collapse_unique(standard_relation),
    maximum_pchembl_value = safe_max(pchembl_value_numeric),
    confidence_score = safe_max(confidence_score),
    assay_type = dplyr::first(as.character(final_assay_type)),
    assay_description = dplyr::first(as.character(assay_description_resolved)),
    assay_relationship_type = collapse_unique(relationship_type),
    bao_format = collapse_unique(bao_format),
    bao_label = collapse_unique(final_bao_label),
    doi = collapse_unique(doi),
    pubmed_id = collapse_unique(pubmed_id),
    data_validity_comments = collapse_unique(data_validity_comment),
    potential_duplicate_flags = collapse_unique(potential_duplicate),
    raw_record_count = dplyr::n(),
    duplicate_kd_ed50_collapsed = any(duplicate_kd_ed50),
    .groups = "drop"
  ) |>
  dplyr::left_join(
    annotations |>
      dplyr::select(
        molecule_chembl_id, display_name, molecule_type,
        max_phase, development_status
      ) |>
      dplyr::distinct(molecule_chembl_id, .keep_all = TRUE),
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    display_name = dplyr::coalesce(display_name, molecule_chembl_id),
    source = molecule_chembl_id,
    target = queried_protein,
    source_node_id = paste0("compound::", molecule_chembl_id),
    target_node_id = paste0("protein::", queried_protein)
  ) |>
  dplyr::ungroup()

corrected_edges <- corrected_evidence |>
  dplyr::mutate(
    evidence_rank = dplyr::case_when(
      corrected_relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ 1L,
      corrected_relationship == "BINDS_TO" ~ 2L,
      TRUE ~ 3L
    )
  ) |>
  dplyr::group_by(
    molecule_chembl_id,
    display_name,
    queried_protein,
    source_node_id,
    target_node_id
  ) |>
  dplyr::arrange(evidence_rank, standard_value_numeric, .by_group = TRUE) |>
  dplyr::summarise(
    source = dplyr::first(molecule_chembl_id),
    target = dplyr::first(queried_protein),
    relationship = dplyr::first(corrected_relationship),
    primary_evidence_class = dplyr::first(corrected_evidence_class),
    corrected_activity_types = collapse_unique(corrected_standard_type),
    original_activity_types = collapse_unique(original_activity_types),
    minimum_activity_nM = safe_min(standard_value_numeric),
    median_activity_nM = safe_median(standard_value_numeric),
    maximum_pchembl_value = safe_max(maximum_pchembl_value),
    unique_assay_count = count_distinct_nonmissing(assay_chembl_id),
    assay_ids = collapse_unique(assay_chembl_id),
    unique_document_count = count_distinct_nonmissing(final_document_chembl_id),
    document_ids = collapse_unique(final_document_chembl_id),
    confidence_score = safe_max(confidence_score),
    assay_types = collapse_unique(assay_type),
    assay_descriptions = collapse_unique(assay_description),
    dois = collapse_unique(doi),
    pubmed_ids = collapse_unique(pubmed_id),
    raw_record_count = sum(raw_record_count, na.rm = TRUE),
    duplicate_kd_ed50_collapsed = any(duplicate_kd_ed50_collapsed),
    molecule_type = dplyr::first(molecule_type),
    max_phase = dplyr::first(max_phase),
    development_status = dplyr::first(development_status),
    .groups = "drop"
  ) |>
  dplyr::arrange(relationship, minimum_activity_nM, display_name, queried_protein)

missing_proteins <- setdiff(unique(corrected_edges$queried_protein), proteins$protein)
if (length(missing_proteins) > 0) {
  stop(paste0(
    "Corrected edges reference proteins absent from nodes_proteins.csv:\n",
    paste(missing_proteins, collapse = "\n")
  ))
}

summary_table <- dplyr::bind_rows(
  tibble::tibble(
    metric = "Quantitative records before correction",
    count = nrow(activities_filtered)
  ),
  tibble::tibble(
    metric = "Corrected evidence records",
    count = nrow(corrected_evidence)
  ),
  tibble::tibble(
    metric = "Corrected compound-protein graph edges",
    count = nrow(corrected_edges)
  ),
  tibble::tibble(
    metric = "Kd/ED50 duplicate groups collapsed",
    count = sum(corrected_evidence$duplicate_kd_ed50_collapsed, na.rm = TRUE)
  ),
  corrected_edges |>
    dplyr::count(relationship, name = "count") |>
    dplyr::transmute(metric = paste0("Relationship: ", relationship), count)
)

# The result tables are plain data frames before export.
corrected_evidence_export <- as.data.frame(dplyr::ungroup(corrected_evidence), stringsAsFactors = FALSE)
corrected_edges_export <- as.data.frame(dplyr::ungroup(corrected_edges), stringsAsFactors = FALSE)
summary_export <- as.data.frame(dplyr::ungroup(summary_table), stringsAsFactors = FALSE)

# Use absolute paths and base write.csv for reliable Windows export.
utils::write.csv(
  corrected_evidence_export,
  corrected_evidence_file,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  corrected_edges_export,
  corrected_edge_file,
  row.names = FALSE,
  na = ""
)
utils::write.csv(
  summary_export,
  summary_file,
  row.names = FALSE,
  na = ""
)

output_files <- c(corrected_evidence_file, corrected_edge_file, summary_file)
file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(file.info(output_files)$size)
)

if (!all(file_check$exists)) {
  stop(paste0(
    "One or more A12B output files were not created.\n",
    paste(capture.output(print(file_check)), collapse = "\n")
  ))
}
if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) {
  stop("One or more A12B output files are empty or invalid.")
}

cat("\n")
cat("A12B V3 evidence correction completed successfully.\n")
cat("--------------------------------------------------\n")
cat("Project root: ", project_root, "\n\n", sep = "")
print(summary_table, n = Inf)
cat("\nOutput verification:\n")
print(file_check, n = Inf)
