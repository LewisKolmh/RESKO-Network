# ============================================================
# A12 CREATE COMPOUND-PROTEIN GRAPH RELATIONSHIPS
#
# Purpose:
# - Convert quality-controlled ChEMBL evidence into graph edges
# - Keep IC50, Kd and ED50 evidence distinct
# - Add molecule names from A11 annotations
# - Create one edge per compound-protein-activity-type group
# - Create a simplified one-edge-per-compound-protein table
#
# Inputs:
#   results/A10_compound_protein_summary.csv
#   results/A11_candidate_molecule_annotations.csv
#   results/nodes_proteins.csv
#
# Outputs:
#   results/A12_edges_compound_protein_evidence.csv
#   results/edges_compound_protein.csv
#   results/A12_compound_protein_edge_summary.csv
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
# 2. FILE PATHS
# ============================================================

compound_summary_file <-
  "results/A10_compound_protein_summary.csv"

annotation_file <-
  "results/A11_candidate_molecule_annotations.csv"

protein_nodes_file <-
  "results/nodes_proteins.csv"

evidence_edges_file <-
  "results/A12_edges_compound_protein_evidence.csv"

simple_edges_file <-
  "results/edges_compound_protein.csv"

edge_summary_file <-
  "results/A12_compound_protein_edge_summary.csv"

input_files <- c(
  compound_summary_file,
  annotation_file,
  protein_nodes_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n"),
      "\n\nCurrent working directory:\n",
      getwd(),
      "\n\nComplete A10 and A11 before running A12."
    )
  )
}

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# ============================================================
# 3. LOAD DATA
# ============================================================

compound_summary <- readr::read_csv(
  compound_summary_file,
  show_col_types = FALSE,
  progress = FALSE
)

annotations <- readr::read_csv(
  annotation_file,
  show_col_types = FALSE,
  progress = FALSE
)

protein_nodes <- readr::read_csv(
  protein_nodes_file,
  show_col_types = FALSE,
  progress = FALSE
)

# ============================================================
# 4. VALIDATE COLUMNS
# ============================================================

check_columns <- function(data, required_columns, file_name) {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing columns in ",
        file_name,
        ": ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}

check_columns(
  compound_summary,
  c(
    "molecule_chembl_id",
    "queried_protein",
    "queried_target_chembl_id",
    "standard_type",
    "standard_units",
    "evidence_class",
    "evidence_tier",
    "measurement_count",
    "minimum_standard_value",
    "median_standard_value",
    "maximum_standard_value",
    "maximum_pchembl_value",
    "unique_assay_count",
    "validity_warning_count",
    "potential_duplicate_count"
  ),
  compound_summary_file
)

check_columns(
  annotations,
  c(
    "molecule_chembl_id",
    "display_name",
    "molecule_type",
    "max_phase",
    "development_status"
  ),
  annotation_file
)

check_columns(
  protein_nodes,
  "protein",
  protein_nodes_file
)

# ============================================================
# 5. CLEAN INPUT DATA
# ============================================================

compound_summary <- compound_summary |>
  dplyr::filter(
    !is.na(molecule_chembl_id),
    molecule_chembl_id != "",
    !is.na(queried_protein),
    queried_protein != ""
  ) |>
  dplyr::mutate(
    measurement_count = suppressWarnings(as.integer(measurement_count)),
    minimum_standard_value = suppressWarnings(as.numeric(minimum_standard_value)),
    median_standard_value = suppressWarnings(as.numeric(median_standard_value)),
    maximum_standard_value = suppressWarnings(as.numeric(maximum_standard_value)),
    maximum_pchembl_value = suppressWarnings(as.numeric(maximum_pchembl_value)),
    unique_assay_count = suppressWarnings(as.integer(unique_assay_count)),
    validity_warning_count = suppressWarnings(as.integer(validity_warning_count)),
    potential_duplicate_count = suppressWarnings(as.integer(potential_duplicate_count))
  ) |>
  dplyr::distinct(
    molecule_chembl_id,
    queried_protein,
    standard_type,
    standard_units,
    .keep_all = TRUE
  )

annotations <- annotations |>
  dplyr::filter(
    !is.na(molecule_chembl_id),
    molecule_chembl_id != ""
  ) |>
  dplyr::distinct(
    molecule_chembl_id,
    .keep_all = TRUE
  )

protein_nodes <- protein_nodes |>
  dplyr::filter(
    !is.na(protein),
    protein != ""
  ) |>
  dplyr::distinct(protein)

# ============================================================
# 6. CHECK THAT ALL TARGET PROTEINS EXIST IN THE GRAPH
# ============================================================

missing_graph_proteins <- setdiff(
  unique(compound_summary$queried_protein),
  protein_nodes$protein
)

if (length(missing_graph_proteins) > 0) {
  stop(
    paste0(
      "The following compound targets are missing from nodes_proteins.csv:\n",
      paste(missing_graph_proteins, collapse = "\n")
    )
  )
}

# ============================================================
# 7. ASSIGN GRAPH RELATIONSHIP TYPES
#
# IC50 is represented as inhibitory evidence.
# Kd is represented as binding evidence.
# ED50 is represented as functional-response evidence.
# ============================================================

compound_evidence_edges <- compound_summary |>
  dplyr::left_join(
    annotations |>
      dplyr::select(
        molecule_chembl_id,
        display_name,
        molecule_type,
        max_phase,
        development_status
      ),
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    display_name = dplyr::coalesce(
      display_name,
      molecule_chembl_id
    ),
    relationship = dplyr::case_when(
      standard_type == "IC50" ~ "INHIBITS",
      standard_type == "Kd" ~ "BINDS_TO",
      standard_type == "ED50" ~ "HAS_FUNCTIONAL_ACTIVITY_AGAINST",
      TRUE ~ "HAS_ACTIVITY_AGAINST"
    ),
    evidence_rank = dplyr::case_when(
      standard_type == "IC50" ~ 1L,
      standard_type == "Kd" ~ 2L,
      standard_type == "ED50" ~ 3L,
      TRUE ~ 4L
    ),
    potency_band = dplyr::case_when(
      is.na(minimum_standard_value) ~ "Unknown",
      minimum_standard_value <= 10 ~ "<=10 nM",
      minimum_standard_value <= 100 ~ ">10 to 100 nM",
      minimum_standard_value <= 1000 ~ ">100 to 1000 nM",
      minimum_standard_value <= 10000 ~ ">1 to 10 uM",
      TRUE ~ ">10 uM"
    ),
    source = molecule_chembl_id,
    target = queried_protein,
    source_node_id = paste0("compound::", molecule_chembl_id),
    target_node_id = paste0("protein::", queried_protein)
  ) |>
  dplyr::transmute(
    source,
    target,
    source_node_id,
    target_node_id,
    relationship,
    molecule_chembl_id,
    compound_name = display_name,
    protein = queried_protein,
    target_chembl_id = queried_target_chembl_id,
    standard_type,
    standard_units,
    evidence_class,
    evidence_tier,
    evidence_rank,
    potency_band,
    measurement_count,
    minimum_standard_value,
    median_standard_value,
    maximum_standard_value,
    maximum_pchembl_value,
    unique_assay_count,
    validity_warning_count,
    potential_duplicate_count,
    molecule_type,
    max_phase,
    development_status
  ) |>
  dplyr::arrange(
    evidence_rank,
    minimum_standard_value,
    compound_name,
    protein
  )

# ============================================================
# 8. CREATE A SIMPLIFIED ONE-EDGE-PER-COMPOUND-PROTEIN TABLE
#
# If a compound has several evidence types against one protein,
# the strongest evidence relationship is used for display while
# every activity type is retained in the supporting properties.
# ============================================================

simple_edges <- compound_evidence_edges |>
  dplyr::group_by(
    molecule_chembl_id,
    compound_name,
    protein,
    source_node_id,
    target_node_id
  ) |>
  dplyr::arrange(
    evidence_rank,
    minimum_standard_value,
    .by_group = TRUE
  ) |>
  dplyr::summarise(
    source = dplyr::first(molecule_chembl_id),
    target = dplyr::first(protein),
    relationship = dplyr::first(relationship),
    strongest_evidence_class = dplyr::first(evidence_class),
    strongest_evidence_tier = dplyr::first(evidence_tier),
    activity_types = paste(
      sort(unique(standard_type)),
      collapse = "; "
    ),
    evidence_classes = paste(
      sort(unique(evidence_class)),
      collapse = "; "
    ),
    minimum_activity_nM = {
      values <- minimum_standard_value[is.finite(minimum_standard_value)]
      if (length(values) == 0) NA_real_ else min(values)
    },
    median_activity_nM = {
      values <- median_standard_value[is.finite(median_standard_value)]
      if (length(values) == 0) NA_real_ else stats::median(values)
    },
    maximum_pchembl_value = {
      values <- maximum_pchembl_value[is.finite(maximum_pchembl_value)]
      if (length(values) == 0) NA_real_ else max(values)
    },
    total_measurements = sum(measurement_count, na.rm = TRUE),
    total_assays = sum(unique_assay_count, na.rm = TRUE),
    validity_warning_count = sum(validity_warning_count, na.rm = TRUE),
    potential_duplicate_count = sum(potential_duplicate_count, na.rm = TRUE),
    molecule_type = dplyr::first(molecule_type),
    max_phase = dplyr::first(max_phase),
    development_status = dplyr::first(development_status),
    .groups = "drop"
  ) |>
  dplyr::arrange(
    relationship,
    minimum_activity_nM,
    compound_name,
    protein
  )

# ============================================================
# 9. CREATE EDGE SUMMARY
# ============================================================

edge_summary <- dplyr::bind_rows(
  tibble::tibble(
    metric = "Evidence-level compound-protein edges",
    count = nrow(compound_evidence_edges)
  ),
  tibble::tibble(
    metric = "Simplified compound-protein edges",
    count = nrow(simple_edges)
  ),
  tibble::tibble(
    metric = "Unique compounds",
    count = dplyr::n_distinct(simple_edges$molecule_chembl_id)
  ),
  tibble::tibble(
    metric = "Unique target proteins",
    count = dplyr::n_distinct(simple_edges$protein)
  ),
  simple_edges |>
    dplyr::count(
      relationship,
      name = "count"
    ) |>
    dplyr::transmute(
      metric = paste0("Relationship: ", relationship),
      count
    )
)

# ============================================================
# 10. WRITE OUTPUTS
# ============================================================

readr::write_csv(
  compound_evidence_edges,
  evidence_edges_file,
  na = ""
)

readr::write_csv(
  simple_edges,
  simple_edges_file,
  na = ""
)

readr::write_csv(
  edge_summary,
  edge_summary_file,
  na = ""
)

output_files <- c(
  evidence_edges_file,
  simple_edges_file,
  edge_summary_file
)

file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(file.info(output_files)$size)
)

if (!all(file_check$exists)) {
  stop("At least one A12 output file was not created.")
}

if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) {
  stop("At least one A12 output file is empty or invalid.")
}

# ============================================================
# 11. PRINT SUMMARY
# ============================================================

cat("\n")
cat("A12 compound-protein edge creation completed.\n")
cat("---------------------------------------------\n")

print(edge_summary, n = Inf)

cat("\nSimplified compound-protein relationships:\n")

print(
  simple_edges |>
    dplyr::select(
      molecule_chembl_id,
      compound_name,
      protein,
      relationship,
      activity_types,
      minimum_activity_nM,
      total_measurements,
      development_status
    ),
  n = Inf,
  width = Inf
)

cat("\nOutput verification:\n")
print(file_check, n = Inf)
