#!/usr/bin/env Rscript

# ============================================================
# RESKO A17C: Screen SIDER drugs by chemical similarity
# ============================================================
#
# Run from the RESKO project root:
#
# conda activate resko-a16
#
# RETICULATE_PYTHON="$(which python)" \
# Rscript scripts/A17C_screen_sider_drugs_by_chemical_similarity.R
#
# Inputs:
#
# results/A17B_sider_identifier_crosswalk.csv
# results/A17A_sider_side_effects_pt.csv
# results/A16_structures_combined.csv
#
# Outputs:
#
# results/A17C_sider_query_similarity_all.csv
# results/A17C_sider_nearest_query.csv
# results/A17C_query_nearest_sider_drugs.csv
# results/A17C_sider_similarity_candidates.csv
# results/A17C_query_similarity_summary.csv
# results/A17C_similarity_plot.png
# results/A17C_sider_similarity_report.html
# results/A17C_similarity_summary.csv
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# ============================================================
# 1. Package validation
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tibble",
  "ggplot2",
  "reticulate",
  "base64enc"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    quietly = TRUE,
    FUN.VALUE = logical(1)
  )
]

if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ",
    paste(missing_packages, collapse = ", "),
    call. = FALSE
  )
}

# ============================================================
# 2. Project paths
# ============================================================

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

scripts_dir <- file.path(
  project_root,
  "scripts"
)

results_dir <- file.path(
  project_root,
  "results"
)

if (!dir.exists(scripts_dir)) {
  stop(
    "The scripts directory was not found. ",
    "Run A17C from the RESKO project root.",
    call. = FALSE
  )
}

if (!dir.exists(results_dir)) {
  stop(
    "The results directory was not found. ",
    "Run A17C from the RESKO project root.",
    call. = FALSE
  )
}

input_sider <- file.path(
  results_dir,
  "A17B_sider_identifier_crosswalk.csv"
)

input_side_effects <- file.path(
  results_dir,
  "A17A_sider_side_effects_pt.csv"
)

input_queries <- file.path(
  results_dir,
  "A16_structures_combined.csv"
)

output_all_similarity <- file.path(
  results_dir,
  "A17C_sider_query_similarity_all.csv"
)

output_sider_nearest <- file.path(
  results_dir,
  "A17C_sider_nearest_query.csv"
)

output_query_nearest <- file.path(
  results_dir,
  "A17C_query_nearest_sider_drugs.csv"
)

output_candidates <- file.path(
  results_dir,
  "A17C_sider_similarity_candidates.csv"
)

output_query_summary <- file.path(
  results_dir,
  "A17C_query_similarity_summary.csv"
)

output_plot <- file.path(
  results_dir,
  "A17C_similarity_plot.png"
)

output_report <- file.path(
  results_dir,
  "A17C_sider_similarity_report.html"
)

output_summary <- file.path(
  results_dir,
  "A17C_similarity_summary.csv"
)

required_inputs <- c(
  input_sider,
  input_side_effects,
  input_queries
)

for (path in required_inputs) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Required input is missing or empty: ",
      path,
      call. = FALSE
    )
  }
}

# ============================================================
# 3. Helper functions
# ============================================================

clean_text <- function(x) {
  output <- as.character(x)

  output[
    is.na(output) |
      trimws(output) == ""
  ] <- NA_character_

  output
}

assign_similarity_band <- function(x) {
  dplyr::case_when(
    x >= 0.70 ~ "close_structural_neighbour",
    x >= 0.50 ~ "moderate_similarity",
    x >= 0.30 ~ "weak_or_partial_similarity",
    TRUE ~ "structurally_distant_under_this_fingerprint"
  )
}

html_escape <- function(x) {
  output <- as.character(x)

  output <- gsub(
    "&",
    "&amp;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    "<",
    "&lt;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    ">",
    "&gt;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    "\"",
    "&quot;",
    output,
    fixed = TRUE
  )

  output
}

backup_existing_file <- function(path) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  timestamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )

  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)

  backup_path <- paste0(
    stem,
    "_previous_",
    timestamp,
    if (extension == "") {
      ""
    } else {
      paste0(".", extension)
    }
  )

  copied <- file.copy(
    from = path,
    to = backup_path,
    overwrite = FALSE
  )

  if (!copied) {
    stop(
      "Could not preserve previous output: ",
      path,
      call. = FALSE
    )
  }

  if (
    !file.exists(backup_path) ||
    is.na(file.info(backup_path)$size) ||
    file.info(backup_path)$size <= 0L
  ) {
    stop(
      "Backup verification failed: ",
      backup_path,
      call. = FALSE
    )
  }

  invisible(backup_path)
}

safe_write_csv <- function(data, path) {
  temporary_path <- paste0(
    path,
    ".tmp"
  )

  if (file.exists(temporary_path)) {
    unlink(temporary_path)
  }

  data_to_write <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  list_columns <- names(data_to_write)[
    vapply(
      data_to_write,
      is.list,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(list_columns) > 0L) {
    stop(
      "Refusing to write list-column(s): ",
      paste(list_columns, collapse = ", "),
      call. = FALSE
    )
  }

  readr::write_csv(
    data_to_write,
    temporary_path,
    na = ""
  )

  if (
    !file.exists(temporary_path) ||
    is.na(file.info(temporary_path)$size) ||
    file.info(temporary_path)$size <= 0L
  ) {
    stop(
      "Failed to create temporary output: ",
      temporary_path,
      call. = FALSE
    )
  }

  backup_existing_file(path)

  moved <- file.rename(
    from = temporary_path,
    to = path
  )

  if (!moved) {
    stop(
      "Could not move temporary output into place: ",
      path,
      call. = FALSE
    )
  }

  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Output verification failed: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}

safe_save_plot <- function(
  plot,
  path,
  width,
  height,
  dpi = 220
) {
  temporary_path <- tempfile(
    pattern = "A17C_plot_",
    tmpdir = results_dir,
    fileext = ".png"
  )

  ggplot2::ggsave(
    filename = temporary_path,
    plot = plot,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )

  if (
    !file.exists(temporary_path) ||
    is.na(file.info(temporary_path)$size) ||
    file.info(temporary_path)$size <= 0L
  ) {
    stop(
      "Failed to create plot: ",
      path,
      call. = FALSE
    )
  }

  backup_existing_file(path)

  moved <- file.rename(
    from = temporary_path,
    to = path
  )

  if (!moved) {
    stop(
      "Could not move plot into place: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}

table_to_html <- function(
  data,
  digits = 3L
) {
  table_data <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  for (column_name in names(table_data)) {
    if (is.numeric(table_data[[column_name]])) {
      table_data[[column_name]] <- ifelse(
        is.na(table_data[[column_name]]),
        "",
        format(
          round(
            table_data[[column_name]],
            digits
          ),
          nsmall = digits,
          trim = TRUE
        )
      )
    } else {
      table_data[[column_name]] <- ifelse(
        is.na(table_data[[column_name]]),
        "",
        html_escape(
          table_data[[column_name]]
        )
      )
    }
  }

  header_cells <- paste0(
    "<th>",
    html_escape(names(table_data)),
    "</th>",
    collapse = ""
  )

  header <- paste0(
    "<tr>",
    header_cells,
    "</tr>"
  )

  rows <- vapply(
    seq_len(nrow(table_data)),
    function(row_index) {
      values <- unlist(
        table_data[
          row_index,
          ,
          drop = FALSE
        ]
      )

      paste0(
        "<tr>",
        paste0(
          "<td>",
          values,
          "</td>",
          collapse = ""
        ),
        "</tr>"
      )
    },
    FUN.VALUE = character(1)
  )

  paste0(
    "<table class='data'>",
    "<thead>",
    header,
    "</thead>",
    "<tbody>",
    paste(rows, collapse = "\n"),
    "</tbody>",
    "</table>"
  )
}

# ============================================================
# 4. Read input files
# ============================================================

sider_drugs <- tryCatch(
  readr::read_csv(
    input_sider,
    show_col_types = FALSE,
    progress = FALSE
  ),
  error = function(error_condition) {
    stop(
      "Failed to read the SIDER crosswalk: ",
      conditionMessage(error_condition),
      call. = FALSE
    )
  }
)

side_effects <- tryCatch(
  readr::read_csv(
    input_side_effects,
    show_col_types = FALSE,
    progress = FALSE
  ),
  error = function(error_condition) {
    stop(
      "Failed to read the SIDER side-effect table: ",
      conditionMessage(error_condition),
      call. = FALSE
    )
  }
)

query_compounds <- tryCatch(
  readr::read_csv(
    input_queries,
    show_col_types = FALSE,
    progress = FALSE
  ),
  error = function(error_condition) {
    stop(
      "Failed to read the A16 compound table: ",
      conditionMessage(error_condition),
      call. = FALSE
    )
  }
)

# ============================================================
# 5. Validate input columns
# ============================================================

required_sider_columns <- c(
  "stitch_flat_id",
  "drug_name",
  "atc_codes",
  "pubchem_cid",
  "canonical_smiles",
  "isomeric_smiles",
  "inchi_key",
  "structure_available"
)

required_side_effect_columns <- c(
  "stitch_flat_id",
  "meddra_umls_id"
)

required_query_columns <- c(
  "compound_id",
  "compound_name",
  "compound_class",
  "analysis_smiles",
  "inchi_key"
)

missing_sider_columns <- setdiff(
  required_sider_columns,
  names(sider_drugs)
)

missing_side_effect_columns <- setdiff(
  required_side_effect_columns,
  names(side_effects)
)

missing_query_columns <- setdiff(
  required_query_columns,
  names(query_compounds)
)

if (length(missing_sider_columns) > 0L) {
  stop(
    "SIDER crosswalk is missing column(s): ",
    paste(
      missing_sider_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(missing_side_effect_columns) > 0L) {
  stop(
    "SIDER side-effect table is missing column(s): ",
    paste(
      missing_side_effect_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(missing_query_columns) > 0L) {
  stop(
    "A16 compound table is missing column(s): ",
    paste(
      missing_query_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

# ============================================================
# 6. Validate input dimensions
# ============================================================

if (nrow(sider_drugs) != 1430L) {
  stop(
    "Expected 1,430 SIDER drugs, found ",
    nrow(sider_drugs),
    ".",
    call. = FALSE
  )
}

if (nrow(query_compounds) != 10L) {
  stop(
    "Expected 10 A16 query compounds, found ",
    nrow(query_compounds),
    ".",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    sider_drugs$stitch_flat_id
  ) > 0L
) {
  stop(
    "Duplicate SIDER flat identifiers detected.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    query_compounds$compound_id
  ) > 0L
) {
  stop(
    "Duplicate A16 query identifiers detected.",
    call. = FALSE
  )
}

# ============================================================
# 7. Prepare structures and side-effect counts
# ============================================================

sider_drugs <- sider_drugs |>
  dplyr::mutate(
    analysis_smiles = dplyr::coalesce(
      clean_text(.data$isomeric_smiles),
      clean_text(.data$canonical_smiles)
    )
  )

query_compounds <- query_compounds |>
  dplyr::mutate(
    analysis_smiles = clean_text(
      .data$analysis_smiles
    )
  )

if (
  any(
    is.na(
      sider_drugs$analysis_smiles
    )
  )
) {
  stop(
    "One or more SIDER drugs lack analysis SMILES.",
    call. = FALSE
  )
}

if (
  any(
    is.na(
      query_compounds$analysis_smiles
    )
  )
) {
  stop(
    "One or more A16 compounds lack analysis SMILES.",
    call. = FALSE
  )
}

side_effect_counts <- side_effects |>
  dplyr::distinct(
    .data$stitch_flat_id,
    .data$meddra_umls_id
  ) |>
  dplyr::count(
    .data$stitch_flat_id,
    name = "unique_side_effect_count"
  )

sider_drugs <- sider_drugs |>
  dplyr::left_join(
    side_effect_counts,
    by = "stitch_flat_id"
  ) |>
  dplyr::mutate(
    unique_side_effect_count = dplyr::coalesce(
      .data$unique_side_effect_count,
      0L
    )
  )

# ============================================================
# 8. Connect R to RDKit
# ============================================================

active_python <- Sys.getenv(
  "RETICULATE_PYTHON"
)

if (active_python == "") {
  active_python <- Sys.which(
    "python"
  )
}

if (active_python == "") {
  stop(
    "Python was not found. Activate the resko-a16 environment.",
    call. = FALSE
  )
}

Sys.setenv(
  RETICULATE_PYTHON = active_python
)

if (
  !reticulate::py_module_available(
    "rdkit"
  )
) {
  stop(
    "RDKit is not available through reticulate using: ",
    active_python,
    call. = FALSE
  )
}

Chem <- reticulate::import(
  "rdkit.Chem",
  convert = FALSE
)

DataStructs <- reticulate::import(
  "rdkit.DataStructs",
  convert = FALSE
)

FingerprintGenerator <- reticulate::import(
  "rdkit.Chem.rdFingerprintGenerator",
  convert = FALSE
)

# ============================================================
# 9. Parse molecular structures
# ============================================================

cat(
  "Parsing 1,430 SIDER structures with RDKit...\n"
)

sider_molecules <- lapply(
  sider_drugs$analysis_smiles,
  function(smiles) {
    Chem$MolFromSmiles(smiles)
  }
)

query_molecules <- lapply(
  query_compounds$analysis_smiles,
  function(smiles) {
    Chem$MolFromSmiles(smiles)
  }
)

sider_valid <- !vapply(
  sider_molecules,
  reticulate::py_is_null_xptr,
  FUN.VALUE = logical(1)
)

query_valid <- !vapply(
  query_molecules,
  reticulate::py_is_null_xptr,
  FUN.VALUE = logical(1)
)

if (any(!sider_valid)) {
  stop(
    "RDKit could not parse ",
    sum(!sider_valid),
    " SIDER structure(s).",
    call. = FALSE
  )
}

if (any(!query_valid)) {
  failed_queries <- query_compounds$compound_name[
    !query_valid
  ]

  stop(
    "RDKit could not parse query structure(s): ",
    paste(
      failed_queries,
      collapse = "; "
    ),
    call. = FALSE
  )
}

# ============================================================
# 10. Generate Morgan fingerprints
# ============================================================

cat(
  "Generating Morgan fingerprints...\n"
)

fingerprint_generator <- FingerprintGenerator$GetMorganGenerator(
  radius = 2L,
  fpSize = 2048L,
  includeChirality = TRUE
)

sider_fingerprints <- lapply(
  sider_molecules,
  function(molecule) {
    fingerprint_generator$GetFingerprint(
      molecule
    )
  }
)

query_fingerprints <- lapply(
  query_molecules,
  function(molecule) {
    fingerprint_generator$GetFingerprint(
      molecule
    )
  }
)

# ============================================================
# 11. Calculate Tanimoto similarities
# ============================================================

cat(
  "Calculating 14,300 SIDER-query comparisons...\n"
)

similarity_matrix <- matrix(
  NA_real_,
  nrow = nrow(sider_drugs),
  ncol = nrow(query_compounds)
)

for (
  query_index in seq_len(
    nrow(query_compounds)
  )
) {
  cat(
    "  Query ",
    query_index,
    " of ",
    nrow(query_compounds),
    ": ",
    query_compounds$compound_name[query_index],
    "\n",
    sep = ""
  )

  for (
    sider_index in seq_len(
      nrow(sider_drugs)
    )
  ) {
    similarity_value <- DataStructs$TanimotoSimilarity(
      sider_fingerprints[[sider_index]],
      query_fingerprints[[query_index]]
    )

    similarity_matrix[
      sider_index,
      query_index
    ] <- as.numeric(
      reticulate::py_to_r(
        similarity_value
      )
    )
  }
}

if (
  any(
    is.na(
      similarity_matrix
    )
  )
) {
  stop(
    "The similarity matrix contains missing values.",
    call. = FALSE
  )
}

if (
  any(
    similarity_matrix < 0 |
      similarity_matrix > 1
  )
) {
  stop(
    "Similarity values outside 0 to 1 were detected.",
    call. = FALSE
  )
}

# ============================================================
# 12. Build all 14,300 comparisons
# ============================================================

all_similarity_rows <- vector(
  "list",
  nrow(query_compounds)
)

for (
  query_index in seq_len(
    nrow(query_compounds)
  )
) {
  query_table <- tibble::tibble(
    query_id = query_compounds$compound_id[
      query_index
    ],
    query_name = query_compounds$compound_name[
      query_index
    ],
    query_class = query_compounds$compound_class[
      query_index
    ],
    sider_stitch_flat_id = sider_drugs$stitch_flat_id,
    sider_drug_name = sider_drugs$drug_name,
    sider_pubchem_cid = sider_drugs$pubchem_cid,
    sider_inchi_key = sider_drugs$inchi_key,
    atc_codes = sider_drugs$atc_codes,
    unique_side_effect_count = sider_drugs$unique_side_effect_count,
    tanimoto_similarity = similarity_matrix[
      ,
      query_index
    ]
  ) |>
    dplyr::mutate(
      tanimoto_distance = 1 - .data$tanimoto_similarity,
      similarity_band = assign_similarity_band(
        .data$tanimoto_similarity
      )
    ) |>
    dplyr::arrange(
      dplyr::desc(
        .data$tanimoto_similarity
      ),
      .data$sider_drug_name
    ) |>
    dplyr::mutate(
      rank_for_query = dplyr::row_number()
    )

  all_similarity_rows[[
    query_index
  ]] <- query_table
}

all_similarity <- dplyr::bind_rows(
  all_similarity_rows
) |>
  dplyr::arrange(
    .data$query_name,
    .data$rank_for_query
  )

if (nrow(all_similarity) != 14300L) {
  stop(
    "Expected 14,300 similarity rows, found ",
    nrow(all_similarity),
    ".",
    call. = FALSE
  )
}

# ============================================================
# 13. Find nearest query for each SIDER drug
# ============================================================

sider_nearest_rows <- vector(
  "list",
  nrow(sider_drugs)
)

for (
  sider_index in seq_len(
    nrow(sider_drugs)
  )
) {
  similarities <- similarity_matrix[
    sider_index,
  ]

  best_query_index <- order(
    -similarities,
    query_compounds$compound_name
  )[1]

  best_similarity <- similarities[
    best_query_index
  ]

  sider_nearest_rows[[
    sider_index
  ]] <- tibble::tibble(
    sider_stitch_flat_id = sider_drugs$stitch_flat_id[
      sider_index
    ],
    sider_drug_name = sider_drugs$drug_name[
      sider_index
    ],
    sider_pubchem_cid = sider_drugs$pubchem_cid[
      sider_index
    ],
    sider_inchi_key = sider_drugs$inchi_key[
      sider_index
    ],
    atc_codes = sider_drugs$atc_codes[
      sider_index
    ],
    unique_side_effect_count = sider_drugs$unique_side_effect_count[
      sider_index
    ],
    nearest_query_id = query_compounds$compound_id[
      best_query_index
    ],
    nearest_query_name = query_compounds$compound_name[
      best_query_index
    ],
    nearest_query_class = query_compounds$compound_class[
      best_query_index
    ],
    nearest_query_tanimoto = best_similarity,
    similarity_band = assign_similarity_band(
      best_similarity
    )
  )
}

sider_nearest <- dplyr::bind_rows(
  sider_nearest_rows
) |>
  dplyr::arrange(
    dplyr::desc(
      .data$nearest_query_tanimoto
    ),
    .data$sider_drug_name
  )

if (nrow(sider_nearest) != 1430L) {
  stop(
    "Expected 1,430 nearest-query rows, found ",
    nrow(sider_nearest),
    ".",
    call. = FALSE
  )
}

# ============================================================
# 14. Select the top 20 neighbours for each query
# ============================================================

query_nearest <- all_similarity |>
  dplyr::filter(
    .data$rank_for_query <= 20L
  ) |>
  dplyr::arrange(
    .data$query_name,
    .data$rank_for_query
  )

if (nrow(query_nearest) != 200L) {
  stop(
    "Expected 200 top-neighbour rows, found ",
    nrow(query_nearest),
    ".",
    call. = FALSE
  )
}

# ============================================================
# 15. Select SIDER similarity candidates
# ============================================================

candidate_threshold <- 0.30

sider_candidates <- sider_nearest |>
  dplyr::filter(
    .data$nearest_query_tanimoto >= candidate_threshold
  ) |>
  dplyr::mutate(
    candidate_scope = "chemical_similarity_only",
    biological_evidence_status = "not_evaluated_in_A17C",
    interpretation = paste(
      "Structurally similar under the selected Morgan fingerprint;",
      "eEF1A binding or inhibition is not established"
    )
  )

# ============================================================
# 16. Summarise each query
# ============================================================

query_summary_rows <- vector(
  "list",
  nrow(query_compounds)
)

for (
  query_index in seq_len(
    nrow(query_compounds)
  )
) {
  query_results <- all_similarity |>
    dplyr::filter(
      .data$query_id ==
        query_compounds$compound_id[
          query_index
        ]
    )

  nearest_row <- query_results |>
    dplyr::slice_min(
      order_by = .data$rank_for_query,
      n = 1L,
      with_ties = FALSE
    )

  query_summary_rows[[
    query_index
  ]] <- tibble::tibble(
    query_id = query_compounds$compound_id[
      query_index
    ],
    query_name = query_compounds$compound_name[
      query_index
    ],
    query_class = query_compounds$compound_class[
      query_index
    ],
    nearest_sider_drug = nearest_row$sider_drug_name[
      1
    ],
    nearest_sider_stitch_id = nearest_row$sider_stitch_flat_id[
      1
    ],
    nearest_tanimoto = nearest_row$tanimoto_similarity[
      1
    ],
    sider_drugs_tanimoto_ge_0_70 = sum(
      query_results$tanimoto_similarity >= 0.70
    ),
    sider_drugs_tanimoto_ge_0_50 = sum(
      query_results$tanimoto_similarity >= 0.50
    ),
    sider_drugs_tanimoto_ge_0_30 = sum(
      query_results$tanimoto_similarity >= 0.30
    ),
    median_tanimoto = stats::median(
      query_results$tanimoto_similarity
    )
  )
}

query_summary <- dplyr::bind_rows(
  query_summary_rows
) |>
  dplyr::arrange(
    dplyr::desc(
      .data$nearest_tanimoto
    ),
    .data$query_name
  )

if (nrow(query_summary) != 10L) {
  stop(
    "Expected 10 query-summary rows.",
    call. = FALSE
  )
}

# ============================================================
# 17. Create similarity plot
# ============================================================

plot_data <- query_nearest |>
  dplyr::mutate(
    plot_label = paste0(
      .data$rank_for_query,
      ". ",
      .data$sider_drug_name
    )
  )

similarity_plot <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$tanimoto_similarity,
    y = stats::reorder(
      .data$plot_label,
      .data$tanimoto_similarity
    )
  )
) +
  ggplot2::geom_col(
    fill = "#E67E22"
  ) +
  ggplot2::facet_wrap(
    ggplot2::vars(
      query_name
    ),
    scales = "free_y",
    ncol = 2
  ) +
  ggplot2::coord_cartesian(
    xlim = c(
      0,
      1
    )
  ) +
  ggplot2::labs(
    title = paste(
      "Top 20 SIDER structural neighbours",
      "for each RESKO query"
    ),
    x = "Morgan fingerprint Tanimoto similarity",
    y = NULL
  ) +
  ggplot2::theme_minimal(
    base_size = 10
  ) +
  ggplot2::theme(
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.y = ggplot2::element_text(
      size = 6
    ),
    plot.title = ggplot2::element_text(
      face = "bold"
    )
  )

safe_save_plot(
  plot = similarity_plot,
  path = output_plot,
  width = 13,
  height = 18
)

# ============================================================
# 18. Create summary table
# ============================================================

summary_table <- tibble::tibble(
  metric = c(
    "sider_drug_count",
    "query_compound_count",
    "total_similarity_comparisons",
    "candidate_threshold",
    "sider_candidates_tanimoto_ge_0_30",
    "sider_candidates_tanimoto_ge_0_50",
    "sider_candidates_tanimoto_ge_0_70",
    "maximum_observed_tanimoto",
    "queries_with_close_neighbour_ge_0_70",
    "queries_with_moderate_neighbour_ge_0_50"
  ),
  value = c(
    as.character(
      nrow(sider_drugs)
    ),
    as.character(
      nrow(query_compounds)
    ),
    as.character(
      length(similarity_matrix)
    ),
    as.character(
      candidate_threshold
    ),
    as.character(
      sum(
        sider_nearest$nearest_query_tanimoto >= 0.30
      )
    ),
    as.character(
      sum(
        sider_nearest$nearest_query_tanimoto >= 0.50
      )
    ),
    as.character(
      sum(
        sider_nearest$nearest_query_tanimoto >= 0.70
      )
    ),
    as.character(
      max(similarity_matrix)
    ),
    as.character(
      sum(
        query_summary$nearest_tanimoto >= 0.70
      )
    ),
    as.character(
      sum(
        query_summary$nearest_tanimoto >= 0.50
      )
    )
  )
)

# ============================================================
# 19. Write CSV outputs
# ============================================================

safe_write_csv(
  all_similarity,
  output_all_similarity
)

safe_write_csv(
  sider_nearest,
  output_sider_nearest
)

safe_write_csv(
  query_nearest,
  output_query_nearest
)

safe_write_csv(
  sider_candidates,
  output_candidates
)

safe_write_csv(
  query_summary,
  output_query_summary
)

safe_write_csv(
  summary_table,
  output_summary
)

# ============================================================
# 20. Create standalone HTML report
# ============================================================

plot_uri <- base64enc::dataURI(
  file = output_plot,
  mime = "image/png"
)

candidate_report_table <- sider_candidates |>
  dplyr::select(
    "sider_drug_name",
    "nearest_query_name",
    "nearest_query_tanimoto",
    "similarity_band",
    "unique_side_effect_count",
    "atc_codes"
  )

query_report_table <- query_summary |>
  dplyr::select(
    "query_name",
    "query_class",
    "nearest_sider_drug",
    "nearest_tanimoto",
    "sider_drugs_tanimoto_ge_0_70",
    "sider_drugs_tanimoto_ge_0_50",
    "sider_drugs_tanimoto_ge_0_30"
  )

report_html <- paste0(
  "<!doctype html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>RESKO A17C SIDER Similarity Report</title>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;",
  "max-width:1250px;margin:0 auto;padding:28px;color:#202124;line-height:1.5}",
  "h1,h2{color:#263238}",
  ".note{background:#fff4e5;border-left:5px solid #e67e22;",
  "padding:14px;margin:16px 0}",
  "table.data{border-collapse:collapse;width:100%;font-size:.88rem}",
  "table.data th,table.data td{border:1px solid #ddd;padding:7px;text-align:left}",
  "table.data th{background:#f3f4f6}",
  "img{max-width:100%;height:auto;border:1px solid #ddd}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>RESKO A17C SIDER Chemical Similarity Report</h1>",
  "<p>All 1,430 SIDER drugs were compared with six RESKO network candidates ",
  "and four established eEF1A reference ligands.</p>",
  "<div class='note'>",
  "<strong>Interpretive boundary:</strong> ",
  "Chemical similarity identifies structural neighbours only. ",
  "It does not establish eEF1A binding, inhibition, translation effects, ",
  "biological efficacy, selectivity, or safety.",
  "</div>",
  "<h2>Query summary</h2>",
  table_to_html(
    query_report_table
  ),
  "<h2>SIDER similarity candidates</h2>",
  "<p>The compounds below had a Tanimoto similarity of at least 0.30 ",
  "to their nearest A16 query compound.</p>",
  table_to_html(
    candidate_report_table
  ),
  "<h2>Top 20 SIDER neighbours per query</h2>",
  "",
  plot_uri,
  "",
  "<h2>Method</h2>",
  "<p>RDKit Morgan fingerprints used radius 2, 2,048 bits, and chirality. ",
  "Similarity was calculated with the Tanimoto coefficient. ",
  "Side-effect counts were derived from unique SIDER flat-drug and ",
  "MedDRA preferred-term pairs.</p>",
  "<h2>Limitations</h2>",
  "<p>SIDER 4.1 is a historical marketed-drug resource. ",
  "A17C assessed chemical similarity only. ",
  "The similarity thresholds are operational prioritisation categories. ",
  "Biological relevance must be assessed using target, assay, pathway, ",
  "network, and provenance evidence in later stages.</p>",
  "</body>",
  "</html>"
)

temporary_report <- tempfile(
  pattern = "A17C_report_",
  tmpdir = results_dir,
  fileext = ".html"
)

writeLines(
  report_html,
  temporary_report,
  useBytes = TRUE
)

if (
  !file.exists(temporary_report) ||
  is.na(file.info(temporary_report)$size) ||
  file.info(temporary_report)$size <= 0L
) {
  stop(
    "The A17C HTML report was not created.",
    call. = FALSE
  )
}

backup_existing_file(
  output_report
)

report_moved <- file.rename(
  from = temporary_report,
  to = output_report
)

if (!report_moved) {
  stop(
    "Could not move the HTML report into place.",
    call. = FALSE
  )
}

# ============================================================
# 21. Verify all outputs
# ============================================================

required_output_files <- c(
  output_all_similarity,
  output_sider_nearest,
  output_query_nearest,
  output_candidates,
  output_query_summary,
  output_plot,
  output_report,
  output_summary
)

for (path in required_output_files) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Output verification failed: ",
      path,
      call. = FALSE
    )
  }
}

all_similarity_check <- readr::read_csv(
  output_all_similarity,
  show_col_types = FALSE,
  progress = FALSE
)

sider_nearest_check <- readr::read_csv(
  output_sider_nearest,
  show_col_types = FALSE,
  progress = FALSE
)

query_nearest_check <- readr::read_csv(
  output_query_nearest,
  show_col_types = FALSE,
  progress = FALSE
)

candidate_check <- readr::read_csv(
  output_candidates,
  show_col_types = FALSE,
  progress = FALSE
)

query_summary_check <- readr::read_csv(
  output_query_summary,
  show_col_types = FALSE,
  progress = FALSE
)

summary_check <- readr::read_csv(
  output_summary,
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(all_similarity_check) != 14300L) {
  stop(
    "The all-similarity output does not contain 14,300 rows.",
    call. = FALSE
  )
}

if (nrow(sider_nearest_check) != 1430L) {
  stop(
    "The SIDER nearest-query output does not contain 1,430 rows.",
    call. = FALSE
  )
}

if (nrow(query_nearest_check) != 200L) {
  stop(
    "The query nearest-neighbour output does not contain 200 rows.",
    call. = FALSE
  )
}

if (nrow(candidate_check) != nrow(sider_candidates)) {
  stop(
    "The candidate output row count is incorrect.",
    call. = FALSE
  )
}

if (nrow(query_summary_check) != 10L) {
  stop(
    "The query summary does not contain 10 rows.",
    call. = FALSE
  )
}

if (nrow(summary_check) != 10L) {
  stop(
    "The A17C summary does not contain 10 rows.",
    call. = FALSE
  )
}

# ============================================================
# 22. Completion summary
# ============================================================

cat(
  "A17C SIDER chemical-similarity screening completed.\n"
)

cat(
  "RDKit Python: ",
  active_python,
  "\n",
  sep = ""
)

cat(
  "SIDER drugs screened: ",
  nrow(sider_drugs),
  "\n",
  sep = ""
)

cat(
  "Query compounds: ",
  nrow(query_compounds),
  "\n",
  sep = ""
)

cat(
  "Similarity comparisons: ",
  length(similarity_matrix),
  "\n",
  sep = ""
)

cat(
  "SIDER candidates with Tanimoto at least 0.30: ",
  nrow(sider_candidates),
  "\n",
  sep = ""
)

cat(
  "SIDER candidates with Tanimoto at least 0.50: ",
  sum(
    sider_nearest$nearest_query_tanimoto >= 0.50
  ),
  "\n",
  sep = ""
)

cat(
  "SIDER candidates with Tanimoto at least 0.70: ",
  sum(
    sider_nearest$nearest_query_tanimoto >= 0.70
  ),
  "\n",
  sep = ""
)

cat(
  "Maximum observed Tanimoto similarity: ",
  signif(
    max(similarity_matrix),
    4
  ),
  "\n",
  sep = ""
)

cat(
  "Verified outputs: ",
  length(required_output_files),
  "\n",
  sep = ""
)