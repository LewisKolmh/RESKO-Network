#!/usr/bin/env Rscript

# ============================================================
# RESKO A17D: Combined SIDER candidate evaluation
# ============================================================
#
# Run from the RESKO project root:
#
# Rscript scripts/A17D_evaluate_sider_candidates.R
#
# Inputs:
#
# results/A17C_sider_similarity_candidates.csv
# results/A17A_sider_side_effects_pt.csv
# results/A17A_sider_frequencies.csv
# results/A17A_sider_indications_pt.csv
# results/A8_target_proteins.csv
#
# Alternative network input:
#
# results/nodes_proteins.csv
#
# Outputs:
#
# results/A17D_candidate_annotations.csv
# results/A17D_candidate_activities.csv
# results/A17D_candidate_target_evidence.csv
# results/A17D_candidate_side_effects.csv
# results/A17D_candidate_frequencies.csv
# results/A17D_candidate_indications.csv
# results/A17D_side_effect_similarity.csv
# results/A17D_candidate_classification.csv
# results/A17D_summary.csv
# results/A17D_candidate_report.html
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# ============================================================
# 1. Validate R packages
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "tibble",
  "httr2",
  "jsonlite"
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
# 2. Define project paths
# ============================================================

project_root <- normalizePath(
  getwd(),
  winslash = "/",
  mustWork = TRUE
)

results_dir <- file.path(
  project_root,
  "results"
)

scripts_dir <- file.path(
  project_root,
  "scripts"
)

cache_dir <- file.path(
  results_dir,
  "A17D_api_cache"
)

if (!dir.exists(results_dir)) {
  stop(
    "The results directory was not found. ",
    "Run A17D from the RESKO project root.",
    call. = FALSE
  )
}

if (!dir.exists(scripts_dir)) {
  stop(
    "The scripts directory was not found. ",
    "Run A17D from the RESKO project root.",
    call. = FALSE
  )
}

dir.create(
  cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_candidates <- file.path(
  results_dir,
  "A17C_sider_similarity_candidates.csv"
)

input_side_effects <- file.path(
  results_dir,
  "A17A_sider_side_effects_pt.csv"
)

input_frequencies <- file.path(
  results_dir,
  "A17A_sider_frequencies.csv"
)

input_indications <- file.path(
  results_dir,
  "A17A_sider_indications_pt.csv"
)

network_options <- c(
  file.path(
    results_dir,
    "A8_target_proteins.csv"
  ),
  file.path(
    results_dir,
    "nodes_proteins.csv"
  )
)

available_network_files <- network_options[
  file.exists(network_options)
]

if (length(available_network_files) == 0L) {
  stop(
    paste(
      "Neither results/A8_target_proteins.csv",
      "nor results/nodes_proteins.csv was found."
    ),
    call. = FALSE
  )
}

input_network <- available_network_files[1]

outputs <- c(
  annotations = file.path(
    results_dir,
    "A17D_candidate_annotations.csv"
  ),
  activities = file.path(
    results_dir,
    "A17D_candidate_activities.csv"
  ),
  target_evidence = file.path(
    results_dir,
    "A17D_candidate_target_evidence.csv"
  ),
  side_effects = file.path(
    results_dir,
    "A17D_candidate_side_effects.csv"
  ),
  frequencies = file.path(
    results_dir,
    "A17D_candidate_frequencies.csv"
  ),
  indications = file.path(
    results_dir,
    "A17D_candidate_indications.csv"
  ),
  side_effect_similarity = file.path(
    results_dir,
    "A17D_side_effect_similarity.csv"
  ),
  classification = file.path(
    results_dir,
    "A17D_candidate_classification.csv"
  ),
  summary = file.path(
    results_dir,
    "A17D_summary.csv"
  ),
  report = file.path(
    results_dir,
    "A17D_candidate_report.html"
  )
)

required_inputs <- c(
  input_candidates,
  input_side_effects,
  input_frequencies,
  input_indications,
  input_network
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

normalise_name <- function(x) {
  output <- tolower(
    clean_text(x)
  )

  output <- gsub(
    "[^a-z0-9]+",
    "",
    output
  )

  output[
    output == ""
  ] <- NA_character_

  output
}

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

collapse_unique <- function(x) {
  values <- clean_text(x)
  values <- values[!is.na(values)]
  values <- sort(unique(values))

  if (length(values) == 0L) {
    return(NA_character_)
  }

  paste(
    values,
    collapse = "; "
  )
}

get_column <- function(
  data,
  column_name,
  default = NA
) {
  if (column_name %in% names(data)) {
    return(data[[column_name]])
  }

  rep(
    default,
    nrow(data)
  )
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

safe_write_csv <- function(
  data,
  path
) {
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
      "Failed to create output: ",
      path,
      call. = FALSE
    )
  }

  backup_existing_file(path)

  moved <- file.rename(
    temporary_path,
    path
  )

  if (!moved) {
    stop(
      "Could not move output into place: ",
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

request_json <- function(
  url,
  cache_name,
  attempts = 5L
) {
  cache_path <- file.path(
    cache_dir,
    cache_name
  )

  if (
    file.exists(cache_path) &&
    !is.na(file.info(cache_path)$size) &&
    file.info(cache_path)$size > 0L
  ) {
    cached <- tryCatch(
      jsonlite::fromJSON(
        cache_path,
        simplifyDataFrame = TRUE
      ),
      error = function(error) error
    )

    if (!inherits(cached, "error")) {
      return(cached)
    }

    unlink(cache_path)
  }

  last_error <- "No request completed"

  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_user_agent(
          "RESKO-A17D/4.0"
        ) |>
        httr2::req_timeout(
          seconds = 120
        ) |>
        httr2::req_perform(),
      error = function(error) error
    )

    if (!inherits(response, "error")) {
      status <- httr2::resp_status(
        response
      )

      if (
        status >= 200L &&
        status < 300L
      ) {
        body <- httr2::resp_body_string(
          response
        )

        parsed <- tryCatch(
          jsonlite::fromJSON(
            body,
            simplifyDataFrame = TRUE
          ),
          error = function(error) error
        )

        if (!inherits(parsed, "error")) {
          writeLines(
            body,
            cache_path,
            useBytes = TRUE
          )

          Sys.sleep(0.25)

          return(parsed)
        }

        last_error <- conditionMessage(
          parsed
        )
      } else {
        last_error <- paste0(
          "HTTP ",
          status
        )
      }
    } else {
      last_error <- conditionMessage(
        response
      )
    }

    if (attempt < attempts) {
      Sys.sleep(
        min(
          2 ^ (attempt - 1L),
          16
        )
      )
    }
  }

  stop(
    "API request failed for ",
    url,
    ". Last error: ",
    last_error,
    call. = FALSE
  )
}

fetch_activities <- function(
  molecule_chembl_id,
  candidate_name
) {
  limit <- 1000L
  offset <- 0L
  page_number <- 1L
  pages <- list()

  repeat {
    url <- paste0(
      "https://www.ebi.ac.uk/chembl/api/data/activity.json",
      "?molecule_chembl_id=",
      utils::URLencode(
        molecule_chembl_id,
        reserved = TRUE
      ),
      "&limit=",
      limit,
      "&offset=",
      offset
    )

    cache_name <- paste0(
      "activities_",
      molecule_chembl_id,
      "_",
      sprintf("%05d", offset),
      ".json"
    )

    parsed <- request_json(
      url,
      cache_name
    )

    if (
      is.null(parsed$activities) ||
      !is.data.frame(parsed$activities)
    ) {
      break
    }

    page <- tibble::as_tibble(
      parsed$activities
    )

    if (nrow(page) == 0L) {
      break
    }

    page$candidate_name <- candidate_name
    page$queried_molecule_chembl_id <-
      molecule_chembl_id

    pages[[page_number]] <- page

    total_count <- if (
      !is.null(parsed$page_meta$total_count)
    ) {
      as.integer(
        parsed$page_meta$total_count
      )
    } else {
      NA_integer_
    }

    offset <- offset + nrow(page)
    page_number <- page_number + 1L

    if (
      !is.na(total_count) &&
      offset >= total_count
    ) {
      break
    }

    if (nrow(page) < limit) {
      break
    }
  }

  dplyr::bind_rows(
    pages
  )
}

make_html_table <- function(data) {
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
            3
          ),
          nsmall = 3,
          trim = TRUE
        )
      )
    } else {
      table_data[[column_name]] <- ifelse(
        is.na(table_data[[column_name]]),
        "",
        as.character(
          table_data[[column_name]]
        )
      )
    }
  }

  header <- paste0(
    "<tr>",
    paste0(
      "<th>",
      names(table_data),
      "</th>",
      collapse = ""
    ),
    "</tr>"
  )

  table_rows <- vapply(
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
    "<table>",
    "<thead>",
    header,
    "</thead>",
    "<tbody>",
    paste(
      table_rows,
      collapse = "\n"
    ),
    "</tbody>",
    "</table>"
  )
}

# ============================================================
# 4. Read and validate inputs
# ============================================================

cat(
  "Reading A17D inputs...\n"
)

candidates_raw <- readr::read_csv(
  input_candidates,
  show_col_types = FALSE,
  progress = FALSE
)

side_effects_all <- readr::read_csv(
  input_side_effects,
  show_col_types = FALSE,
  progress = FALSE
)

frequencies_all <- readr::read_csv(
  input_frequencies,
  show_col_types = FALSE,
  progress = FALSE
)

indications_all <- readr::read_csv(
  input_indications,
  show_col_types = FALSE,
  progress = FALSE
)

network_raw <- readr::read_csv(
  input_network,
  show_col_types = FALSE,
  progress = FALSE
)

required_candidate_columns <- c(
  "sider_stitch_flat_id",
  "sider_drug_name",
  "sider_pubchem_cid",
  "sider_inchi_key",
  "nearest_query_name",
  "nearest_query_tanimoto",
  "similarity_band",
  "unique_side_effect_count",
  "atc_codes"
)

missing_candidate_columns <- setdiff(
  required_candidate_columns,
  names(candidates_raw)
)

if (length(missing_candidate_columns) > 0L) {
  stop(
    "Candidate input is missing column(s): ",
    paste(
      missing_candidate_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

expected_candidate_names <- c(
  "nilotinib",
  "imatinib",
  "alprazolam",
  "triazolam",
  "ponatinib",
  "temazepam"
)

candidates <- candidates_raw |>
  dplyr::mutate(
    normalised_name = normalise_name(
      .data$sider_drug_name
    )
  ) |>
  dplyr::filter(
    .data$normalised_name %in%
      expected_candidate_names
  ) |>
  dplyr::distinct(
    .data$sider_stitch_flat_id,
    .keep_all = TRUE
  ) |>
  dplyr::arrange(
    match(
      .data$normalised_name,
      expected_candidate_names
    )
  )

if (nrow(candidates) != 6L) {
  stop(
    "Expected six A17C candidates, found ",
    nrow(candidates),
    ".",
    call. = FALSE
  )
}

missing_expected_candidates <- setdiff(
  expected_candidate_names,
  candidates$normalised_name
)

if (
  length(
    missing_expected_candidates
  ) > 0L
) {
  stop(
    "Missing expected candidate(s): ",
    paste(
      missing_expected_candidates,
      collapse = ", "
    ),
    call. = FALSE
  )
}

protein_columns <- intersect(
  c(
    "protein",
    "gene_symbol",
    "symbol",
    "protein_name",
    "display_name"
  ),
  names(network_raw)
)

if (length(protein_columns) == 0L) {
  stop(
    "Could not identify a protein column in ",
    input_network,
    ".",
    call. = FALSE
  )
}

network_proteins <- unique(
  clean_text(
    network_raw[[
      protein_columns[1]
    ]]
  )
)

network_proteins <- network_proteins[
  !is.na(network_proteins)
]

network_proteins <- toupper(
  network_proteins
)

# ============================================================
# 5. Extract local SIDER evidence
# ============================================================

cat(
  "Extracting SIDER evidence...\n"
)

candidate_ids <- candidates$sider_stitch_flat_id

candidate_name_table <- candidates |>
  dplyr::select(
    stitch_flat_id =
      "sider_stitch_flat_id",
    candidate_name =
      "sider_drug_name"
  )

candidate_side_effects <- side_effects_all |>
  dplyr::filter(
    .data$stitch_flat_id %in%
      candidate_ids
  ) |>
  dplyr::left_join(
    candidate_name_table,
    by = "stitch_flat_id"
  ) |>
  dplyr::distinct(
    .data$stitch_flat_id,
    .data$meddra_umls_id,
    .keep_all = TRUE
  ) |>
  dplyr::arrange(
    .data$candidate_name,
    .data$side_effect_name
  )

candidate_frequencies <- frequencies_all |>
  dplyr::filter(
    .data$stitch_flat_id %in%
      candidate_ids
  ) |>
  dplyr::left_join(
    candidate_name_table,
    by = "stitch_flat_id"
  ) |>
  dplyr::arrange(
    .data$candidate_name,
    .data$side_effect_name
  )

candidate_indications <- indications_all |>
  dplyr::filter(
    .data$stitch_flat_id %in%
      candidate_ids
  ) |>
  dplyr::left_join(
    candidate_name_table,
    by = "stitch_flat_id"
  ) |>
  dplyr::distinct(
    .data$stitch_flat_id,
    .data$meddra_umls_id,
    .keep_all = TRUE
  ) |>
  dplyr::arrange(
    .data$candidate_name,
    .data$indication_name
  )

# ============================================================
# 6. Map the six candidates to ChEMBL
# ============================================================

cat(
  "Mapping candidates to ChEMBL...\n"
)

annotation_rows <- vector(
  "list",
  nrow(candidates)
)

for (
  candidate_index in seq_len(
    nrow(candidates)
  )
) {
  candidate_name <-
    candidates$sider_drug_name[
      candidate_index
    ]

  inchi_key <-
    candidates$sider_inchi_key[
      candidate_index
    ]

  encoded_inchi_key <- utils::URLencode(
    inchi_key,
    reserved = TRUE
  )

  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule.json",
    "?molecule_structures__standard_inchi_key=",
    encoded_inchi_key
  )

  cache_name <- paste0(
    "molecule_",
    candidates$normalised_name[
      candidate_index
    ],
    ".json"
  )

  parsed <- request_json(
    url,
    cache_name
  )

  molecules <- if (
    !is.null(parsed$molecules) &&
    is.data.frame(parsed$molecules)
  ) {
    tibble::as_tibble(
      parsed$molecules
    )
  } else {
    tibble::tibble()
  }

  if (nrow(molecules) == 0L) {
    annotation_rows[[
      candidate_index
    ]] <- tibble::tibble(
      candidate_name = candidate_name,
      sider_stitch_flat_id =
        candidates$sider_stitch_flat_id[
          candidate_index
        ],
      pubchem_cid =
        candidates$sider_pubchem_cid[
          candidate_index
        ],
      inchi_key = inchi_key,
      molecule_chembl_id = NA_character_,
      pref_name = NA_character_,
      molecule_type = NA_character_,
      max_phase = NA_real_,
      first_approval = NA_real_,
      mapping_status = "not_mapped"
    )

    next
  }

  molecule <- molecules[
    1,
    ,
    drop = FALSE
  ]

  molecule_chembl_id <- clean_text(
    get_column(
      molecule,
      "molecule_chembl_id"
    )
  )[1]

  pref_name <- clean_text(
    get_column(
      molecule,
      "pref_name"
    )
  )[1]

  molecule_type <- clean_text(
    get_column(
      molecule,
      "molecule_type"
    )
  )[1]

  max_phase <- safe_numeric(
    get_column(
      molecule,
      "max_phase"
    )
  )[1]

  first_approval <- safe_numeric(
    get_column(
      molecule,
      "first_approval"
    )
  )[1]

  mapping_status <- if (
    nrow(molecules) == 1L
  ) {
    "mapped_unique"
  } else {
    "multiple_records_first_retained"
  }

  annotation_rows[[
    candidate_index
  ]] <- tibble::tibble(
    candidate_name = candidate_name,
    sider_stitch_flat_id =
      candidates$sider_stitch_flat_id[
        candidate_index
      ],
    pubchem_cid =
      candidates$sider_pubchem_cid[
        candidate_index
      ],
    inchi_key = inchi_key,
    molecule_chembl_id =
      molecule_chembl_id,
    pref_name = pref_name,
    molecule_type =
      molecule_type,
    max_phase = max_phase,
    first_approval =
      first_approval,
    mapping_status =
      mapping_status
  )
}

candidate_annotations <- dplyr::bind_rows(
  annotation_rows
) |>
  dplyr::left_join(
    candidates |>
      dplyr::select(
        candidate_name =
          "sider_drug_name",
        "nearest_query_name",
        "nearest_query_tanimoto",
        "similarity_band",
        "unique_side_effect_count",
        "atc_codes"
      ),
    by = "candidate_name"
  )

mapped_candidate_count <- sum(
  !is.na(
    candidate_annotations$molecule_chembl_id
  )
)

if (mapped_candidate_count == 0L) {
  stop(
    "None of the six candidates mapped to ChEMBL.",
    call. = FALSE
  )
}

# ============================================================
# 7. Retrieve paginated ChEMBL activities
# ============================================================

cat(
  "Retrieving ChEMBL activities...\n"
)

activity_tables <- list()
activity_position <- 1L

for (
  candidate_index in seq_len(
    nrow(candidate_annotations)
  )
) {
  molecule_chembl_id <-
    candidate_annotations$molecule_chembl_id[
      candidate_index
    ]

  candidate_name <-
    candidate_annotations$candidate_name[
      candidate_index
    ]

  if (
    is.na(molecule_chembl_id) ||
    molecule_chembl_id == ""
  ) {
    next
  }

  cat(
    "  ",
    candidate_name,
    " ",
    molecule_chembl_id,
    "\n",
    sep = ""
  )

  activity_table <- fetch_activities(
    molecule_chembl_id,
    candidate_name
  )

  if (nrow(activity_table) > 0L) {
    activity_tables[[
      activity_position
    ]] <- activity_table

    activity_position <-
      activity_position + 1L
  }
}

activities_raw <- dplyr::bind_rows(
  activity_tables
)

required_activity_columns <- c(
  "candidate_name",
  "molecule_chembl_id",
  "queried_molecule_chembl_id",
  "activity_id",
  "target_chembl_id",
  "target_pref_name",
  "target_organism",
  "assay_chembl_id",
  "document_chembl_id",
  "standard_type",
  "standard_relation",
  "standard_value",
  "standard_units",
  "pchembl_value",
  "data_validity_comment",
  "potential_duplicate"
)

for (
  column_name in required_activity_columns
) {
  if (
    !column_name %in%
      names(activities_raw)
  ) {
    activities_raw[[
      column_name
    ]] <- NA
  }
}

if (nrow(activities_raw) == 0L) {
  candidate_activities <- tibble::tibble(
    candidate_name = character(),
    molecule_chembl_id = character(),
    activity_id = character(),
    target_chembl_id = character(),
    target_pref_name = character(),
    target_organism = character(),
    assay_chembl_id = character(),
    document_chembl_id = character(),
    standard_type = character(),
    standard_relation = character(),
    standard_value = numeric(),
    standard_units = character(),
    pchembl_value = numeric(),
    data_validity_comment = character(),
    potential_duplicate = logical()
  )
} else {
  resolved_molecule_id <- dplyr::coalesce(
    clean_text(
      activities_raw$molecule_chembl_id
    ),
    clean_text(
      activities_raw$queried_molecule_chembl_id
    )
  )

  candidate_activities <- tibble::tibble(
    candidate_name = clean_text(
      activities_raw$candidate_name
    ),
    molecule_chembl_id =
      resolved_molecule_id,
    activity_id = clean_text(
      activities_raw$activity_id
    ),
    target_chembl_id = clean_text(
      activities_raw$target_chembl_id
    ),
    target_pref_name = clean_text(
      activities_raw$target_pref_name
    ),
    target_organism = clean_text(
      activities_raw$target_organism
    ),
    assay_chembl_id = clean_text(
      activities_raw$assay_chembl_id
    ),
    document_chembl_id = clean_text(
      activities_raw$document_chembl_id
    ),
    standard_type = clean_text(
      activities_raw$standard_type
    ),
    standard_relation = clean_text(
      activities_raw$standard_relation
    ),
    standard_value = safe_numeric(
      activities_raw$standard_value
    ),
    standard_units = clean_text(
      activities_raw$standard_units
    ),
    pchembl_value = safe_numeric(
      activities_raw$pchembl_value
    ),
    data_validity_comment = clean_text(
      activities_raw$data_validity_comment
    ),
    potential_duplicate = as.logical(
      activities_raw$potential_duplicate
    )
  ) |>
    dplyr::distinct()
}

# ============================================================
# 8. Classify ChEMBL target evidence
# ============================================================

cat(
  "Classifying target evidence...\n"
)

if (nrow(candidate_activities) == 0L) {
  candidate_target_evidence <-
    candidate_activities |>
    dplyr::mutate(
      matched_network_proteins =
        character(),
      target_evidence_class =
        character()
    )
} else {
  target_text <- toupper(
    dplyr::coalesce(
      candidate_activities$target_pref_name,
      ""
    )
  )

  matched_network_proteins <- vapply(
    target_text,
    function(current_target) {
      matches <- network_proteins[
        vapply(
          network_proteins,
          function(protein) {
            pattern <- paste0(
              "(^|[^A-Z0-9])",
              protein,
              "([^A-Z0-9]|$)"
            )

            grepl(
              pattern,
              current_target
            )
          },
          FUN.VALUE = logical(1)
        )
      ]

      if (length(matches) == 0L) {
        NA_character_
      } else {
        paste(
          sort(unique(matches)),
          collapse = "; "
        )
      }
    },
    FUN.VALUE = character(1)
  )

  target_evidence_class <- dplyr::case_when(
    grepl(
      "EEF1A1|EEF1A2",
      target_text
    ) ~ "direct_eef1a",

    grepl(
      "EEF1B2|EEF1D|EEF1G",
      target_text
    ) ~ "eef1_complex",

    !is.na(
      matched_network_proteins
    ) ~ "translation_network",

    TRUE ~ "other_target"
  )

  candidate_target_evidence <-
    candidate_activities |>
    dplyr::mutate(
      matched_network_proteins =
        matched_network_proteins,
      target_evidence_class =
        target_evidence_class
    )
}

# ============================================================
# 9. Calculate pairwise side-effect similarities
# ============================================================

cat(
  "Calculating side-effect similarities...\n"
)

drug_names <- candidates$sider_drug_name

side_effect_sets <- split(
  candidate_side_effects$meddra_umls_id,
  candidate_side_effects$candidate_name
)

for (drug_name in drug_names) {
  if (is.null(side_effect_sets[[drug_name]])) {
    side_effect_sets[[drug_name]] <-
      character()
  }
}

side_effect_sets <- lapply(
  side_effect_sets,
  function(values) {
    values <- clean_text(values)
    values <- values[!is.na(values)]
    sort(unique(values))
  }
)

all_side_effect_terms <- sort(
  unique(
    unlist(
      side_effect_sets,
      use.names = FALSE
    )
  )
)

if (
  length(all_side_effect_terms) == 0L
) {
  stop(
    "No SIDER side effects were found for the six candidates.",
    call. = FALSE
  )
}

presence_matrix <- matrix(
  0,
  nrow = length(drug_names),
  ncol = length(all_side_effect_terms)
)

rownames(presence_matrix) <-
  drug_names

colnames(presence_matrix) <-
  all_side_effect_terms

for (
  drug_index in seq_along(
    drug_names
  )
) {
  drug_name <- drug_names[
    drug_index
  ]

  drug_terms <- side_effect_sets[[
    drug_name
  ]]

  if (length(drug_terms) > 0L) {
    term_positions <- match(
      drug_terms,
      all_side_effect_terms
    )

    presence_matrix[
      drug_index,
      term_positions
    ] <- 1
  }
}

document_frequency <- colSums(
  presence_matrix > 0
)

inverse_document_frequency <- log(
  (
    1 + nrow(presence_matrix)
  ) /
    (
      1 + document_frequency
    )
) + 1

tfidf_matrix <- sweep(
  presence_matrix,
  2,
  inverse_document_frequency,
  "*"
)

pair_indices <- utils::combn(
  seq_along(drug_names),
  2L
)

similarity_rows <- vector(
  "list",
  ncol(pair_indices)
)

kinase_group <- c(
  "nilotinib",
  "imatinib",
  "ponatinib"
)

benzodiazepine_group <- c(
  "alprazolam",
  "triazolam",
  "temazepam"
)

for (
  pair_index in seq_len(
    ncol(pair_indices)
  )
) {
  first_index <- pair_indices[
    1,
    pair_index
  ]

  second_index <- pair_indices[
    2,
    pair_index
  ]

  first_binary <- presence_matrix[
    first_index,
  ]

  second_binary <- presence_matrix[
    second_index,
  ]

  shared_count <- sum(
    first_binary == 1 &
      second_binary == 1
  )

  union_count <- sum(
    first_binary == 1 |
      second_binary == 1
  )

  jaccard_similarity <- if (
    union_count == 0L
  ) {
    NA_real_
  } else {
    shared_count / union_count
  }

  first_tfidf <- tfidf_matrix[
    first_index,
  ]

  second_tfidf <- tfidf_matrix[
    second_index,
  ]

  cosine_denominator <- sqrt(
    sum(first_tfidf ^ 2)
  ) * sqrt(
    sum(second_tfidf ^ 2)
  )

  tfidf_cosine_similarity <- if (
    cosine_denominator == 0
  ) {
    NA_real_
  } else {
    sum(
      first_tfidf *
        second_tfidf
    ) /
      cosine_denominator
  }

  first_name <- drug_names[
    first_index
  ]

  second_name <- drug_names[
    second_index
  ]

  comparison_group <- if (
    first_name %in% kinase_group &&
    second_name %in% kinase_group
  ) {
    "kinase_inhibitor_group"
  } else if (
    first_name %in%
      benzodiazepine_group &&
    second_name %in%
      benzodiazepine_group
  ) {
    "benzodiazepine_related_group"
  } else {
    "cross_group"
  }

  similarity_rows[[
    pair_index
  ]] <- tibble::tibble(
    drug_1 = first_name,
    drug_2 = second_name,
    comparison_group =
      comparison_group,
    drug_1_side_effect_count =
      sum(first_binary),
    drug_2_side_effect_count =
      sum(second_binary),
    shared_side_effect_count =
      shared_count,
    union_side_effect_count =
      union_count,
    jaccard_similarity =
      jaccard_similarity,
    tfidf_cosine_similarity =
      tfidf_cosine_similarity
  )
}

side_effect_similarity <- dplyr::bind_rows(
  similarity_rows
) |>
  dplyr::arrange(
    dplyr::desc(
      .data$tfidf_cosine_similarity
    ),
    dplyr::desc(
      .data$jaccard_similarity
    )
  )

if (
  nrow(side_effect_similarity) != 15L
) {
  stop(
    "Expected 15 pairwise side-effect comparisons, found ",
    nrow(side_effect_similarity),
    ".",
    call. = FALSE
  )
}

# ============================================================
# 10. Summarise and classify biological evidence
# ============================================================

if (
  nrow(candidate_target_evidence) == 0L
) {
  target_summary <- tibble::tibble(
    candidate_name = character(),
    activity_record_count = integer(),
    unique_target_count = integer(),
    direct_eef1a_records = integer(),
    eef1_complex_records = integer(),
    translation_network_records = integer(),
    matched_network_proteins =
      character()
  )
} else {
  target_summary <-
    candidate_target_evidence |>
    dplyr::group_by(
      .data$candidate_name
    ) |>
    dplyr::summarise(
      activity_record_count =
        dplyr::n(),

      unique_target_count =
        dplyr::n_distinct(
          .data$target_chembl_id,
          na.rm = TRUE
        ),

      direct_eef1a_records = sum(
        .data$target_evidence_class ==
          "direct_eef1a",
        na.rm = TRUE
      ),

      eef1_complex_records = sum(
        .data$target_evidence_class ==
          "eef1_complex",
        na.rm = TRUE
      ),

      translation_network_records =
        sum(
          .data$target_evidence_class ==
            "translation_network",
          na.rm = TRUE
        ),

      matched_network_proteins =
        collapse_unique(
          .data$matched_network_proteins
        ),

      .groups = "drop"
    )
}

candidate_classification <- candidates |>
  dplyr::select(
    candidate_name =
      "sider_drug_name",
    "sider_stitch_flat_id",
    "sider_pubchem_cid",
    "nearest_query_name",
    "nearest_query_tanimoto",
    "similarity_band",
    "unique_side_effect_count",
    "atc_codes"
  ) |>
  dplyr::left_join(
    candidate_annotations |>
      dplyr::select(
        "candidate_name",
        "molecule_chembl_id",
        "pref_name",
        "max_phase",
        "first_approval",
        "mapping_status"
      ),
    by = "candidate_name"
  ) |>
  dplyr::left_join(
    target_summary,
    by = "candidate_name"
  ) |>
  dplyr::mutate(
    activity_record_count =
      dplyr::coalesce(
        .data$activity_record_count,
        0L
      ),

    unique_target_count =
      dplyr::coalesce(
        .data$unique_target_count,
        0L
      ),

    direct_eef1a_records =
      dplyr::coalesce(
        .data$direct_eef1a_records,
        0L
      ),

    eef1_complex_records =
      dplyr::coalesce(
        .data$eef1_complex_records,
        0L
      ),

    translation_network_records =
      dplyr::coalesce(
        .data$translation_network_records,
        0L
      ),

    biological_classification =
      dplyr::case_when(
        .data$direct_eef1a_records > 0L ~
          "Direct eEF1A evidence",

        .data$eef1_complex_records > 0L ~
          "eEF1-complex evidence",

        .data$translation_network_records > 0L ~
          "Translation-network evidence",

        .data$activity_record_count > 0L ~
          paste(
            "Chemical similarity only;",
            "ChEMBL targets are outside",
            "the current eEF1A network"
          ),

        TRUE ~
          paste(
            "Chemical similarity only;",
            "no ChEMBL activities retrieved"
          )
      ),

    progression_status =
      dplyr::case_when(
        .data$direct_eef1a_records > 0L ~
          "advance_for_manual_review",

        .data$eef1_complex_records > 0L ~
          "advance_for_manual_review",

        .data$translation_network_records > 0L ~
          "advance_for_manual_review",

        TRUE ~
          "retain_as_chemical_comparator"
      )
  )

if (
  nrow(candidate_classification) != 6L
) {
  stop(
    "Candidate classification must contain six rows.",
    call. = FALSE
  )
}

# ============================================================
# 11. Create summary
# ============================================================

unique_target_count <- dplyr::n_distinct(
  candidate_activities$target_chembl_id,
  na.rm = TRUE
)

summary_table <- tibble::tibble(
  metric = c(
    "candidate_count",
    "chembl_mapped_candidate_count",
    "activity_record_count",
    "unique_chembl_target_count",
    "candidate_side_effect_record_count",
    "candidate_frequency_record_count",
    "candidate_indication_record_count",
    "side_effect_pairwise_comparison_count",
    "direct_eef1a_candidate_count",
    "eef1_complex_candidate_count",
    "translation_network_candidate_count",
    "chemical_only_candidate_count"
  ),
  value = c(
    as.character(
      nrow(candidates)
    ),
    as.character(
      mapped_candidate_count
    ),
    as.character(
      nrow(candidate_activities)
    ),
    as.character(
      unique_target_count
    ),
    as.character(
      nrow(candidate_side_effects)
    ),
    as.character(
      nrow(candidate_frequencies)
    ),
    as.character(
      nrow(candidate_indications)
    ),
    as.character(
      nrow(side_effect_similarity)
    ),
    as.character(
      sum(
        candidate_classification$
          biological_classification ==
          "Direct eEF1A evidence"
      )
    ),
    as.character(
      sum(
        candidate_classification$
          biological_classification ==
          "eEF1-complex evidence"
      )
    ),
    as.character(
      sum(
        candidate_classification$
          biological_classification ==
          "Translation-network evidence"
      )
    ),
    as.character(
      sum(
        grepl(
          "Chemical similarity only",
          candidate_classification$
            biological_classification
        )
      )
    )
  )
)

# ============================================================
# 12. Write CSV outputs
# ============================================================

cat(
  "Writing A17D outputs...\n"
)

safe_write_csv(
  candidate_annotations,
  outputs[["annotations"]]
)

safe_write_csv(
  candidate_activities,
  outputs[["activities"]]
)

safe_write_csv(
  candidate_target_evidence,
  outputs[["target_evidence"]]
)

safe_write_csv(
  candidate_side_effects,
  outputs[["side_effects"]]
)

safe_write_csv(
  candidate_frequencies,
  outputs[["frequencies"]]
)

safe_write_csv(
  candidate_indications,
  outputs[["indications"]]
)

safe_write_csv(
  side_effect_similarity,
  outputs[["side_effect_similarity"]]
)

safe_write_csv(
  candidate_classification,
  outputs[["classification"]]
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

# ============================================================
# 13. Create HTML report
# ============================================================

classification_report <- candidate_classification |>
  dplyr::select(
    "candidate_name",
    "nearest_query_name",
    "nearest_query_tanimoto",
    "activity_record_count",
    "unique_target_count",
    "matched_network_proteins",
    "biological_classification",
    "progression_status"
  )

similarity_report <- side_effect_similarity |>
  dplyr::select(
    "drug_1",
    "drug_2",
    "comparison_group",
    "shared_side_effect_count",
    "jaccard_similarity",
    "tfidf_cosine_similarity"
  )

report_html <- paste0(
  "<!doctype html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>RESKO A17D Candidate Evaluation</title>",
  "<style>",
  "body{font-family:Arial,sans-serif;max-width:1250px;",
  "margin:0 auto;padding:28px;line-height:1.5;color:#222}",
  "h1,h2{color:#263238}",
  ".note{background:#fff4e5;border-left:5px solid #e67e22;",
  "padding:14px;margin:16px 0}",
  "table{border-collapse:collapse;width:100%;font-size:13px}",
  "th,td{border:1px solid #ddd;padding:7px;text-align:left}",
  "th{background:#f3f4f6}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>RESKO A17D Biological and Side-Effect Evaluation</h1>",
  "<p>Six SIDER structural neighbours were evaluated using ",
  "ChEMBL activity records and local SIDER evidence.</p>",
  "<div class='note'>",
  "<strong>Interpretive boundary:</strong> ",
  "Chemical similarity and database annotations do not prove ",
  "direct eEF1A binding or inhibition.",
  "</div>",
  "<h2>Candidate classification</h2>",
  make_html_table(
    classification_report
  ),
  "<h2>Pairwise side-effect similarity</h2>",
  make_html_table(
    similarity_report
  ),
  "<h2>Method</h2>",
  "<p>ChEMBL compounds were mapped using standard InChIKeys. ",
  "Activities were retrieved using paginated API requests. ",
  "SIDER profiles used unique flat-drug and MedDRA preferred-term pairs. ",
  "Side-effect similarity used Jaccard and TF-IDF-weighted cosine scores.</p>",
  "<h2>Limitations</h2>",
  "<p>ChEMBL activities include heterogeneous assays and targets. ",
  "Candidate classifications require assay-level and document-level review. ",
  "SIDER 4.1 is historical and frequency coverage is incomplete.</p>",
  "</body>",
  "</html>"
)

temporary_report <- tempfile(
  pattern = "A17D_report_",
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
    "A17D HTML report was not created.",
    call. = FALSE
  )
}

backup_existing_file(
  outputs[["report"]]
)

report_moved <- file.rename(
  temporary_report,
  outputs[["report"]]
)

if (!report_moved) {
  stop(
    "Could not move the HTML report into place.",
    call. = FALSE
  )
}

# ============================================================
# 14. Verify outputs
# ============================================================

for (path in outputs) {
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

annotation_check <- readr::read_csv(
  outputs[["annotations"]],
  show_col_types = FALSE,
  progress = FALSE
)

classification_check <- readr::read_csv(
  outputs[["classification"]],
  show_col_types = FALSE,
  progress = FALSE
)

similarity_check <- readr::read_csv(
  outputs[["side_effect_similarity"]],
  show_col_types = FALSE,
  progress = FALSE
)

summary_check <- readr::read_csv(
  outputs[["summary"]],
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(annotation_check) != 6L) {
  stop(
    "Candidate annotations must contain six rows.",
    call. = FALSE
  )
}

if (
  nrow(classification_check) != 6L
) {
  stop(
    "Candidate classification must contain six rows.",
    call. = FALSE
  )
}

if (nrow(similarity_check) != 15L) {
  stop(
    "Side-effect similarity must contain 15 rows.",
    call. = FALSE
  )
}

if (nrow(summary_check) != 12L) {
  stop(
    "A17D summary must contain 12 rows.",
    call. = FALSE
  )
}

# ============================================================
# 15. Completion summary
# ============================================================

cat(
  "A17D combined candidate evaluation completed.\n"
)

cat(
  "Candidates evaluated: ",
  nrow(candidates),
  "\n",
  sep = ""
)

cat(
  "Candidates mapped to ChEMBL: ",
  mapped_candidate_count,
  "\n",
  sep = ""
)

cat(
  "ChEMBL activity records: ",
  nrow(candidate_activities),
  "\n",
  sep = ""
)

cat(
  "Unique ChEMBL targets: ",
  unique_target_count,
  "\n",
  sep = ""
)

cat(
  "SIDER side-effect records: ",
  nrow(candidate_side_effects),
  "\n",
  sep = ""
)

cat(
  "SIDER frequency records: ",
  nrow(candidate_frequencies),
  "\n",
  sep = ""
)

cat(
  "SIDER indication records: ",
  nrow(candidate_indications),
  "\n",
  sep = ""
)

cat(
  "Side-effect comparisons: ",
  nrow(side_effect_similarity),
  "\n",
  sep = ""
)

cat(
  "Candidates advancing for manual review: ",
  sum(
    candidate_classification$
      progression_status ==
      "advance_for_manual_review"
  ),
  "\n",
  sep = ""
)

cat(
  "Verified outputs: ",
  length(outputs),
  "\n",
  sep = ""
)