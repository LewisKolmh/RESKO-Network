#!/usr/bin/env Rscript

# ============================================================
# RESKO A17E: Provenance correction of SIDER candidates
# ============================================================
#
# Run from the RESKO project root:
#
# Rscript scripts/A17E_correct_sider_candidate_provenance.R
#
# Inputs:
#
# results/A17D_candidate_target_evidence.csv
# results/A17D_candidate_classification.csv
# results/A8_target_proteins.csv
#
# Alternative network input:
#
# results/nodes_proteins.csv
#
# Outputs:
#
# results/A17E_relevant_activity_records.csv
# results/A17E_assay_metadata.csv
# results/A17E_document_metadata.csv
# results/A17E_provenance_corrected_evidence.csv
# results/A17E_candidate_classification_corrected.csv
# results/A17E_correction_summary.csv
# results/A17E_provenance_report.html
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

# ============================================================
# 1. Validate packages
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
# 2. Project paths
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
  "A17E_api_cache"
)

if (
  !dir.exists(results_dir) ||
  !dir.exists(scripts_dir)
) {
  stop(
    "Run A17E from the RESKO project root.",
    call. = FALSE
  )
}

dir.create(
  cache_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

input_evidence <- file.path(
  results_dir,
  "A17D_candidate_target_evidence.csv"
)

input_classification <- file.path(
  results_dir,
  "A17D_candidate_classification.csv"
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

network_files <- network_options[
  file.exists(network_options)
]

if (length(network_files) == 0L) {
  stop(
    paste(
      "Neither results/A8_target_proteins.csv",
      "nor results/nodes_proteins.csv was found."
    ),
    call. = FALSE
  )
}

input_network <- network_files[1]

outputs <- c(
  relevant_records = file.path(
    results_dir,
    "A17E_relevant_activity_records.csv"
  ),
  assay_metadata = file.path(
    results_dir,
    "A17E_assay_metadata.csv"
  ),
  document_metadata = file.path(
    results_dir,
    "A17E_document_metadata.csv"
  ),
  corrected_evidence = file.path(
    results_dir,
    "A17E_provenance_corrected_evidence.csv"
  ),
  corrected_classification = file.path(
    results_dir,
    "A17E_candidate_classification_corrected.csv"
  ),
  summary = file.path(
    results_dir,
    "A17E_correction_summary.csv"
  ),
  report = file.path(
    results_dir,
    "A17E_provenance_report.html"
  )
)

required_inputs <- c(
  input_evidence,
  input_classification,
  input_network
)

for (path in required_inputs) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Missing or empty input: ",
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

safe_numeric <- function(x) {
  suppressWarnings(
    as.numeric(x)
  )
}

safe_logical <- function(x) {
  value <- tolower(
    clean_text(x)
  )

  ifelse(
    value %in% c(
      "true",
      "t",
      "1",
      "yes"
    ),
    TRUE,
    ifelse(
      value %in% c(
        "false",
        "f",
        "0",
        "no"
      ),
      FALSE,
      NA
    )
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

normalise_protein <- function(x) {
  output <- toupper(
    clean_text(x)
  )

  output <- gsub(
    "[^A-Z0-9]+",
    "",
    output
  )

  output[
    output == ""
  ] <- NA_character_

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
    path,
    backup_path,
    overwrite = FALSE
  )

  if (!copied) {
    stop(
      "Could not preserve previous output: ",
      path,
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

  output_data <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  list_columns <- names(output_data)[
    vapply(
      output_data,
      is.list,
      FUN.VALUE = logical(1)
    )
  ]

  if (length(list_columns) > 0L) {
    stop(
      "Output contains list column(s): ",
      paste(list_columns, collapse = ", "),
      call. = FALSE
    )
  }

  readr::write_csv(
    output_data,
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

  invisible(path)
}

request_json <- function(
  url,
  cache_name,
  attempts = 3L
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
      error = function(error) NULL
    )

    if (!is.null(cached)) {
      return(cached)
    }

    unlink(cache_path)
  }

  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_user_agent(
          "RESKO-A17E/3.0"
        ) |>
        httr2::req_timeout(
          seconds = 120
        ) |>
        httr2::req_perform(),
      error = function(error) NULL
    )

    if (!is.null(response)) {
      status <- httr2::resp_status(
        response
      )

      if (
        status >= 200L &&
        status < 300L
      ) {
        response_body <- httr2::resp_body_string(
          response
        )

        parsed <- tryCatch(
          jsonlite::fromJSON(
            response_body,
            simplifyDataFrame = TRUE
          ),
          error = function(error) NULL
        )

        if (!is.null(parsed)) {
          writeLines(
            response_body,
            cache_path,
            useBytes = TRUE
          )

          Sys.sleep(0.25)

          return(parsed)
        }
      }
    }

    if (attempt < attempts) {
      Sys.sleep(
        min(
          2 ^ attempt,
          8
        )
      )
    }
  }

  NULL
}

get_value <- function(
  object,
  field,
  default = NA_character_
) {
  if (is.null(object)) {
    return(default)
  }

  if (
    is.data.frame(object) &&
    field %in% names(object)
  ) {
    value <- object[[field]][1]
  } else if (
    is.list(object) &&
    !is.null(object[[field]])
  ) {
    value <- object[[field]][1]
  } else {
    return(default)
  }

  if (
    length(value) == 0L ||
    is.null(value) ||
    is.na(value)
  ) {
    return(default)
  }

  as.character(value)
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
      paste0(
        "<tr>",
        paste0(
          "<td>",
          unlist(
            table_data[
              row_index,
              ,
              drop = FALSE
            ]
          ),
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
# 4. Read A17D evidence
# ============================================================

cat(
  "Reading A17E inputs...\n"
)

evidence <- readr::read_csv(
  input_evidence,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

classification <- readr::read_csv(
  input_classification,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

network <- readr::read_csv(
  input_network,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

required_evidence_columns <- c(
  "candidate_name",
  "molecule_chembl_id",
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
  "potential_duplicate",
  "matched_network_proteins",
  "target_evidence_class"
)

missing_evidence_columns <- setdiff(
  required_evidence_columns,
  names(evidence)
)

if (length(missing_evidence_columns) > 0L) {
  stop(
    "A17D evidence is missing column(s): ",
    paste(
      missing_evidence_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

expected_candidates <- c(
  "nilotinib",
  "imatinib",
  "alprazolam",
  "triazolam",
  "ponatinib",
  "temazepam"
)

observed_candidates <- sort(
  unique(
    tolower(
      clean_text(
        classification$candidate_name
      )
    )
  )
)

missing_candidates <- setdiff(
  expected_candidates,
  observed_candidates
)

if (length(missing_candidates) > 0L) {
  stop(
    "A17D classification is missing candidate(s): ",
    paste(
      missing_candidates,
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
  names(network)
)

if (length(protein_columns) == 0L) {
  stop(
    "Could not identify the network protein column.",
    call. = FALSE
  )
}

network_proteins <- unique(
  normalise_protein(
    network[[
      protein_columns[1]
    ]]
  )
)

network_proteins <- network_proteins[
  !is.na(network_proteins)
]

# ============================================================
# 5. Recalculate network membership
# ============================================================

cat(
  "Recalculating candidate-network matches...\n"
)

target_names <- clean_text(
  evidence$target_pref_name
)

normalised_targets <- normalise_protein(
  target_names
)

match_network_target <- function(target) {
  if (is.na(target)) {
    return(NA_character_)
  }

  matches <- network_proteins[
    vapply(
      network_proteins,
      function(protein) {
        identical(
          protein,
          target
        )
      },
      FUN.VALUE = logical(1)
    )
  ]

  collapse_unique(matches)
}

matched_network <- vapply(
  normalised_targets,
  match_network_target,
  FUN.VALUE = character(1)
)

target_upper <- toupper(
  dplyr::coalesce(
    target_names,
    ""
  )
)

recalculated_class <- dplyr::case_when(
  grepl(
    "EEF1A1|EEF1A2",
    target_upper
  ) ~ "direct_eef1a",

  grepl(
    "EEF1B2|EEF1D|EEF1G",
    target_upper
  ) ~ "eef1_complex",

  !is.na(
    matched_network
  ) ~ "translation_network",

  TRUE ~ "other_target"
)

evidence_corrected <- evidence |>
  dplyr::mutate(
    standard_value =
      safe_numeric(
        .data$standard_value
      ),

    pchembl_value =
      safe_numeric(
        .data$pchembl_value
      ),

    potential_duplicate =
      safe_logical(
        .data$potential_duplicate
      ),

    original_network_match =
      clean_text(
        .data$matched_network_proteins
      ),

    original_evidence_class =
      clean_text(
        .data$target_evidence_class
      ),

    matched_network_proteins =
      matched_network,

    target_evidence_class =
      recalculated_class,

    classification_changed =
      dplyr::coalesce(
        .data$original_evidence_class,
        ""
      ) !=
      dplyr::coalesce(
        .data$target_evidence_class,
        ""
      )
  )

relevant_records <- evidence_corrected |>
  dplyr::filter(
    .data$original_evidence_class !=
      "other_target" |
      .data$target_evidence_class !=
        "other_target"
  ) |>
  dplyr::arrange(
    .data$candidate_name,
    .data$target_chembl_id,
    .data$assay_chembl_id
  )

if (nrow(relevant_records) == 0L) {
  stop(
    "No potentially relevant A17D records were found.",
    call. = FALSE
  )
}

# ============================================================
# 6. Retrieve assay metadata
# ============================================================

cat(
  "Retrieving assay metadata...\n"
)

assay_ids <- sort(
  unique(
    clean_text(
      relevant_records$assay_chembl_id
    )
  )
)

assay_ids <- assay_ids[
  !is.na(assay_ids)
]

assay_rows <- vector(
  "list",
  length(assay_ids)
)

for (assay_index in seq_along(assay_ids)) {
  assay_id <- assay_ids[
    assay_index
  ]

  cat(
    "  Assay ",
    assay_index,
    " of ",
    length(assay_ids),
    ": ",
    assay_id,
    "\n",
    sep = ""
  )

  filtered_url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/assay.json",
    "?assay_chembl_id=",
    utils::URLencode(
      assay_id,
      reserved = TRUE
    ),
    "&limit=10"
  )

  filtered <- request_json(
    filtered_url,
    paste0(
      "assay_",
      assay_id,
      "_filtered.json"
    ),
    attempts = 3L
  )

  assay_object <- NULL
  retrieval_status <- "metadata_unavailable"

  if (
    !is.null(filtered) &&
    !is.null(filtered$assays) &&
    is.data.frame(filtered$assays) &&
    nrow(filtered$assays) > 0L
  ) {
    assay_object <- filtered$assays[
      1,
      ,
      drop = FALSE
    ]

    retrieval_status <- "filtered_endpoint"
  }

  assay_rows[[
    assay_index
  ]] <- tibble::tibble(
    assay_chembl_id = assay_id,

    assay_type =
      get_value(
        assay_object,
        "assay_type"
      ),

    assay_description =
      get_value(
        assay_object,
        "description"
      ),

    assay_organism =
      get_value(
        assay_object,
        "assay_organism"
      ),

    assay_tax_id =
      safe_numeric(
        get_value(
          assay_object,
          "assay_tax_id"
        )
      ),

    assay_test_type =
      get_value(
        assay_object,
        "assay_test_type"
      ),

    assay_category =
      get_value(
        assay_object,
        "assay_category"
      ),

    bao_format =
      get_value(
        assay_object,
        "bao_format"
      ),

    bao_label =
      get_value(
        assay_object,
        "bao_label"
      ),

    confidence_score =
      safe_numeric(
        get_value(
          assay_object,
          "confidence_score"
        )
      ),

    relationship_type =
      get_value(
        assay_object,
        "relationship_type"
      ),

    target_chembl_id =
      get_value(
        assay_object,
        "target_chembl_id"
      ),

    document_chembl_id =
      get_value(
        assay_object,
        "document_chembl_id"
      ),

    assay_metadata_status =
      retrieval_status
  )
}

assay_metadata <- dplyr::bind_rows(
  assay_rows
)

if (nrow(assay_metadata) == 0L) {
  stop(
    "No assay metadata rows were created.",
    call. = FALSE
  )
}

# ============================================================
# 7. Retrieve document metadata
# ============================================================

cat(
  "Retrieving document metadata...\n"
)

document_ids <- sort(
  unique(
    c(
      clean_text(
        relevant_records$document_chembl_id
      ),
      clean_text(
        assay_metadata$document_chembl_id
      )
    )
  )
)

document_ids <- document_ids[
  !is.na(document_ids)
]

document_rows <- vector(
  "list",
  length(document_ids)
)

for (
  document_index in seq_along(
    document_ids
  )
) {
  document_id <- document_ids[
    document_index
  ]

  cat(
    "  Document ",
    document_index,
    " of ",
    length(document_ids),
    ": ",
    document_id,
    "\n",
    sep = ""
  )

  filtered_url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/document.json",
    "?document_chembl_id=",
    utils::URLencode(
      document_id,
      reserved = TRUE
    ),
    "&limit=10"
  )

  filtered <- request_json(
    filtered_url,
    paste0(
      "document_",
      document_id,
      "_filtered.json"
    ),
    attempts = 3L
  )

  document_object <- NULL
  retrieval_status <- "metadata_unavailable"

  if (
    !is.null(filtered) &&
    !is.null(filtered$documents) &&
    is.data.frame(filtered$documents) &&
    nrow(filtered$documents) > 0L
  ) {
    document_object <- filtered$documents[
      1,
      ,
      drop = FALSE
    ]

    retrieval_status <- "filtered_endpoint"
  }

  document_rows[[
    document_index
  ]] <- tibble::tibble(
    document_chembl_id =
      document_id,

    title =
      get_value(
        document_object,
        "title"
      ),

    journal =
      get_value(
        document_object,
        "journal"
      ),

    year =
      safe_numeric(
        get_value(
          document_object,
          "year"
        )
      ),

    volume =
      get_value(
        document_object,
        "volume"
      ),

    first_page =
      get_value(
        document_object,
        "first_page"
      ),

    last_page =
      get_value(
        document_object,
        "last_page"
      ),

    doi =
      get_value(
        document_object,
        "doi"
      ),

    pubmed_id =
      get_value(
        document_object,
        "pubmed_id"
      ),

    document_type =
      get_value(
        document_object,
        "doc_type"
      ),

    authors =
      get_value(
        document_object,
        "authors"
      ),

    document_metadata_status =
      retrieval_status
  )
}

document_metadata <- dplyr::bind_rows(
  document_rows
)

if (nrow(document_metadata) == 0L) {
  stop(
    "No document metadata rows were created.",
    call. = FALSE
  )
}

# ============================================================
# 8. Construct record-level provenance
# ============================================================

cat(
  "Constructing provenance-corrected evidence...\n"
)

record_level <- relevant_records |>
  dplyr::left_join(
    assay_metadata,
    by = "assay_chembl_id",
    suffix = c("", "_assay")
  ) |>
  dplyr::left_join(
    document_metadata,
    by = "document_chembl_id"
  ) |>
  dplyr::mutate(
    human_assay =
      grepl(
        "Homo sapiens|Human",
        dplyr::coalesce(
          .data$assay_organism,
          ""
        ),
        ignore.case = TRUE
      ),

    weak_30_micromolar_record =
      .data$standard_units == "nM" &
      !is.na(
        .data$standard_value
      ) &
      .data$standard_value >= 30000,

    validity_warning_present =
      !is.na(
        clean_text(
          .data$data_validity_comment
        )
      ),

    potential_duplicate =
      dplyr::coalesce(
        .data$potential_duplicate,
        FALSE
      ),

    assay_metadata_available =
      dplyr::coalesce(
        .data$assay_metadata_status,
        "metadata_unavailable"
      ) != "metadata_unavailable",

    document_metadata_available =
      dplyr::coalesce(
        .data$document_metadata_status,
        "metadata_unavailable"
      ) != "metadata_unavailable"
  )

# ============================================================
# 9. Collapse duplicate provenance representations
# ============================================================

corrected_evidence <- record_level |>
  dplyr::group_by(
    .data$candidate_name,
    .data$molecule_chembl_id,
    .data$target_chembl_id,
    .data$target_pref_name,
    .data$target_organism,
    .data$matched_network_proteins,
    .data$target_evidence_class,
    .data$assay_chembl_id,
    .data$document_chembl_id,
    .data$standard_type,
    .data$standard_relation,
    .data$standard_value,
    .data$standard_units
  ) |>
  dplyr::summarise(
    raw_activity_record_count =
      dplyr::n(),

    activity_ids =
      collapse_unique(
        .data$activity_id
      ),

    pchembl_values =
      collapse_unique(
        .data$pchembl_value
      ),

    assay_type =
      collapse_unique(
        .data$assay_type
      ),

    assay_description =
      collapse_unique(
        .data$assay_description
      ),

    assay_organism =
      collapse_unique(
        .data$assay_organism
      ),

    confidence_score =
      suppressWarnings(
        max(
          .data$confidence_score,
          na.rm = TRUE
        )
      ),

    document_title =
      collapse_unique(
        .data$title
      ),

    doi =
      collapse_unique(
        .data$doi
      ),

    pubmed_id =
      collapse_unique(
        .data$pubmed_id
      ),

    potential_duplicate_present =
      any(
        .data$potential_duplicate,
        na.rm = TRUE
      ),

    validity_warning_present =
      any(
        .data$validity_warning_present,
        na.rm = TRUE
      ),

    human_assay =
      any(
        .data$human_assay,
        na.rm = TRUE
      ),

    weak_30_micromolar_record =
      any(
        .data$weak_30_micromolar_record,
        na.rm = TRUE
      ),

    assay_metadata_available =
      any(
        .data$assay_metadata_available,
        na.rm = TRUE
      ),

    document_metadata_available =
      any(
        .data$document_metadata_available,
        na.rm = TRUE
      ),

    .groups = "drop"
  ) |>
  dplyr::mutate(
    confidence_score =
      ifelse(
        is.infinite(
          .data$confidence_score
        ),
        NA_real_,
        .data$confidence_score
      ),

    independent_assay_count =
      ifelse(
        is.na(
          .data$assay_chembl_id
        ),
        0L,
        1L
      ),

    independent_document_count =
      ifelse(
        is.na(
          .data$document_chembl_id
        ),
        0L,
        1L
      ),

    provenance_interpretation =
      dplyr::case_when(
        .data$target_evidence_class ==
          "direct_eef1a" ~
          paste(
            "Direct eEF1A annotation;",
            "manual assay review required"
          ),

        .data$target_evidence_class ==
          "eef1_complex" ~
          paste(
            "eEF1-complex annotation;",
            "manual assay review required"
          ),

        .data$target_evidence_class ==
          "translation_network" &
          .data$weak_30_micromolar_record ~
          paste(
            "Weak translation-network binding annotation;",
            "not direct eEF1A evidence"
          ),

        .data$target_evidence_class ==
          "translation_network" ~
          paste(
            "Translation-network annotation;",
            "not direct eEF1A evidence"
          ),

        TRUE ~
          paste(
            "Chemical similarity only;",
            "network relationship not confirmed"
          )
      )
  )

if (nrow(corrected_evidence) == 0L) {
  stop(
    "No corrected evidence relationships were created.",
    call. = FALSE
  )
}

# ============================================================
# 10. Candidate-level corrected classification
# ============================================================

candidate_evidence_summary <- corrected_evidence |>
  dplyr::group_by(
    .data$candidate_name
  ) |>
  dplyr::summarise(
    corrected_relationship_count =
      dplyr::n(),

    independent_assay_count =
      dplyr::n_distinct(
        .data$assay_chembl_id,
        na.rm = TRUE
      ),

    independent_document_count =
      dplyr::n_distinct(
        .data$document_chembl_id,
        na.rm = TRUE
      ),

    direct_eef1a_relationships =
      sum(
        .data$target_evidence_class ==
          "direct_eef1a",
        na.rm = TRUE
      ),

    eef1_complex_relationships =
      sum(
        .data$target_evidence_class ==
          "eef1_complex",
        na.rm = TRUE
      ),

    translation_network_relationships =
      sum(
        .data$target_evidence_class ==
          "translation_network",
        na.rm = TRUE
      ),

    weak_30_micromolar_relationships =
      sum(
        .data$weak_30_micromolar_record,
        na.rm = TRUE
      ),

    matched_network_proteins =
      collapse_unique(
        .data$matched_network_proteins
      ),

    .groups = "drop"
  )

classification_for_join <- classification |>
  dplyr::rename(
    matched_network_proteins_A17D =
      "matched_network_proteins"
  )

corrected_classification <- classification_for_join |>
  dplyr::left_join(
    candidate_evidence_summary,
    by = "candidate_name"
  ) |>
  dplyr::mutate(
    corrected_relationship_count =
      dplyr::coalesce(
        safe_numeric(
          .data$corrected_relationship_count
        ),
        0
      ),

    independent_assay_count =
      dplyr::coalesce(
        safe_numeric(
          .data$independent_assay_count
        ),
        0
      ),

    independent_document_count =
      dplyr::coalesce(
        safe_numeric(
          .data$independent_document_count
        ),
        0
      ),

    direct_eef1a_relationships =
      dplyr::coalesce(
        safe_numeric(
          .data$direct_eef1a_relationships
        ),
        0
      ),

    eef1_complex_relationships =
      dplyr::coalesce(
        safe_numeric(
          .data$eef1_complex_relationships
        ),
        0
      ),

    translation_network_relationships =
      dplyr::coalesce(
        safe_numeric(
          .data$translation_network_relationships
        ),
        0
      ),

    weak_30_micromolar_relationships =
      dplyr::coalesce(
        safe_numeric(
          .data$weak_30_micromolar_relationships
        ),
        0
      ),

    matched_network_proteins =
      clean_text(
        .data$matched_network_proteins
      ),

    corrected_biological_classification =
      dplyr::case_when(
        .data$direct_eef1a_relationships > 0 ~
          paste(
            "Direct eEF1A evidence;",
            "manual provenance review required"
          ),

        .data$eef1_complex_relationships > 0 ~
          paste(
            "eEF1-complex evidence;",
            "manual provenance review required"
          ),

        .data$translation_network_relationships > 0 &
          .data$weak_30_micromolar_relationships > 0 ~
          paste(
            "Weak translation-network binding annotation;",
            "not direct eEF1A evidence"
          ),

        .data$translation_network_relationships > 0 ~
          paste(
            "Translation-network annotation;",
            "not direct eEF1A evidence"
          ),

        TRUE ~
          "Chemical and side-effect comparator only"
      ),

    corrected_progression_status =
      dplyr::case_when(
        .data$direct_eef1a_relationships > 0 ~
          "advance_for_manual_provenance_review",

        .data$eef1_complex_relationships > 0 ~
          "advance_for_manual_provenance_review",

        .data$translation_network_relationships > 0 ~
          "retain_as_low_priority_network_comparator",

        TRUE ~
          "retain_as_chemical_comparator"
      )
  )

if (nrow(corrected_classification) != 6L) {
  stop(
    "Corrected classification must contain six candidates.",
    call. = FALSE
  )
}

# ============================================================
# 11. Summary calculations
# ============================================================

independent_assays <- dplyr::n_distinct(
  corrected_evidence$assay_chembl_id,
  na.rm = TRUE
)

independent_documents <- dplyr::n_distinct(
  corrected_evidence$document_chembl_id,
  na.rm = TRUE
)

classification_changes <- sum(
  evidence_corrected$classification_changed,
  na.rm = TRUE
)

weak_relationship_count <- sum(
  corrected_evidence$weak_30_micromolar_record,
  na.rm = TRUE
)

direct_candidate_count <- sum(
  corrected_classification$
    direct_eef1a_relationships > 0,
  na.rm = TRUE
)

complex_candidate_count <- sum(
  corrected_classification$
    eef1_complex_relationships > 0,
  na.rm = TRUE
)

network_candidate_count <- sum(
  corrected_classification$
    translation_network_relationships > 0,
  na.rm = TRUE
)

chemical_comparator_count <- sum(
  corrected_classification$
    corrected_progression_status ==
    "retain_as_chemical_comparator",
  na.rm = TRUE
)

low_priority_comparator_count <- sum(
  corrected_classification$
    corrected_progression_status ==
    "retain_as_low_priority_network_comparator",
  na.rm = TRUE
)

assay_metadata_unavailable <- sum(
  assay_metadata$assay_metadata_status ==
    "metadata_unavailable",
  na.rm = TRUE
)

document_metadata_unavailable <- sum(
  document_metadata$
    document_metadata_status ==
    "metadata_unavailable",
  na.rm = TRUE
)

summary_table <- tibble::tibble(
  metric = c(
    "input_activity_record_count",
    "relevant_input_record_count",
    "recalculated_network_record_count",
    "classification_changed_record_count",
    "corrected_evidence_relationship_count",
    "candidate_count_with_corrected_evidence",
    "independent_assay_count",
    "independent_document_count",
    "weak_30000_nM_relationship_count",
    "direct_eef1a_candidate_count",
    "eef1_complex_candidate_count",
    "translation_network_candidate_count",
    "chemical_comparator_count",
    "low_priority_network_comparator_count",
    "assay_metadata_unavailable_count",
    "document_metadata_unavailable_count"
  ),
  value = c(
    as.character(
      nrow(evidence)
    ),

    as.character(
      nrow(relevant_records)
    ),

    as.character(
      sum(
        evidence_corrected$
          target_evidence_class !=
          "other_target",
        na.rm = TRUE
      )
    ),

    as.character(
      classification_changes
    ),

    as.character(
      nrow(corrected_evidence)
    ),

    as.character(
      dplyr::n_distinct(
        corrected_evidence$candidate_name
      )
    ),

    as.character(
      independent_assays
    ),

    as.character(
      independent_documents
    ),

    as.character(
      weak_relationship_count
    ),

    as.character(
      direct_candidate_count
    ),

    as.character(
      complex_candidate_count
    ),

    as.character(
      network_candidate_count
    ),

    as.character(
      chemical_comparator_count
    ),

    as.character(
      low_priority_comparator_count
    ),

    as.character(
      assay_metadata_unavailable
    ),

    as.character(
      document_metadata_unavailable
    )
  )
)

# ============================================================
# 12. Write CSV outputs
# ============================================================

cat(
  "Writing A17E outputs...\n"
)

safe_write_csv(
  relevant_records,
  outputs[["relevant_records"]]
)

safe_write_csv(
  assay_metadata,
  outputs[["assay_metadata"]]
)

safe_write_csv(
  document_metadata,
  outputs[["document_metadata"]]
)

safe_write_csv(
  corrected_evidence,
  outputs[["corrected_evidence"]]
)

safe_write_csv(
  corrected_classification,
  outputs[["corrected_classification"]]
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

# ============================================================
# 13. Create HTML report
# ============================================================

classification_report <- corrected_classification |>
  dplyr::select(
    "candidate_name",
    "nearest_query_name",
    "nearest_query_tanimoto",
    "independent_assay_count",
    "independent_document_count",
    "matched_network_proteins",
    "corrected_biological_classification",
    "corrected_progression_status"
  )

evidence_report <- corrected_evidence |>
  dplyr::select(
    "candidate_name",
    "target_pref_name",
    "matched_network_proteins",
    "standard_type",
    "standard_value",
    "standard_units",
    "assay_chembl_id",
    "document_chembl_id",
    "assay_type",
    "confidence_score",
    "weak_30_micromolar_record",
    "provenance_interpretation"
  )

report_html <- paste0(
  "<!doctype html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width,initial-scale=1'>",
  "<title>RESKO A17E Provenance Correction</title>",
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
  "<h1>RESKO A17E Provenance Correction</h1>",
  "<p>Provisional A17D network relationships were re-evaluated ",
  "using exact network membership, assay metadata, document metadata, ",
  "and independent evidence counts.</p>",
  "<div class='note'>",
  "<strong>Interpretive boundary:</strong> ",
  "A weak network-protein binding annotation is not direct evidence ",
  "of EEF1A1 or EEF1A2 engagement.",
  "</div>",
  "<h2>Corrected candidate classification</h2>",
  make_html_table(
    classification_report
  ),
  "<h2>Corrected evidence relationships</h2>",
  make_html_table(
    evidence_report
  ),
  "<h2>Assay metadata</h2>",
  make_html_table(
    assay_metadata
  ),
  "<h2>Document metadata</h2>",
  make_html_table(
    document_metadata
  ),
  "<h2>Method</h2>",
  "<p>A17D evidence was read using character-only parsing. ",
  "Network membership was recalculated using exact normalised protein names. ",
  "Assay and document metadata were retrieved using filtered ChEMBL endpoints. ",
  "Evidence was grouped by compound, target, assay, document, activity type, ",
  "value, and units.</p>",
  "<h2>Limitations</h2>",
  "<p>Database annotations do not demonstrate cellular target engagement, ",
  "functional inhibition, selectivity, clinical efficacy, or causal ",
  "relationships with reported side effects.</p>",
  "</body>",
  "</html>"
)

temporary_report <- tempfile(
  pattern = "A17E_report_",
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
    "A17E HTML report was not created.",
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
    "Could not move the A17E report into place.",
    call. = FALSE
  )
}

# ============================================================
# 14. Validate all outputs
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

relevant_check <- readr::read_csv(
  outputs[["relevant_records"]],
  show_col_types = FALSE,
  progress = FALSE
)

assay_check <- readr::read_csv(
  outputs[["assay_metadata"]],
  show_col_types = FALSE,
  progress = FALSE
)

document_check <- readr::read_csv(
  outputs[["document_metadata"]],
  show_col_types = FALSE,
  progress = FALSE
)

evidence_check <- readr::read_csv(
  outputs[["corrected_evidence"]],
  show_col_types = FALSE,
  progress = FALSE
)

classification_check <- readr::read_csv(
  outputs[["corrected_classification"]],
  show_col_types = FALSE,
  progress = FALSE
)

summary_check <- readr::read_csv(
  outputs[["summary"]],
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(relevant_check) == 0L) {
  stop(
    "Relevant-record output is empty.",
    call. = FALSE
  )
}

if (nrow(assay_check) == 0L) {
  stop(
    "Assay-metadata output is empty.",
    call. = FALSE
  )
}

if (nrow(document_check) == 0L) {
  stop(
    "Document-metadata output is empty.",
    call. = FALSE
  )
}

if (nrow(evidence_check) == 0L) {
  stop(
    "Corrected-evidence output is empty.",
    call. = FALSE
  )
}

if (nrow(classification_check) != 6L) {
  stop(
    "Corrected classification must contain six candidates.",
    call. = FALSE
  )
}

if (nrow(summary_check) != 16L) {
  stop(
    "A17E summary must contain 16 rows.",
    call. = FALSE
  )
}

# ============================================================
# 15. Display corrected candidate classifications
# ============================================================

classification_display <- classification_check |>
  dplyr::select(
    "candidate_name",
    "nearest_query_name",
    "nearest_query_tanimoto",
    "matched_network_proteins",
    "independent_assay_count",
    "independent_document_count",
    "corrected_biological_classification",
    "corrected_progression_status"
  ) |>
  dplyr::arrange(
    .data$candidate_name
  )

cat(
  "\nCorrected candidate classifications:\n"
)

print(
  classification_display,
  n = Inf,
  width = Inf
)

# ============================================================
# 16. Completion summary
# ============================================================

cat(
  "\nA17E provenance correction completed.\n"
)

cat(
  "Relevant input records: ",
  nrow(relevant_records),
  "\n",
  sep = ""
)

cat(
  "Corrected evidence relationships: ",
  nrow(corrected_evidence),
  "\n",
  sep = ""
)

cat(
  "Independent assays: ",
  independent_assays,
  "\n",
  sep = ""
)

cat(
  "Independent documents: ",
  independent_documents,
  "\n",
  sep = ""
)

cat(
  "Classification changes: ",
  classification_changes,
  "\n",
  sep = ""
)

cat(
  "Weak 30 micromolar relationships: ",
  weak_relationship_count,
  "\n",
  sep = ""
)

cat(
  "Direct eEF1A candidates: ",
  direct_candidate_count,
  "\n",
  sep = ""
)

cat(
  "eEF1-complex candidates: ",
  complex_candidate_count,
  "\n",
  sep = ""
)

cat(
  "Translation-network candidates: ",
  network_candidate_count,
  "\n",
  sep = ""
)

cat(
  "Low-priority network comparators: ",
  low_priority_comparator_count,
  "\n",
  sep = ""
)

cat(
  "Chemical comparators: ",
  chemical_comparator_count,
  "\n",
  sep = ""
)

cat(
  "Assays with unavailable metadata: ",
  assay_metadata_unavailable,
  "\n",
  sep = ""
)

cat(
  "Documents with unavailable metadata: ",
  document_metadata_unavailable,
  "\n",
  sep = ""
)

cat(
  "Verified outputs: ",
  length(outputs),
  "\n",
  sep = ""
)

cat(
  "A17E report: ",
  outputs[["report"]],
  "\n",
  sep = ""
)