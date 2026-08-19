#!/usr/bin/env Rscript

# ============================================================
# RESKO A16A: Retrieve reference structures
# ============================================================
# Purpose:
#   Combine six A11 ChEMBL candidate structures with structures
#   for four established eEF1A reference ligands retrieved from
#   PubChem PUG REST.
#
# Run from the RESKO project root with:
#   Rscript scripts/A16A_retrieve_reference_structures.R
#
# Inputs:
#   results/A11_candidate_molecule_annotations.csv
#
# Outputs:
#   results/A16_reference_compounds.csv
#   results/A16_structures_combined.csv
#   results/A16_structure_retrieval_summary.csv
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

# -----------------------------
# 1. Package checks
# -----------------------------
required_packages <- c("readr", "dplyr", "tibble", "httr2", "jsonlite")

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
    ". Install them from the project terminal with: Rscript -e ",
    "'install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "), repos=\"https://cloud.r-project.org\")'",
    call. = FALSE
  )
}

# -----------------------------
# 2. Project paths
# -----------------------------
project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
input_file <- file.path(results_dir, "A11_candidate_molecule_annotations.csv")

output_reference <- file.path(results_dir, "A16_reference_compounds.csv")
output_combined <- file.path(results_dir, "A16_structures_combined.csv")
output_summary <- file.path(results_dir, "A16_structure_retrieval_summary.csv")

if (!dir.exists(file.path(project_root, "scripts"))) {
  stop(
    "The scripts directory was not found. Run this script from the RESKO project root.",
    call. = FALSE
  )
}

dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(input_file)) {
  stop("Required input file not found: ", input_file, call. = FALSE)
}

input_size <- file.info(input_file)$size
if (is.na(input_size) || input_size <= 0L) {
  stop("Required input file is empty: ", input_file, call. = FALSE)
}

# -----------------------------
# 3. Helper functions
# -----------------------------
clean_text <- function(x) {
  output <- as.character(x)
  output[is.na(output) | trimws(output) == ""] <- NA_character_
  output
}

get_property <- function(record, names_to_try) {
  if (is.null(record)) {
    return(NA_character_)
  }

  for (property_name in names_to_try) {
    if (property_name %in% names(record)) {
      value <- record[[property_name]][1]
      if (length(value) > 0L && !is.na(value) && trimws(as.character(value)) != "") {
        return(as.character(value))
      }
    }
  }

  NA_character_
}

backup_existing_file <- function(path) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  backup_path <- paste0(stem, "_previous_", timestamp, ".", extension)

  copied <- file.copy(path, backup_path, overwrite = FALSE)
  if (!copied || !file.exists(backup_path) || file.info(backup_path)$size <= 0L) {
    stop("Could not preserve previous output: ", path, call. = FALSE)
  }

  invisible(backup_path)
}

safe_write_csv <- function(data, path) {
  temp_path <- paste0(path, ".tmp")

  if (file.exists(temp_path)) {
    unlink(temp_path)
  }

  data_to_write <- as.data.frame(data, stringsAsFactors = FALSE)

  list_columns <- names(data_to_write)[vapply(data_to_write, is.list, logical(1))]
  if (length(list_columns) > 0L) {
    stop(
      "Refusing to write list-column(s) to CSV: ",
      paste(list_columns, collapse = ", "),
      call. = FALSE
    )
  }

  readr::write_csv(data_to_write, temp_path, na = "")

  if (!file.exists(temp_path) || is.na(file.info(temp_path)$size) || file.info(temp_path)$size <= 0L) {
    stop("Failed to create non-empty temporary output: ", temp_path, call. = FALSE)
  }

  backup_existing_file(path)

  if (!file.rename(temp_path, path)) {
    stop("Could not move temporary output into place: ", path, call. = FALSE)
  }

  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Output verification failed: ", path, call. = FALSE)
  }

  invisible(path)
}

request_pubchem <- function(query_name, max_attempts = 5L) {
  property_names <- paste(
    c(
      "Title",
      "IUPACName",
      "CanonicalSMILES",
      "IsomericSMILES",
      "InChI",
      "InChIKey",
      "MolecularFormula",
      "MolecularWeight"
    ),
    collapse = ","
  )

  encoded_name <- utils::URLencode(query_name, reserved = TRUE)
  request_url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/",
    encoded_name,
    "/property/",
    property_names,
    "/JSON"
  )

  last_message <- "Request was not attempted"

  for (attempt in seq_len(max_attempts)) {
    response <- tryCatch(
      httr2::request(request_url) |>
        httr2::req_user_agent("RESKO-A16A/1.1 scientific-structure-retrieval") |>
        httr2::req_timeout(seconds = 45) |>
        httr2::req_perform(),
      error = function(error_condition) error_condition
    )

    if (!inherits(response, "error")) {
      status_code <- httr2::resp_status(response)

      if (status_code >= 200L && status_code < 300L) {
        parsed <- tryCatch(
          jsonlite::fromJSON(
            httr2::resp_body_string(response),
            simplifyDataFrame = TRUE
          ),
          error = function(error_condition) error_condition
        )

        if (
          !inherits(parsed, "error") &&
          !is.null(parsed$PropertyTable$Properties) &&
          is.data.frame(parsed$PropertyTable$Properties) &&
          nrow(parsed$PropertyTable$Properties) >= 1L
        ) {
          return(list(
            ok = TRUE,
            data = parsed$PropertyTable$Properties[1, , drop = FALSE],
            query = query_name,
            message = "retrieved"
          ))
        }

        last_message <- "PubChem returned no usable property record"
      } else if (status_code == 404L) {
        return(list(
          ok = FALSE,
          data = NULL,
          query = query_name,
          message = "not_found"
        ))
      } else {
        last_message <- paste0("HTTP ", status_code)
      }
    } else {
      last_message <- conditionMessage(response)
    }

    if (attempt < max_attempts) {
      Sys.sleep(min(2^(attempt - 1L), 16L))
    }
  }

  list(
    ok = FALSE,
    data = NULL,
    query = query_name,
    message = last_message
  )
}

# -----------------------------
# 4. Read and validate A11 candidates
# -----------------------------
candidates_raw <- tryCatch(
  readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE),
  error = function(error_condition) {
    stop("Failed to read A11 input: ", conditionMessage(error_condition), call. = FALSE)
  }
)

required_candidate_columns <- c(
  "molecule_chembl_id",
  "display_name",
  "canonical_smiles",
  "standard_inchi_key"
)

missing_columns <- setdiff(required_candidate_columns, names(candidates_raw))
if (length(missing_columns) > 0L) {
  stop(
    "A11 input is missing required column(s): ",
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

expected_candidate_ids <- c(
  "CHEMBL1232461",
  "CHEMBL1802814",
  "CHEMBL1802815",
  "CHEMBL1802973",
  "CHEMBL3752910",
  "CHEMBL5653589"
)

candidates_selected <- candidates_raw |>
  dplyr::filter(.data$molecule_chembl_id %in% expected_candidate_ids)

found_candidate_ids <- unique(candidates_selected$molecule_chembl_id)
missing_candidate_ids <- setdiff(expected_candidate_ids, found_candidate_ids)

if (length(missing_candidate_ids) > 0L) {
  stop(
    "A11 input is missing expected candidate(s): ",
    paste(missing_candidate_ids, collapse = ", "),
    call. = FALSE
  )
}

if (anyDuplicated(candidates_selected$molecule_chembl_id) > 0L) {
  duplicate_ids <- unique(
    candidates_selected$molecule_chembl_id[
      duplicated(candidates_selected$molecule_chembl_id)
    ]
  )

  stop(
    "A11 input contains duplicate candidate rows for: ",
    paste(duplicate_ids, collapse = ", "),
    call. = FALSE
  )
}

if (nrow(candidates_selected) != 6L) {
  stop(
    "Expected exactly 6 selected candidates, found ",
    nrow(candidates_selected),
    ".",
    call. = FALSE
  )
}

candidate_table <- candidates_selected |>
  dplyr::transmute(
    compound_id = clean_text(.data$molecule_chembl_id),
    compound_name = dplyr::if_else(
      is.na(clean_text(.data$display_name)),
      clean_text(.data$molecule_chembl_id),
      clean_text(.data$display_name)
    ),
    compound_class = "network_candidate",
    structure_source = "ChEMBL A11 candidate annotation",
    pubchem_cid = NA_character_,
    pubchem_title = NA_character_,
    iupac_name = NA_character_,
    canonical_smiles = clean_text(.data$canonical_smiles),
    isomeric_smiles = clean_text(.data$canonical_smiles),
    analysis_smiles = clean_text(.data$canonical_smiles),
    inchi = NA_character_,
    inchi_key = clean_text(.data$standard_inchi_key),
    molecular_formula = NA_character_,
    molecular_weight = NA_real_,
    structure_available = !is.na(clean_text(.data$canonical_smiles)),
    retrieval_status = dplyr::if_else(
      !is.na(clean_text(.data$canonical_smiles)),
      "available_from_A11",
      "missing_from_A11"
    ),
    retrieval_query = NA_character_
  )

# -----------------------------
# 5. Define reference ligands
# -----------------------------
reference_queries <- tibble::tribble(
  ~compound_id, ~compound_name, ~primary_query, ~alternative_query,
  "REF_PLITIDEPSIN", "Plitidepsin", "Plitidepsin", "Aplidin",
  "REF_DIDEMNIN_B", "Didemnin B", "Didemnin B", NA_character_,
  "REF_TERNATIN_4", "Ternatin-4", "Ternatin-4", "Ternatin 4",
  "REF_NANNOCYSTIN_A", "Nannocystin A", "Nannocystin A", NA_character_
)

# -----------------------------
# 6. Retrieve reference structures
# -----------------------------
reference_rows <- vector("list", nrow(reference_queries))

for (row_index in seq_len(nrow(reference_queries))) {
  query_names <- c(
    reference_queries$primary_query[row_index],
    reference_queries$alternative_query[row_index]
  )
  query_names <- query_names[!is.na(query_names) & trimws(query_names) != ""]

  retrieval_result <- NULL
  attempted_queries <- character(0)

  for (query_name in query_names) {
    attempted_queries <- c(attempted_queries, query_name)
    retrieval_result <- request_pubchem(query_name)

    # PubChem asks clients not to exceed five requests per second.
    Sys.sleep(0.25)

    if (isTRUE(retrieval_result$ok)) {
      break
    }
  }

  record <- if (
    !is.null(retrieval_result) && isTRUE(retrieval_result$ok)
  ) {
    retrieval_result$data
  } else {
    NULL
  }

  canonical_smiles <- get_property(
    record,
    c("CanonicalSMILES", "ConnectivitySMILES")
  )

  isomeric_smiles <- get_property(
    record,
    c("IsomericSMILES", "SMILES")
  )

  analysis_smiles <- if (!is.na(isomeric_smiles)) {
    isomeric_smiles
  } else {
    canonical_smiles
  }

  molecular_weight_text <- get_property(record, "MolecularWeight")
  molecular_weight <- suppressWarnings(as.numeric(molecular_weight_text))

  retrieval_succeeded <- (
    !is.null(retrieval_result) &&
    isTRUE(retrieval_result$ok) &&
    !is.na(analysis_smiles)
  )

  retrieval_status <- if (retrieval_succeeded) {
    "retrieved"
  } else if (is.null(retrieval_result)) {
    "request_not_attempted"
  } else {
    retrieval_result$message
  }

  retrieval_query <- if (
    !is.null(retrieval_result) && isTRUE(retrieval_result$ok)
  ) {
    retrieval_result$query
  } else {
    paste(attempted_queries, collapse = " | ")
  }

  reference_rows[[row_index]] <- tibble::tibble(
    compound_id = reference_queries$compound_id[row_index],
    compound_name = reference_queries$compound_name[row_index],
    compound_class = "reference_ligand",
    structure_source = "PubChem PUG REST",
    pubchem_cid = get_property(record, "CID"),
    pubchem_title = get_property(record, "Title"),
    iupac_name = get_property(record, "IUPACName"),
    canonical_smiles = canonical_smiles,
    isomeric_smiles = isomeric_smiles,
    analysis_smiles = analysis_smiles,
    inchi = get_property(record, "InChI"),
    inchi_key = get_property(record, "InChIKey"),
    molecular_formula = get_property(record, "MolecularFormula"),
    molecular_weight = molecular_weight,
    structure_available = !is.na(analysis_smiles),
    retrieval_status = retrieval_status,
    retrieval_query = retrieval_query
  )
}

reference_table <- dplyr::bind_rows(reference_rows)

if (nrow(reference_table) != 4L) {
  stop(
    "Internal error: reference table does not contain exactly 4 rows.",
    call. = FALSE
  )
}

# -----------------------------
# 7. Combine and validate structures
# -----------------------------
combined_table <- dplyr::bind_rows(candidate_table, reference_table) |>
  dplyr::mutate(
    structure_available = as.logical(.data$structure_available),
    analysis_smiles = clean_text(.data$analysis_smiles)
  )

if (nrow(combined_table) != 10L) {
  stop(
    "Expected 10 combined compounds, found ",
    nrow(combined_table),
    ".",
    call. = FALSE
  )
}

if (anyDuplicated(combined_table$compound_id) > 0L) {
  stop("Combined table contains duplicate compound_id values.", call. = FALSE)
}

missing_structures <- combined_table |>
  dplyr::filter(
    !.data$structure_available |
      is.na(.data$analysis_smiles)
  ) |>
  dplyr::pull(.data$compound_name)

summary_table <- tibble::tibble(
  metric = c(
    "candidate_count",
    "candidate_structures_available",
    "reference_structures_requested",
    "reference_structures_retrieved",
    "combined_compound_count",
    "combined_structures_available",
    "missing_structure_count",
    "missing_structures"
  ),
  value = c(
    as.character(nrow(candidate_table)),
    as.character(sum(candidate_table$structure_available)),
    as.character(nrow(reference_table)),
    as.character(sum(reference_table$structure_available)),
    as.character(nrow(combined_table)),
    as.character(sum(combined_table$structure_available)),
    as.character(length(missing_structures)),
    if (length(missing_structures) == 0L) {
      "None"
    } else {
      paste(missing_structures, collapse = "; ")
    }
  )
)

# -----------------------------
# 8. Write and verify outputs
# -----------------------------
safe_write_csv(reference_table, output_reference)
safe_write_csv(combined_table, output_combined)
safe_write_csv(summary_table, output_summary)

output_files <- c(
  output_reference,
  output_combined,
  output_summary
)

for (output_path in output_files) {
  if (
    !file.exists(output_path) ||
    is.na(file.info(output_path)$size) ||
    file.info(output_path)$size <= 0L
  ) {
    stop("Output verification failed: ", output_path, call. = FALSE)
  }
}

# Confirm written row counts by reading the outputs back from disk.
reference_check <- readr::read_csv(
  output_reference,
  show_col_types = FALSE,
  progress = FALSE
)
combined_check <- readr::read_csv(
  output_combined,
  show_col_types = FALSE,
  progress = FALSE
)
summary_check <- readr::read_csv(
  output_summary,
  show_col_types = FALSE,
  progress = FALSE
)

if (nrow(reference_check) != 4L) {
  stop("Written reference output does not contain 4 rows.", call. = FALSE)
}

if (nrow(combined_check) != 10L) {
  stop("Written combined output does not contain 10 rows.", call. = FALSE)
}

if (nrow(summary_check) != 8L) {
  stop("Written summary output does not contain 8 rows.", call. = FALSE)
}

# -----------------------------
# 9. Concise completion summary
# -----------------------------
cat("A16A structure retrieval completed.\n")
cat("Candidate count: ", nrow(candidate_table), "\n", sep = "")
cat(
  "Candidate structures available: ",
  sum(candidate_table$structure_available),
  "\n",
  sep = ""
)
cat("Reference structures requested: ", nrow(reference_table), "\n", sep = "")
cat(
  "Reference structures retrieved: ",
  sum(reference_table$structure_available),
  "\n",
  sep = ""
)
cat("Combined compound count: ", nrow(combined_table), "\n", sep = "")
cat(
  "Combined structures available: ",
  sum(combined_table$structure_available),
  "\n",
  sep = ""
)
cat(
  "Missing structures: ",
  if (length(missing_structures) == 0L) {
    "None"
  } else {
    paste(missing_structures, collapse = "; ")
  },
  "\n",
  sep = ""
)

if (length(missing_structures) > 0L) {
  stop(
    "A16A outputs were written, but one or more structures are missing. ",
    "Resolve the missing structures before running A16B.",
    call. = FALSE
  )
}
