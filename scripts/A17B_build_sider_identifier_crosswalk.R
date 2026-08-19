#!/usr/bin/env Rscript

# ============================================================
# RESKO A17B V2: SIDER identifier harmonisation and current mapping
# ============================================================
# Run from the RESKO project root:
#   Rscript scripts/A17B_build_sider_identifier_crosswalk.R
#
# Inputs:
#   results/A17A_sider_drugs.csv
#   results/A16_structures_combined.csv
#
# Outputs:
#   results/A17B_sider_identifier_crosswalk.csv
#   results/A17B_current_compound_sider_mapping.csv
#   results/A17B_ambiguous_matches_review.csv
#   results/A17B_unmapped_compounds.csv
#   results/A17B_mapping_summary.csv
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("readr", "dplyr", "tibble", "httr2", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_packages, collapse = ", "),
    ". Install with: Rscript -e 'install.packages(c(",
    paste(sprintf("\"%s\"", missing_packages), collapse = ", "),
    "), repos=\"https://cloud.r-project.org\")'",
    call. = FALSE
  )
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(project_root, "scripts")) || !dir.exists(file.path(project_root, "results"))) {
  stop("Run this script from the RESKO project root containing scripts/ and results/.", call. = FALSE)
}

results_dir <- file.path(project_root, "results")
input_sider <- file.path(results_dir, "A17A_sider_drugs.csv")
input_current <- file.path(results_dir, "A16_structures_combined.csv")
cache_file <- file.path(results_dir, "A17B_pubchem_property_cache.csv")

outputs <- c(
  crosswalk = file.path(results_dir, "A17B_sider_identifier_crosswalk.csv"),
  current_mapping = file.path(results_dir, "A17B_current_compound_sider_mapping.csv"),
  ambiguous = file.path(results_dir, "A17B_ambiguous_matches_review.csv"),
  unmapped = file.path(results_dir, "A17B_unmapped_compounds.csv"),
  summary = file.path(results_dir, "A17B_mapping_summary.csv")
)

for (path in c(input_sider, input_current)) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Required input is missing or empty: ", path, call. = FALSE)
  }
}

clean_text <- function(x) {
  output <- as.character(x)
  output[is.na(output) | trimws(output) == ""] <- NA_character_
  output
}

normalise_name <- function(x) {
  x <- tolower(clean_text(x))
  x <- gsub("[^a-z0-9]+", "", x)
  x[x == ""] <- NA_character_
  x
}

connectivity_key <- function(inchi_key) {
  key <- clean_text(inchi_key)
  ifelse(!is.na(key) & nchar(key) >= 14L, substr(key, 1L, 14L), NA_character_)
}

stitch_to_pubchem_cid <- function(stitch_id) {
  stitch_id <- clean_text(stitch_id)

  # SIDER 4.1 uses STITCH flat identifiers such as CID100002524.
  # The first digit after CID is the STITCH flat/stereo prefix; the
  # remaining digits encode the PubChem CID, including leading zeros.
  valid <- !is.na(stitch_id) & grepl("^CID[01][0-9]+$", stitch_id)
  numeric_part <- ifelse(valid, sub("^CID[01]", "", stitch_id), NA_character_)
  output <- suppressWarnings(as.numeric(numeric_part))
  output[!valid] <- NA_real_
  output
}

backup_existing <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  backup_path <- paste0(stem, "_previous_", stamp, if (extension == "") "" else paste0(".", extension))
  if (!file.copy(path, backup_path, overwrite = FALSE)) {
    stop("Could not preserve previous output: ", path, call. = FALSE)
  }
  if (!file.exists(backup_path) || file.info(backup_path)$size <= 0L) {
    stop("Backup verification failed: ", backup_path, call. = FALSE)
  }
  invisible(backup_path)
}

safe_write_csv <- function(data, path) {
  temp_path <- paste0(path, ".tmp")
  if (file.exists(temp_path)) unlink(temp_path)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  list_columns <- names(data)[vapply(data, is.list, logical(1))]
  if (length(list_columns) > 0L) {
    stop("Refusing to write list-column(s): ", paste(list_columns, collapse = ", "), call. = FALSE)
  }
  readr::write_csv(data, temp_path, na = "")
  if (!file.exists(temp_path) || is.na(file.info(temp_path)$size) || file.info(temp_path)$size <= 0L) {
    stop("Failed to create temporary output: ", temp_path, call. = FALSE)
  }
  backup_existing(path)
  if (!file.rename(temp_path, path)) stop("Could not move output into place: ", path, call. = FALSE)
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Output verification failed: ", path, call. = FALSE)
  }
  invisible(path)
}

request_pubchem_batch <- function(cids, attempts = 5L) {
  cid_string <- paste(format(cids, scientific = FALSE, trim = TRUE), collapse = ",")
  properties <- paste(
    c("Title", "IUPACName", "CanonicalSMILES", "IsomericSMILES", "InChI", "InChIKey", "MolecularFormula", "MolecularWeight"),
    collapse = ","
  )
  url <- paste0(
    "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/", cid_string,
    "/property/", properties, "/JSON"
  )
  last_error <- "No request completed"
  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_user_agent("RESKO-A17B/2.0 SIDER-identifier-harmonisation") |>
        httr2::req_timeout(seconds = 120) |>
        httr2::req_perform(),
      error = function(e) e
    )
    if (!inherits(response, "error")) {
      status <- httr2::resp_status(response)
      if (status >= 200L && status < 300L) {
        parsed <- tryCatch(
          jsonlite::fromJSON(httr2::resp_body_string(response), simplifyDataFrame = TRUE),
          error = function(e) e
        )
        if (!inherits(parsed, "error") && !is.null(parsed$PropertyTable$Properties)) {
          properties_table <- parsed$PropertyTable$Properties
          if (is.data.frame(properties_table) && nrow(properties_table) > 0L) {
            return(tibble::as_tibble(properties_table))
          }
        }
        last_error <- "PubChem returned no usable property table"
      } else {
        last_error <- paste0("HTTP ", status)
      }
    } else {
      last_error <- conditionMessage(response)
    }
    if (attempt < attempts) Sys.sleep(min(2^(attempt - 1L), 16L))
  }
  stop("PubChem batch failed after ", attempts, " attempts: ", last_error, call. = FALSE)
}

sider_drugs <- readr::read_csv(input_sider, show_col_types = FALSE, progress = FALSE)
current_compounds <- readr::read_csv(input_current, show_col_types = FALSE, progress = FALSE)

required_sider_columns <- c("stitch_flat_id", "drug_name", "atc_codes", "sider_version")
required_current_columns <- c(
  "compound_id", "compound_name", "compound_class", "pubchem_cid",
  "canonical_smiles", "analysis_smiles", "inchi_key"
)
missing_sider <- setdiff(required_sider_columns, names(sider_drugs))
missing_current <- setdiff(required_current_columns, names(current_compounds))
if (length(missing_sider) > 0L) stop("SIDER input missing: ", paste(missing_sider, collapse = ", "), call. = FALSE)
if (length(missing_current) > 0L) stop("Current-compound input missing: ", paste(missing_current, collapse = ", "), call. = FALSE)
if (nrow(sider_drugs) != 1430L) stop("Expected 1,430 SIDER drugs, found ", nrow(sider_drugs), ".", call. = FALSE)
if (nrow(current_compounds) != 10L) stop("Expected 10 current compounds, found ", nrow(current_compounds), ".", call. = FALSE)
if (anyDuplicated(sider_drugs$stitch_flat_id) > 0L) stop("Duplicate SIDER flat identifiers detected.", call. = FALSE)
if (anyDuplicated(current_compounds$compound_id) > 0L) stop("Duplicate current compound identifiers detected.", call. = FALSE)

sider_drugs <- sider_drugs |>
  dplyr::mutate(
    pubchem_cid = stitch_to_pubchem_cid(.data$stitch_flat_id),
    normalised_drug_name = normalise_name(.data$drug_name)
  )
if (any(is.na(sider_drugs$pubchem_cid))) {
  bad_ids <- sider_drugs$stitch_flat_id[is.na(sider_drugs$pubchem_cid)]
  stop("Could not derive PubChem CID from STITCH ID(s): ", paste(head(bad_ids, 10L), collapse = ", "), call. = FALSE)
}
if (anyDuplicated(sider_drugs$pubchem_cid) > 0L) {
  stop("Multiple SIDER flat identifiers resolved to the same PubChem CID.", call. = FALSE)
}

required_cids <- sort(unique(sider_drugs$pubchem_cid))
cache <- if (file.exists(cache_file) && file.info(cache_file)$size > 0L) {
  readr::read_csv(cache_file, show_col_types = FALSE, progress = FALSE)
} else {
  tibble::tibble()
}
if (nrow(cache) > 0L && !"pubchem_cid" %in% names(cache)) {
  stop("Existing PubChem cache lacks pubchem_cid: ", cache_file, call. = FALSE)
}

cached_cids <- if (nrow(cache) > 0L) suppressWarnings(as.numeric(cache$pubchem_cid)) else numeric(0)
missing_cids <- setdiff(required_cids, cached_cids)

if (length(missing_cids) > 0L) {
  cat("Retrieving PubChem properties for ", length(missing_cids), " SIDER compounds...\n", sep = "")
  batches <- split(missing_cids, ceiling(seq_along(missing_cids) / 75L))
  retrieved <- vector("list", length(batches))
  for (i in seq_along(batches)) {
    cat("  PubChem batch ", i, " of ", length(batches), "\n", sep = "")
    retrieved[[i]] <- request_pubchem_batch(batches[[i]])
    Sys.sleep(0.25)
  }
  retrieved_table <- dplyr::bind_rows(retrieved)
  optional_pubchem_columns <- c(
    "Title", "IUPACName", "CanonicalSMILES", "ConnectivitySMILES",
    "IsomericSMILES", "SMILES", "InChI", "InChIKey",
    "MolecularFormula", "MolecularWeight"
  )
  for (column_name in optional_pubchem_columns) {
    if (!column_name %in% names(retrieved_table)) {
      retrieved_table[[column_name]] <- NA_character_
    }
  }
  new_cache <- retrieved_table |>
    dplyr::transmute(
      pubchem_cid = suppressWarnings(as.numeric(.data$CID)),
      pubchem_title = clean_text(.data$Title),
      iupac_name = clean_text(.data$IUPACName),
      canonical_smiles = dplyr::coalesce(
        clean_text(.data$CanonicalSMILES),
        clean_text(.data$ConnectivitySMILES)
      ),
      isomeric_smiles = dplyr::coalesce(
        clean_text(.data$IsomericSMILES),
        clean_text(.data$SMILES)
      ),
      inchi = clean_text(.data$InChI),
      inchi_key = clean_text(.data$InChIKey),
      molecular_formula = clean_text(.data$MolecularFormula),
      molecular_weight = suppressWarnings(as.numeric(.data$MolecularWeight)),
      retrieval_timestamp = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
    )
  cache <- dplyr::bind_rows(cache, new_cache) |>
    dplyr::distinct(.data$pubchem_cid, .keep_all = TRUE) |>
    dplyr::arrange(.data$pubchem_cid)
  safe_write_csv(cache, cache_file)
}

pubchem_properties <- cache |>
  dplyr::filter(.data$pubchem_cid %in% required_cids) |>
  dplyr::distinct(.data$pubchem_cid, .keep_all = TRUE)

missing_property_cids <- setdiff(required_cids, pubchem_properties$pubchem_cid)
if (length(missing_property_cids) > 0L) {
  stop(
    "PubChem returned no property record for ", length(missing_property_cids),
    " required CID(s). First missing CID(s): ", paste(head(missing_property_cids, 10L), collapse = ", "),
    call. = FALSE
  )
}

sider_crosswalk <- sider_drugs |>
  dplyr::left_join(pubchem_properties, by = "pubchem_cid") |>
  dplyr::mutate(
    connectivity_inchi_key = connectivity_key(.data$inchi_key),
    structure_available = !is.na(clean_text(.data$canonical_smiles)),
    pubchem_mapping_status = dplyr::if_else(
      !is.na(.data$inchi_key) & .data$structure_available,
      "pubchem_identifier_and_structure_retrieved",
      "pubchem_record_incomplete"
    ),
    identifier_provenance = "SIDER STITCH flat identifier converted to PubChem CID; properties retrieved from PubChem PUG REST"
  ) |>
  dplyr::arrange(.data$drug_name)

if (nrow(sider_crosswalk) != 1430L) stop("Crosswalk does not contain 1,430 rows.", call. = FALSE)
if (anyDuplicated(sider_crosswalk$stitch_flat_id) > 0L) stop("Crosswalk contains duplicate SIDER IDs.", call. = FALSE)

current_prepared <- current_compounds |>
  dplyr::mutate(
    current_pubchem_cid = suppressWarnings(as.numeric(.data$pubchem_cid)),
    current_inchi_key = clean_text(.data$inchi_key),
    current_connectivity_key = connectivity_key(.data$current_inchi_key),
    normalised_compound_name = normalise_name(.data$compound_name)
  )

mapping_rows <- vector("list", nrow(current_prepared))
for (i in seq_len(nrow(current_prepared))) {
  compound <- current_prepared[i, , drop = FALSE]

  cid_matches <- if (!is.na(compound$current_pubchem_cid[[1]])) {
    sider_crosswalk |>
      dplyr::filter(.data$pubchem_cid == compound$current_pubchem_cid[[1]])
  } else sider_crosswalk[0, , drop = FALSE]

  full_key_matches <- if (!is.na(compound$current_inchi_key[[1]])) {
    sider_crosswalk |>
      dplyr::filter(.data$inchi_key == compound$current_inchi_key[[1]])
  } else sider_crosswalk[0, , drop = FALSE]

  connectivity_matches <- if (!is.na(compound$current_connectivity_key[[1]])) {
    sider_crosswalk |>
      dplyr::filter(.data$connectivity_inchi_key == compound$current_connectivity_key[[1]])
  } else sider_crosswalk[0, , drop = FALSE]

  name_matches <- if (!is.na(compound$normalised_compound_name[[1]])) {
    sider_crosswalk |>
      dplyr::filter(
        .data$normalised_drug_name == compound$normalised_compound_name[[1]] |
          normalise_name(.data$pubchem_title) == compound$normalised_compound_name[[1]]
      )
  } else sider_crosswalk[0, , drop = FALSE]

  if (nrow(cid_matches) > 0L) {
    selected <- cid_matches
    match_type <- "exact_pubchem_cid_match"
    confidence <- "high"
  } else if (nrow(full_key_matches) > 0L) {
    selected <- full_key_matches
    match_type <- "exact_inchi_key_match"
    confidence <- "high"
  } else if (nrow(connectivity_matches) > 0L) {
    selected <- connectivity_matches
    match_type <- "connectivity_inchi_key_match_requires_form_review"
    confidence <- "moderate"
  } else if (nrow(name_matches) > 0L) {
    selected <- name_matches
    match_type <- "normalised_name_match_requires_review"
    confidence <- "low"
  } else {
    selected <- sider_crosswalk[NA_integer_, , drop = FALSE]
    selected[1, ] <- NA
    match_type <- "not_represented_in_sider_4.1"
    confidence <- "none"
  }

  match_count <- if (match_type == "not_represented_in_sider_4.1") 0L else nrow(selected)
  mapping_rows[[i]] <- tibble::tibble(
    compound_id = compound$compound_id[[1]],
    compound_name = compound$compound_name[[1]],
    compound_class = compound$compound_class[[1]],
    current_pubchem_cid = compound$current_pubchem_cid[[1]],
    current_inchi_key = compound$current_inchi_key[[1]],
    current_connectivity_key = compound$current_connectivity_key[[1]],
    sider_match_type = match_type,
    sider_match_confidence = confidence,
    sider_match_count = match_count,
    sider_stitch_flat_id = if (match_count > 0L) selected$stitch_flat_id else NA_character_,
    sider_drug_name = if (match_count > 0L) selected$drug_name else NA_character_,
    sider_pubchem_cid = if (match_count > 0L) selected$pubchem_cid else NA_real_,
    sider_inchi_key = if (match_count > 0L) selected$inchi_key else NA_character_,
    sider_connectivity_key = if (match_count > 0L) selected$connectivity_inchi_key else NA_character_,
    manual_review_required = match_count > 1L || confidence %in% c("moderate", "low"),
    interpretation = dplyr::case_when(
      match_count == 0L ~ "Not represented in SIDER 4.1; this does not imply absence of side effects",
      confidence == "high" ~ "High-confidence identifier match to a SIDER 4.1 drug",
      confidence == "moderate" ~ "Connectivity match; review salt, stereochemistry, parent and component form",
      TRUE ~ "Name-only match; manual identity review required"
    )
  )
}

current_mapping <- dplyr::bind_rows(mapping_rows) |>
  dplyr::arrange(.data$compound_class, .data$compound_name, .data$sider_stitch_flat_id)

ambiguous_review <- current_mapping |>
  dplyr::filter(.data$manual_review_required | .data$sider_match_count > 1L)

unmapped_compounds <- current_mapping |>
  dplyr::filter(.data$sider_match_type == "not_represented_in_sider_4.1")

compound_level <- current_mapping |>
  dplyr::group_by(.data$compound_id, .data$compound_name) |>
  dplyr::summarise(
    match_type = dplyr::first(.data$sider_match_type),
    match_confidence = dplyr::first(.data$sider_match_confidence),
    match_count = max(.data$sider_match_count),
    manual_review_required = any(.data$manual_review_required),
    .groups = "drop"
  )

summary_table <- tibble::tibble(
  metric = c(
    "sider_drug_count",
    "sider_pubchem_records_retrieved",
    "sider_structures_available",
    "current_compound_count",
    "current_compounds_high_confidence_mapped",
    "current_compounds_connectivity_mapped",
    "current_compounds_name_only_mapped",
    "current_compounds_ambiguous_or_review_required",
    "current_compounds_not_represented",
    "mapping_rows"
  ),
  value = c(
    as.character(nrow(sider_crosswalk)),
    as.character(sum(!is.na(sider_crosswalk$inchi_key))),
    as.character(sum(sider_crosswalk$structure_available)),
    as.character(nrow(current_compounds)),
    as.character(sum(compound_level$match_confidence == "high")),
    as.character(sum(compound_level$match_confidence == "moderate")),
    as.character(sum(compound_level$match_confidence == "low")),
    as.character(sum(compound_level$manual_review_required)),
    as.character(sum(compound_level$match_confidence == "none")),
    as.character(nrow(current_mapping))
  )
)

safe_write_csv(sider_crosswalk, outputs[["crosswalk"]])
safe_write_csv(current_mapping, outputs[["current_mapping"]])
safe_write_csv(ambiguous_review, outputs[["ambiguous"]])
safe_write_csv(unmapped_compounds, outputs[["unmapped"]])
safe_write_csv(summary_table, outputs[["summary"]])

expected_rows <- c(
  crosswalk = nrow(sider_crosswalk),
  current_mapping = nrow(current_mapping),
  ambiguous = nrow(ambiguous_review),
  unmapped = nrow(unmapped_compounds),
  summary = nrow(summary_table)
)

for (name in names(outputs)) {
  path <- outputs[[name]]
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Output verification failed: ", path, call. = FALSE)
  }
  check <- readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
  if (nrow(check) != expected_rows[[name]]) {
    stop("Written row-count mismatch for ", path, ".", call. = FALSE)
  }
}

cat("A17B SIDER identifier harmonisation completed.\n")
cat("SIDER crosswalk rows: ", nrow(sider_crosswalk), "\n", sep = "")
cat("SIDER structures available: ", sum(sider_crosswalk$structure_available), "\n", sep = "")
cat("Current compounds assessed: ", nrow(current_compounds), "\n", sep = "")
cat("High-confidence SIDER mappings: ", sum(compound_level$match_confidence == "high"), "\n", sep = "")
cat("Connectivity mappings requiring review: ", sum(compound_level$match_confidence == "moderate"), "\n", sep = "")
cat("Name-only mappings requiring review: ", sum(compound_level$match_confidence == "low"), "\n", sep = "")
cat("Current compounds not represented in SIDER 4.1: ", sum(compound_level$match_confidence == "none"), "\n", sep = "")
cat("Ambiguous or review-required compounds: ", sum(compound_level$manual_review_required), "\n", sep = "")
cat("Verified outputs: ", length(outputs), "\n", sep = "")
