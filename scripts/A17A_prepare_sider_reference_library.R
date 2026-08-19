#!/usr/bin/env Rscript

# ============================================================
# RESKO A17A: Prepare the SIDER 4.1 reference library
# ============================================================
# Run from the RESKO project root:
#   Rscript scripts/A17A_prepare_sider_reference_library.R
#
# Purpose:
#   Download and freeze the official SIDER 4.1 files, validate their
#   formats, create clean drug, side-effect, frequency and indication
#   tables, and record source/checksum provenance.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("readr", "dplyr", "tibble", "httr2")
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
raw_dir <- file.path(results_dir, "A17A_sider_4.1_raw")
dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)

outputs <- c(
  drugs = file.path(results_dir, "A17A_sider_drugs.csv"),
  side_effects = file.path(results_dir, "A17A_sider_side_effects_pt.csv"),
  frequencies = file.path(results_dir, "A17A_sider_frequencies.csv"),
  indications = file.path(results_dir, "A17A_sider_indications_pt.csv"),
  manifest = file.path(results_dir, "A17A_sider_data_manifest.csv"),
  summary = file.path(results_dir, "A17A_sider_summary.csv")
)

sider_release <- "4.1"
sider_release_date <- "2015-10-21"
retrieval_time <- format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")
base_url <- "http://sideeffects.embl.de/media/download/"

source_files <- tibble::tribble(
  ~file_name, ~required,
  "README", TRUE,
  "drug_names.tsv", TRUE,
  "drug_atc.tsv", TRUE,
  "meddra_all_se.tsv.gz", TRUE,
  "meddra_freq.tsv.gz", TRUE,
  "meddra_all_indications.tsv.gz", TRUE,
  "meddra.tsv.gz", TRUE
) |>
  dplyr::mutate(
    source_url = paste0(base_url, .data$file_name),
    local_path = file.path(raw_dir, .data$file_name)
  )

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

download_with_retry <- function(url, destination, attempts = 5L) {
  if (file.exists(destination) && !is.na(file.info(destination)$size) && file.info(destination)$size > 0L) {
    return("existing_valid_file")
  }
  temp_path <- paste0(destination, ".download")
  if (file.exists(temp_path)) unlink(temp_path)
  last_error <- "No download attempt was completed"
  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_user_agent("RESKO-A17A/1.0 SIDER-4.1-retrieval") |>
        httr2::req_timeout(seconds = 120) |>
        httr2::req_perform(path = temp_path),
      error = function(e) e
    )
    if (!inherits(response, "error")) {
      status <- httr2::resp_status(response)
      if (status >= 200L && status < 300L && file.exists(temp_path) && file.info(temp_path)$size > 0L) {
        if (!file.rename(temp_path, destination)) {
          stop("Could not move downloaded file into place: ", destination, call. = FALSE)
        }
        return(paste0("downloaded_http_", status))
      }
      last_error <- paste0("HTTP ", status)
    } else {
      last_error <- conditionMessage(response)
    }
    if (file.exists(temp_path)) unlink(temp_path)
    if (attempt < attempts) Sys.sleep(min(2^(attempt - 1L), 16L))
  }
  stop("Failed to download ", url, " after ", attempts, " attempts. Last error: ", last_error, call. = FALSE)
}

cat("Downloading or validating SIDER 4.1 source files...\n")
download_status <- character(nrow(source_files))
for (i in seq_len(nrow(source_files))) {
  download_status[i] <- download_with_retry(source_files$source_url[i], source_files$local_path[i])
  cat("  ", source_files$file_name[i], ": ", download_status[i], "\n", sep = "")
}
source_files$download_status <- download_status

for (path in source_files$local_path) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Downloaded source file is missing or empty: ", path, call. = FALSE)
  }
}

read_no_header <- function(path, column_names, column_types = NULL) {
  tryCatch(
    readr::read_tsv(
      path,
      col_names = column_names,
      col_types = column_types,
      quote = "",
      comment = "",
      progress = FALSE,
      show_col_types = FALSE,
      trim_ws = FALSE
    ),
    error = function(e) stop("Failed to parse ", path, ": ", conditionMessage(e), call. = FALSE)
  )
}

drug_names_raw <- read_no_header(
  file.path(raw_dir, "drug_names.tsv"),
  c("stitch_flat_id", "drug_name"),
  readr::cols(.default = readr::col_character())
)

drug_atc_raw <- read_no_header(
  file.path(raw_dir, "drug_atc.tsv"),
  c("stitch_flat_id", "atc_code"),
  readr::cols(.default = readr::col_character())
)

side_effects_raw <- read_no_header(
  file.path(raw_dir, "meddra_all_se.tsv.gz"),
  c("stitch_flat_id", "stitch_stereo_id", "label_umls_id", "meddra_concept_type", "meddra_umls_id", "side_effect_name"),
  readr::cols(.default = readr::col_character())
)

frequencies_raw <- read_no_header(
  file.path(raw_dir, "meddra_freq.tsv.gz"),
  c(
    "stitch_flat_id", "stitch_stereo_id", "label_umls_id", "placebo",
    "frequency_description", "frequency_lower", "frequency_upper",
    "meddra_concept_type", "meddra_umls_id", "side_effect_name"
  ),
  readr::cols(.default = readr::col_character())
)

indications_raw <- read_no_header(
  file.path(raw_dir, "meddra_all_indications.tsv.gz"),
  c(
    "stitch_flat_id", "label_umls_id", "detection_method", "label_concept_name",
    "meddra_concept_type", "meddra_umls_id", "indication_name"
  ),
  readr::cols(.default = readr::col_character())
)

if (ncol(drug_names_raw) != 2L) stop("drug_names.tsv did not contain 2 columns.", call. = FALSE)
if (ncol(side_effects_raw) != 6L) stop("meddra_all_se.tsv.gz did not contain 6 columns.", call. = FALSE)
if (ncol(frequencies_raw) != 10L) stop("meddra_freq.tsv.gz did not contain 10 columns.", call. = FALSE)
if (ncol(indications_raw) != 7L) stop("meddra_all_indications.tsv.gz did not contain 7 columns.", call. = FALSE)

sider_drugs <- drug_names_raw |>
  dplyr::filter(!is.na(.data$stitch_flat_id), .data$stitch_flat_id != "") |>
  dplyr::mutate(
    drug_name = trimws(.data$drug_name),
    sider_version = sider_release
  ) |>
  dplyr::left_join(
    drug_atc_raw |>
      dplyr::filter(!is.na(.data$stitch_flat_id), .data$stitch_flat_id != "") |>
      dplyr::group_by(.data$stitch_flat_id) |>
      dplyr::summarise(atc_codes = paste(sort(unique(.data$atc_code)), collapse = "; "), .groups = "drop"),
    by = "stitch_flat_id"
  ) |>
  dplyr::distinct(.data$stitch_flat_id, .keep_all = TRUE) |>
  dplyr::arrange(.data$drug_name)

if (anyDuplicated(sider_drugs$stitch_flat_id) > 0L) stop("Duplicate SIDER drug identifiers remained after cleaning.", call. = FALSE)
if (nrow(sider_drugs) != 1430L) {
  stop("Expected 1,430 SIDER drugs, found ", nrow(sider_drugs), ". Check the source release and parsing.", call. = FALSE)
}

sider_side_effects <- side_effects_raw |>
  dplyr::filter(
    .data$meddra_concept_type == "PT",
    !is.na(.data$stitch_flat_id), .data$stitch_flat_id != "",
    !is.na(.data$meddra_umls_id), .data$meddra_umls_id != ""
  ) |>
  dplyr::mutate(sider_version = sider_release) |>
  dplyr::distinct(
    .data$stitch_flat_id, .data$stitch_stereo_id, .data$meddra_umls_id,
    .data$side_effect_name, .keep_all = TRUE
  ) |>
  dplyr::arrange(.data$stitch_flat_id, .data$side_effect_name)

sider_frequencies <- frequencies_raw |>
  dplyr::filter(
    .data$meddra_concept_type == "PT",
    !is.na(.data$stitch_flat_id), .data$stitch_flat_id != "",
    !is.na(.data$meddra_umls_id), .data$meddra_umls_id != ""
  ) |>
  dplyr::mutate(
    frequency_lower = suppressWarnings(as.numeric(.data$frequency_lower)),
    frequency_upper = suppressWarnings(as.numeric(.data$frequency_upper)),
    is_placebo = tolower(trimws(dplyr::coalesce(.data$placebo, ""))) == "placebo",
    sider_version = sider_release
  ) |>
  dplyr::distinct() |>
  dplyr::arrange(.data$stitch_flat_id, .data$side_effect_name)

sider_indications <- indications_raw |>
  dplyr::filter(
    .data$meddra_concept_type == "PT",
    !is.na(.data$stitch_flat_id), .data$stitch_flat_id != "",
    !is.na(.data$meddra_umls_id), .data$meddra_umls_id != ""
  ) |>
  dplyr::mutate(sider_version = sider_release) |>
  dplyr::distinct(
    .data$stitch_flat_id, .data$detection_method, .data$meddra_umls_id,
    .data$indication_name, .keep_all = TRUE
  ) |>
  dplyr::arrange(.data$stitch_flat_id, .data$indication_name)

unknown_side_effect_drugs <- setdiff(unique(sider_side_effects$stitch_flat_id), sider_drugs$stitch_flat_id)
if (length(unknown_side_effect_drugs) > 0L) {
  stop(
    "Side-effect data contain ", length(unknown_side_effect_drugs),
    " STITCH flat identifiers absent from drug_names.tsv. Review source identity before proceeding.",
    call. = FALSE
  )
}

file_sizes <- file.info(source_files$local_path)$size
md5_values <- unname(tools::md5sum(source_files$local_path))
manifest <- source_files |>
  dplyr::transmute(
    database = "SIDER",
    database_version = sider_release,
    release_date = sider_release_date,
    retrieval_timestamp = retrieval_time,
    file_name = .data$file_name,
    source_url = .data$source_url,
    local_path = sub(paste0("^", project_root, "/"), "", .data$local_path),
    file_size_bytes = as.numeric(file_sizes),
    md5 = md5_values,
    download_status = .data$download_status,
    license = "See official SIDER 4.1 download page; MedDRA-derived side-effect files CC BY-SA 4.0",
    source_note = "Official SIDER 4.1 download snapshot; database is historical and was released in 2015"
  )

summary_table <- tibble::tibble(
  metric = c(
    "sider_version", "release_date", "drug_count", "pt_side_effect_relationship_count",
    "unique_pt_side_effect_count", "frequency_record_count", "pt_indication_relationship_count",
    "source_file_count", "missing_source_file_count"
  ),
  value = c(
    sider_release,
    sider_release_date,
    as.character(nrow(sider_drugs)),
    as.character(nrow(sider_side_effects)),
    as.character(dplyr::n_distinct(sider_side_effects$meddra_umls_id)),
    as.character(nrow(sider_frequencies)),
    as.character(nrow(sider_indications)),
    as.character(nrow(source_files)),
    "0"
  )
)

safe_write_csv(sider_drugs, outputs[["drugs"]])
safe_write_csv(sider_side_effects, outputs[["side_effects"]])
safe_write_csv(sider_frequencies, outputs[["frequencies"]])
safe_write_csv(sider_indications, outputs[["indications"]])
safe_write_csv(manifest, outputs[["manifest"]])
safe_write_csv(summary_table, outputs[["summary"]])

expected_rows <- c(
  drugs = nrow(sider_drugs),
  side_effects = nrow(sider_side_effects),
  frequencies = nrow(sider_frequencies),
  indications = nrow(sider_indications),
  manifest = nrow(manifest),
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

cat("A17A SIDER 4.1 reference library completed.\n")
cat("SIDER drugs: ", nrow(sider_drugs), "\n", sep = "")
cat("PT drug-side-effect relationships: ", nrow(sider_side_effects), "\n", sep = "")
cat("Unique PT side effects: ", dplyr::n_distinct(sider_side_effects$meddra_umls_id), "\n", sep = "")
cat("PT frequency records: ", nrow(sider_frequencies), "\n", sep = "")
cat("PT indication relationships: ", nrow(sider_indications), "\n", sep = "")
cat("Source files frozen: ", nrow(source_files), "\n", sep = "")
cat("Verified outputs: ", length(outputs), "\n", sep = "")
