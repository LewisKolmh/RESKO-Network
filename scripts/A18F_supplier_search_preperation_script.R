#!/usr/bin/env Rscript

# ============================================================
# RESKO A18F: Prepare external supplier-search manifests
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)
required <- c("readr")
missing <- required[!vapply(required, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing) > 0L) stop("Missing package(s): ", paste(missing, collapse = ", "), call. = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results <- file.path(root, "results")
input_ids <- file.path(results, "A18A_compound_identifiers.csv")
input_details <- file.path(results, "A18A_compound_detail_manifest.csv")
input_pubchem <- file.path(results, "A18B_compound_commercial_summary.csv")
outputs <- c(
  search_csv = file.path(results, "A18F_supplier_search_manifest.csv"),
  search_smi = file.path(results, "A18F_supplier_search_manifest.smi"),
  import_template = file.path(results, "A18F_supplier_results_template.csv"),
  status = file.path(results, "A18F_supplier_lookup_status.csv"),
  summary = file.path(results, "A18F_summary.csv")
)

for (p in c(input_ids, input_details)) {
  if (!file.exists(p) || file.info(p)$size <= 0L) stop("Missing input: ", p, call. = FALSE)
}

clean <- function(x) {
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

write_csv_safe <- function(x, p) {
  tmp <- paste0(p, ".tmp")
  readr::write_csv(as.data.frame(x, stringsAsFactors = FALSE), tmp, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0L) stop("Failed to write: ", p, call. = FALSE)
  if (file.exists(p)) unlink(p)
  if (!file.rename(tmp, p)) stop("Failed to install: ", p, call. = FALSE)
}

ids <- readr::read_csv(input_ids, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
details <- readr::read_csv(input_details, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
if (nrow(ids) != 16L || nrow(details) != 16L) stop("Expected 16 compounds in A18A inputs.", call. = FALSE)

m <- match(ids$compound_id, details$compound_id)
search <- data.frame(
  compound_id = ids$compound_id,
  compound_name = ids$compound_name,
  chembl_id = ids$chembl_id,
  pubchem_cid = ids$pubchem_cid,
  full_inchi_key = ids$full_inchi_key,
  parent_connectivity_key = ids$parent_connectivity_key,
  canonical_smiles = ids$canonical_smiles,
  compound_origin = details$compound_origin[m],
  biological_classification = details$biological_classification[m],
  progression_status = details$progression_status[m],
  stringsAsFactors = FALSE
)

if (file.exists(input_pubchem) && file.info(input_pubchem)$size > 0L) {
  pc <- readr::read_csv(input_pubchem, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
  pm <- match(search$compound_id, pc$compound_id)
  search$resolved_pubchem_cid <- pc$resolved_pubchem_cid[pm]
  search$pubchem_vendor_signal <- ifelse(
    pc$pubchem_vendor_lookup_status[pm] %in% c("vendor_information_retrieved", "supplier_names_retrieved"),
    "pubchem_vendor_category_or_supplier_signal",
    "no_programmatic_supplier_signal"
  )
} else {
  search$resolved_pubchem_cid <- search$pubchem_cid
  search$pubchem_vendor_signal <- "not_checked"
}

search$external_supplier_lookup_required <- TRUE
search$lookup_priority <- ifelse(
  grepl("reference|candidate", tolower(ifelse(is.na(search$biological_classification), "", search$biological_classification))),
  "primary",
  "comparison"
)

write_csv_safe(search, outputs[["search_csv"]])

valid_smi <- !is.na(clean(search$canonical_smiles))
smi_lines <- paste(search$canonical_smiles[valid_smi], search$compound_id[valid_smi], sep = "\t")
writeLines(smi_lines, outputs[["search_smi"]], useBytes = TRUE)

import_template <- data.frame(
  compound_id = search$compound_id,
  supplier_name = NA_character_,
  catalogue_number = NA_character_,
  product_name = NA_character_,
  product_url = NA_character_,
  listed_smiles = NA_character_,
  listed_inchi_key = NA_character_,
  form_description = NA_character_,
  purity = NA_character_,
  pack_size = NA_character_,
  price = NA_real_,
  currency = NA_character_,
  stock_status = NA_character_,
  lead_time = NA_character_,
  date_checked = NA_character_,
  source_platform = NA_character_,
  exact_identity_manually_verified = FALSE,
  reviewer = NA_character_,
  review_notes = NA_character_,
  stringsAsFactors = FALSE
)
write_csv_safe(import_template, outputs[["import_template"]])

status <- data.frame(
  compound_id = search$compound_id,
  compound_name = search$compound_name,
  pubchem_vendor_signal = search$pubchem_vendor_signal,
  supplier_names_recovered_programmatically = 0L,
  supplier_lookup_status = "external_supplier_search_required",
  exact_supplier_product_verified = FALSE,
  stringsAsFactors = FALSE
)
write_csv_safe(status, outputs[["status"]])

summary <- data.frame(
  metric = c("compound_count", "compounds_with_smiles_exported", "compounds_requiring_external_search", "programmatically_recovered_supplier_names", "verified_supplier_products"),
  value = c(nrow(search), sum(valid_smi), nrow(search), 0, 0),
  stringsAsFactors = FALSE
)
write_csv_safe(summary, outputs[["summary"]])

for (p in outputs) if (!file.exists(p) || file.info(p)$size <= 0L) stop("Output validation failed: ", p, call. = FALSE)
cat("\nA18F external supplier-search package completed.\n")
cat("Compounds exported: ", nrow(search), "\n", sep = "")
cat("SMILES records exported: ", sum(valid_smi), "\n", sep = "")
cat("Verified outputs: ", length(outputs), "\n", sep = "")
