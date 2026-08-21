#!/usr/bin/env Rscript

# ============================================================
# RESKO A18G: Import supplier results or complete with unresolved status
# ============================================================
#
# This script always completes when the A18F manifest exists.
# If A18G_supplier_results_completed.csv is absent or contains no populated
# supplier rows, it creates valid zero-record product outputs and marks all
# compounds as requiring external supplier lookup.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("readr", "jsonlite")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_packages) > 0L) {
  stop("Missing required package(s): ", paste(missing_packages, collapse = ", "), call. = FALSE)
}

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(root, "results")
manifest_file <- file.path(results_dir, "A18F_supplier_search_manifest.csv")
import_file <- file.path(results_dir, "A18G_supplier_results_completed.csv")

outputs <- c(
  products = file.path(results_dir, "A18G_verified_supplier_products.csv"),
  compound_summary = file.path(results_dir, "A18G_compound_supplier_summary.csv"),
  review = file.path(results_dir, "A18G_supplier_identity_review.csv"),
  json = file.path(results_dir, "A18G_verified_supplier_products.json"),
  summary = file.path(results_dir, "A18G_summary.csv")
)

if (!file.exists(manifest_file) || file.info(manifest_file)$size <= 0L) {
  stop("Missing A18F supplier-search manifest: ", manifest_file, call. = FALSE)
}

clean_text <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

safe_logical <- function(x) {
  y <- tolower(clean_text(x))
  ifelse(
    y %in% c("true", "t", "1", "yes"),
    TRUE,
    ifelse(y %in% c("false", "f", "0", "no"), FALSE, NA)
  )
}

collapse_unique <- function(x) {
  x <- sort(unique(clean_text(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else paste(x, collapse = "; ")
}

url_encode <- function(x) {
  vapply(
    as.character(x),
    function(value) utils::URLencode(value, reserved = TRUE),
    FUN.VALUE = character(1)
  )
}

safe_write_csv <- function(x, path) {
  temporary <- paste0(path, ".tmp")
  if (file.exists(temporary)) unlink(temporary)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (any(vapply(x, is.list, logical(1)))) {
    stop("Cannot write list column(s) to: ", path, call. = FALSE)
  }
  readr::write_csv(x, temporary, na = "")
  if (!file.exists(temporary) || file.info(temporary)$size <= 0L) {
    stop("Failed to write: ", path, call. = FALSE)
  }
  if (file.exists(path)) unlink(path)
  if (!file.rename(temporary, path)) stop("Failed to install: ", path, call. = FALSE)
}

empty_products <- function() {
  data.frame(
    product_record_id = character(),
    compound_id = character(),
    supplier_name = character(),
    catalogue_number = character(),
    product_name = character(),
    product_url = character(),
    listed_smiles = character(),
    listed_inchi_key = character(),
    form_description = character(),
    purity = character(),
    pack_size = character(),
    price = numeric(),
    currency = character(),
    stock_status = character(),
    lead_time = character(),
    date_checked = character(),
    source_platform = character(),
    exact_identity_manually_verified = logical(),
    reviewer = character(),
    review_notes = character(),
    manual_review_required = logical(),
    identity_status = character(),
    stringsAsFactors = FALSE
  )
}

required_import_columns <- c(
  "compound_id", "supplier_name", "catalogue_number", "product_name",
  "product_url", "listed_smiles", "listed_inchi_key", "form_description",
  "purity", "pack_size", "price", "currency", "stock_status", "lead_time",
  "date_checked", "source_platform", "exact_identity_manually_verified",
  "reviewer", "review_notes"
)

cat("Reading A18F supplier-search manifest...\n")
manifest <- readr::read_csv(
  manifest_file,
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE
)
if (nrow(manifest) != 16L) stop("A18F manifest must contain 16 compounds.", call. = FALSE)
if (anyDuplicated(manifest$compound_id) > 0L) stop("Duplicate compound IDs in A18F manifest.", call. = FALSE)

# Create the expected import file automatically when absent.
if (!file.exists(import_file)) {
  template <- data.frame(compound_id = manifest$compound_id, stringsAsFactors = FALSE)
  for (column_name in setdiff(required_import_columns, "compound_id")) {
    template[[column_name]] <- if (column_name == "exact_identity_manually_verified") FALSE else NA_character_
  }
  readr::write_csv(template, import_file, na = "")
  cat("Created empty A18G supplier-results file from the A18F manifest.\n")
}

cat("Reading available supplier-result rows...\n")
raw_results <- readr::read_csv(
  import_file,
  col_types = readr::cols(.default = readr::col_character()),
  progress = FALSE
)

missing_columns <- setdiff(required_import_columns, names(raw_results))
if (length(missing_columns) > 0L) {
  stop("Supplier-results file is missing column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}

populated <-
  !is.na(clean_text(raw_results$supplier_name)) |
  !is.na(clean_text(raw_results$product_url)) |
  !is.na(clean_text(raw_results$catalogue_number))

raw_results <- raw_results[populated, required_import_columns, drop = FALSE]

if (nrow(raw_results) > 0L) {
  unknown_ids <- setdiff(unique(clean_text(raw_results$compound_id)), manifest$compound_id)
  unknown_ids <- unknown_ids[!is.na(unknown_ids)]
  if (length(unknown_ids) > 0L) {
    stop("Unknown compound ID(s): ", paste(unknown_ids, collapse = ", "), call. = FALSE)
  }

  products <- raw_results
  products$exact_identity_manually_verified <- safe_logical(
    products$exact_identity_manually_verified
  )
  products$exact_identity_manually_verified[is.na(products$exact_identity_manually_verified)] <- FALSE
  products$manual_review_required <- !products$exact_identity_manually_verified
  products$identity_status <- ifelse(
    products$exact_identity_manually_verified,
    "exact_identity_manually_verified",
    "manual_identity_review_required"
  )
  products$product_record_id <- paste0(
    "SUPPLIER_PRODUCT_",
    sprintf("%06d", seq_len(nrow(products)))
  )
  products <- products[, c("product_record_id", required_import_columns,
                           "manual_review_required", "identity_status"), drop = FALSE]
} else {
  products <- empty_products()
}

cat("Building automated compound-level supplier status...\n")
summary_rows <- vector("list", nrow(manifest))

for (index in seq_len(nrow(manifest))) {
  compound_id <- manifest$compound_id[index]
  current <- products[products$compound_id == compound_id, , drop = FALSE]

  pubchem_cid <- clean_text(manifest$resolved_pubchem_cid[index])
  if (is.na(pubchem_cid)) pubchem_cid <- clean_text(manifest$pubchem_cid[index])

  pubchem_url <- if (
    length(pubchem_cid) > 0L && !is.na(pubchem_cid)
  ) {
    paste0("https://pubchem.ncbi.nlm.nih.gov/compound/", pubchem_cid,
           "#section=Chemical-Vendors")
  } else {
    NA_character_
  }

  molport_query <- if (!is.na(clean_text(manifest$full_inchi_key[index]))) {
    manifest$full_inchi_key[index]
  } else {
    manifest$compound_name[index]
  }

  summary_rows[[index]] <- data.frame(
    compound_id = compound_id,
    compound_name = manifest$compound_name[index],
    supplier_record_count = nrow(current),
    unique_supplier_count = length(unique(clean_text(current$supplier_name)[
      !is.na(clean_text(current$supplier_name))
    ])),
    supplier_names = collapse_unique(current$supplier_name),
    catalogue_numbers = collapse_unique(current$catalogue_number),
    product_urls = collapse_unique(current$product_url),
    verified_exact_product_count = sum(current$exact_identity_manually_verified, na.rm = TRUE),
    pubchem_vendor_page = pubchem_url,
    molport_search_url = paste0(
      "https://www.molport.com/shop/molecule-link/", url_encode(molport_query)
    ),
    emolecules_search_url = paste0(
      "https://www.emolecules.com/search?query=", url_encode(molport_query)
    ),
    supplier_lookup_status = if (
      nrow(current) == 0L
    ) {
      "external_supplier_search_required"
    } else if (any(current$exact_identity_manually_verified, na.rm = TRUE)) {
      "verified_supplier_product_available"
    } else {
      "supplier_records_require_identity_review"
    },
    exact_supplier_product_verified = any(
      current$exact_identity_manually_verified,
      na.rm = TRUE
    ),
    stringsAsFactors = FALSE
  )
}

compound_summary <- do.call(rbind, summary_rows)

if (nrow(products) > 0L) {
  review <- products[, c(
    "product_record_id", "compound_id", "supplier_name", "catalogue_number",
    "product_name", "listed_inchi_key", "form_description",
    "exact_identity_manually_verified", "manual_review_required", "reviewer",
    "review_notes"
  ), drop = FALSE]
} else {
  review <- data.frame(
    product_record_id = character(), compound_id = character(),
    supplier_name = character(), catalogue_number = character(),
    product_name = character(), listed_inchi_key = character(),
    form_description = character(), exact_identity_manually_verified = logical(),
    manual_review_required = logical(), reviewer = character(),
    review_notes = character(), stringsAsFactors = FALSE
  )
}

metrics <- data.frame(
  metric = c(
    "compound_count", "supplier_product_record_count",
    "compounds_with_supplier_results", "unique_supplier_count",
    "exact_products_manually_verified",
    "records_requiring_manual_identity_review",
    "compounds_requiring_external_supplier_search"
  ),
  value = c(
    nrow(manifest), nrow(products),
    sum(compound_summary$supplier_record_count > 0L),
    length(unique(clean_text(products$supplier_name)[
      !is.na(clean_text(products$supplier_name))
    ])),
    sum(products$exact_identity_manually_verified, na.rm = TRUE),
    sum(products$manual_review_required, na.rm = TRUE),
    sum(compound_summary$supplier_lookup_status ==
          "external_supplier_search_required")
  ),
  stringsAsFactors = FALSE
)

cat("Writing A18G outputs...\n")
safe_write_csv(products, outputs[["products"]])
safe_write_csv(compound_summary, outputs[["compound_summary"]])
safe_write_csv(review, outputs[["review"]])
jsonlite::write_json(
  products,
  outputs[["json"]],
  pretty = TRUE,
  auto_unbox = TRUE,
  dataframe = "rows",
  na = "null",
  null = "null"
)
safe_write_csv(metrics, outputs[["summary"]])

for (path in outputs) {
  if (!file.exists(path) || file.info(path)$size <= 0L) {
    stop("Output validation failed: ", path, call. = FALSE)
  }
}

summary_check <- readr::read_csv(outputs[["compound_summary"]], show_col_types = FALSE)
metrics_check <- readr::read_csv(outputs[["summary"]], show_col_types = FALSE)
if (nrow(summary_check) != 16L) stop("Compound supplier summary must contain 16 rows.", call. = FALSE)
if (nrow(metrics_check) != 7L) stop("A18G summary must contain seven rows.", call. = FALSE)

cat("\nA18G automated supplier-status import completed.\n")
cat("Supplier product records imported: ", nrow(products), "\n", sep = "")
cat("Compounds with supplier results: ", sum(compound_summary$supplier_record_count > 0L), "\n", sep = "")
cat("Compounds requiring external supplier search: ", sum(compound_summary$supplier_lookup_status == "external_supplier_search_required"), "\n", sep = "")
cat("Exact products manually verified: ", sum(products$exact_identity_manually_verified, na.rm = TRUE), "\n", sep = "")
cat("Verified outputs: ", length(outputs), "\n", sep = "")
