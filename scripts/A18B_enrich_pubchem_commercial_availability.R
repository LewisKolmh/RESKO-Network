#!/usr/bin/env Rscript

# RESKO A18B: PubChem supplier-name enrichment
options(stringsAsFactors = FALSE, warn = 1)

required_packages <- c("readr", "httr2", "jsonlite")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
if (length(missing_packages) > 0L) stop("Missing packages: ", paste(missing_packages, collapse = ", "), call. = FALSE)

root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(root, "results")
cache_dir <- file.path(results_dir, "A18B_api_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

input_file <- file.path(results_dir, "A18A_supplier_lookup_manifest.csv")
outputs <- c(
  summary = file.path(results_dir, "A18B_compound_commercial_summary.csv"),
  products = file.path(results_dir, "A18B_commercial_products.csv"),
  suppliers = file.path(results_dir, "A18B_supplier_directory.csv"),
  review = file.path(results_dir, "A18B_identity_review.csv"),
  metrics = file.path(results_dir, "A18B_summary.csv"),
  report = file.path(results_dir, "A18B_commercial_availability_report.html")
)

if (!file.exists(input_file) || file.info(input_file)$size <= 0L) stop("Missing A18A supplier manifest.", call. = FALSE)

clean_text <- function(x) {
  if (is.null(x) || length(x) == 0L) return(NA_character_)
  x <- as.character(x)
  x[is.na(x) | trimws(x) == ""] <- NA_character_
  x
}

scalar_text <- function(x, default = NA_character_) {
  x <- clean_text(x)
  x <- x[!is.na(x)]
  if (length(x) == 0L) default else x[1]
}

safe_numeric <- function(x) suppressWarnings(as.numeric(x))

collapse_unique <- function(x) {
  x <- sort(unique(clean_text(x)))
  x <- x[!is.na(x)]
  if (length(x) == 0L) NA_character_ else paste(x, collapse = "; ")
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

backup_file <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  ext <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  backup <- paste0(stem, "_previous_", stamp, if (ext == "") "" else paste0(".", ext))
  if (!file.copy(path, backup, overwrite = FALSE)) stop("Could not back up ", path, call. = FALSE)
  invisible(backup)
}

safe_write_csv <- function(x, path) {
  tmp <- paste0(path, ".tmp")
  if (file.exists(tmp)) unlink(tmp)
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (any(vapply(x, is.list, logical(1)))) stop("List column in output: ", path, call. = FALSE)
  readr::write_csv(x, tmp, na = "")
  if (!file.exists(tmp) || file.info(tmp)$size <= 0L) stop("Failed to write ", path, call. = FALSE)
  backup_file(path)
  if (!file.rename(tmp, path)) stop("Failed to install ", path, call. = FALSE)
}

request_json <- function(url, cache_name, attempts = 4L) {
  cache_path <- file.path(cache_dir, cache_name)
  if (file.exists(cache_path) && file.info(cache_path)$size > 0L) {
    cached <- tryCatch(jsonlite::fromJSON(cache_path, simplifyVector = FALSE), error = function(e) NULL)
    if (!is.null(cached)) return(cached)
    unlink(cache_path)
  }
  for (attempt in seq_len(attempts)) {
    response <- tryCatch(
      httr2::request(url) |>
        httr2::req_user_agent("RESKO-A18B-suppliers/1.0") |>
        httr2::req_timeout(seconds = 90) |>
        httr2::req_perform(),
      error = function(e) NULL
    )
    if (!is.null(response)) {
      code <- httr2::resp_status(response)
      if (code >= 200L && code < 300L) {
        body <- httr2::resp_body_string(response)
        parsed <- tryCatch(jsonlite::fromJSON(body, simplifyVector = FALSE), error = function(e) NULL)
        if (!is.null(parsed)) {
          writeLines(body, cache_path, useBytes = TRUE)
          Sys.sleep(0.25)
          return(parsed)
        }
      }
      if (code == 404L) return(NULL)
    }
    if (attempt < attempts) Sys.sleep(min(2^attempt, 8))
  }
  NULL
}

resolve_cid <- function(existing_cid, inchikey, compound_id) {
  cid <- safe_numeric(existing_cid)
  if (length(cid) > 0L && !is.na(cid[1]) && cid[1] > 0) return(cid[1])
  key <- scalar_text(inchikey)
  if (is.na(key)) return(NA_real_)
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/inchikey/", utils::URLencode(key, reserved = TRUE), "/cids/JSON")
  obj <- request_json(url, paste0("cid_", gsub("[^A-Za-z0-9]+", "_", compound_id), ".json"))
  if (is.null(obj) || is.null(obj$IdentifierList$CID)) return(NA_real_)
  values <- safe_numeric(unlist(obj$IdentifierList$CID, recursive = TRUE, use.names = FALSE))
  values <- values[!is.na(values)]
  if (length(values) == 0L) NA_real_ else values[1]
}

empty_products <- function() data.frame(
  commercial_product_id = character(), compound_id = character(), compound_name = character(),
  pubchem_cid = numeric(), supplier_name = character(), supplier_source_id = character(),
  product_name = character(), catalogue_number = character(), product_url = character(),
  commercial_data_source = character(), date_checked = character(), identity_match_class = character(),
  manual_verification_status = character(), stock_status = character(), purity = character(),
  pack_size = character(), price = numeric(), currency = character(), lead_time = character(),
  stringsAsFactors = FALSE
)

extract_vendor_references <- function(obj, compound_id, compound_name, cid) {
  refs <- obj$Record$Reference
  if (is.null(refs) || length(refs) == 0L) return(empty_products())
  rows <- list()
  j <- 1L
  for (ref in refs) {
    category <- scalar_text(ref$Category)
    if (is.na(category) || tolower(category) != "chemical vendors") next
    supplier <- scalar_text(ref$SourceName, "Supplier name not exposed in the public PubChem record")
    source_id <- scalar_text(ref$SourceID)
    url <- scalar_text(ref$URL)
    if (is.na(url)) url <- scalar_text(ref$LicenseURL)
    description <- scalar_text(ref$Description)
    rows[[j]] <- data.frame(
      commercial_product_id = NA_character_, compound_id = compound_id, compound_name = compound_name,
      pubchem_cid = cid, supplier_name = supplier, supplier_source_id = source_id,
      product_name = description, catalogue_number = source_id, product_url = url,
      commercial_data_source = "PubChem Chemical Vendors", date_checked = as.character(Sys.Date()),
      identity_match_class = "pubchem_supplier_association_requires_product_review",
      manual_verification_status = "not_verified", stock_status = "not_independently_verified",
      purity = NA_character_, pack_size = NA_character_, price = NA_real_, currency = NA_character_,
      lead_time = NA_character_, stringsAsFactors = FALSE
    )
    j <- j + 1L
  }
  if (length(rows) == 0L) return(empty_products())
  out <- do.call(rbind, rows)
  out <- unique(out)
  out$commercial_product_id <- paste0("PUBCHEM_VENDOR_", sprintf("%06d", seq_len(nrow(out))))
  out
}

make_html_table <- function(x) {
  x <- as.data.frame(x, stringsAsFactors = FALSE)
  if (nrow(x) == 0L) return("<p>No records available.</p>")
  for (nm in names(x)) x[[nm]] <- ifelse(is.na(x[[nm]]), "", html_escape(x[[nm]]))
  header <- paste0("<tr>", paste0("<th>", html_escape(names(x)), "</th>", collapse = ""), "</tr>")
  rows <- vapply(seq_len(nrow(x)), function(i) paste0("<tr>", paste0("<td>", unlist(x[i, , drop = FALSE]), "</td>", collapse = ""), "</tr>"), character(1))
  paste0("<table><thead>", header, "</thead><tbody>", paste(rows, collapse = "\n"), "</tbody></table>")
}

cat("Reading A18A supplier manifest...\n")
manifest <- readr::read_csv(input_file, col_types = readr::cols(.default = readr::col_character()), progress = FALSE)
if (nrow(manifest) != 16L) stop("Expected 16 compounds.", call. = FALSE)

resolved <- rep(NA_real_, nrow(manifest))
product_tables <- vector("list", nrow(manifest))
lookup_status <- rep("not_started", nrow(manifest))

for (i in seq_len(nrow(manifest))) {
  cat("  Supplier lookup ", i, " of ", nrow(manifest), ": ", manifest$compound_name[i], "\n", sep = "")
  resolved[i] <- resolve_cid(manifest$pubchem_cid[i], manifest$full_inchi_key[i], manifest$compound_id[i])
  if (is.na(resolved[i])) {
    lookup_status[i] <- "pubchem_cid_not_resolved"
    next
  }
  url <- paste0("https://pubchem.ncbi.nlm.nih.gov/rest/pug_view/data/compound/", as.integer(resolved[i]), "/JSON")
  obj <- request_json(url, paste0("full_record_", as.integer(resolved[i]), ".json"))
  if (is.null(obj) || is.null(obj$Record)) {
    lookup_status[i] <- "pubchem_record_unavailable"
    next
  }
  current <- extract_vendor_references(obj, manifest$compound_id[i], manifest$compound_name[i], resolved[i])
  if (nrow(current) == 0L) {
    lookup_status[i] <- "no_supplier_references_found"
  } else {
    lookup_status[i] <- "supplier_names_retrieved"
    product_tables[[i]] <- current
  }
}

valid <- product_tables[vapply(product_tables, function(x) !is.null(x) && nrow(x) > 0L, logical(1))]
products <- if (length(valid) > 0L) do.call(rbind, valid) else empty_products()
if (nrow(products) > 0L) products$commercial_product_id <- paste0("PUBCHEM_VENDOR_", sprintf("%06d", seq_len(nrow(products))))

compound_summary <- data.frame(
  supplier_lookup_id = manifest$supplier_lookup_id, compound_id = manifest$compound_id,
  compound_name = manifest$compound_name, compound_origin = manifest$compound_origin,
  chembl_id = manifest$chembl_id, original_pubchem_cid = safe_numeric(manifest$pubchem_cid),
  resolved_pubchem_cid = resolved, full_inchi_key = manifest$full_inchi_key,
  parent_connectivity_key = manifest$parent_connectivity_key,
  pubchem_vendor_lookup_status = lookup_status, commercial_availability_status = NA_character_,
  public_vendor_information_present = FALSE, supplier_count = 0L,
  commercial_product_record_count = 0L, supplier_names = NA_character_, product_urls = NA_character_,
  exact_product_present = NA, identity_review_required = FALSE,
  date_checked = as.character(Sys.Date()), stringsAsFactors = FALSE
)

for (i in seq_len(nrow(compound_summary))) {
  current <- products[products$compound_id == compound_summary$compound_id[i], , drop = FALSE]
  if (nrow(current) > 0L) {
    compound_summary$public_vendor_information_present[i] <- TRUE
    compound_summary$supplier_count[i] <- length(unique(current$supplier_name))
    compound_summary$commercial_product_record_count[i] <- nrow(current)
    compound_summary$supplier_names[i] <- collapse_unique(current$supplier_name)
    compound_summary$product_urls[i] <- collapse_unique(current$product_url)
    compound_summary$identity_review_required[i] <- TRUE
  }
}
compound_summary$commercial_availability_status <- ifelse(
  compound_summary$public_vendor_information_present,
  "public_supplier_names_available",
  ifelse(is.na(compound_summary$resolved_pubchem_cid), "pubchem_identity_not_resolved", "no_public_supplier_names_found")
)

supplier_names <- sort(unique(clean_text(products$supplier_name)))
supplier_names <- supplier_names[!is.na(supplier_names)]
supplier_directory <- if (length(supplier_names) == 0L) data.frame(
  supplier_id = character(), supplier_name = character(), compound_count = integer(),
  listing_count = integer(), product_urls = character(), data_source = character(),
  last_checked = character(), stringsAsFactors = FALSE
) else do.call(rbind, lapply(seq_along(supplier_names), function(i) {
  current <- products[products$supplier_name == supplier_names[i], , drop = FALSE]
  data.frame(
    supplier_id = paste0("SUPPLIER_", sprintf("%05d", i)), supplier_name = supplier_names[i],
    compound_count = length(unique(current$compound_id)), listing_count = nrow(current),
    product_urls = collapse_unique(current$product_url), data_source = "PubChem Chemical Vendors",
    last_checked = as.character(Sys.Date()), stringsAsFactors = FALSE
  )
}))

identity_review <- data.frame(
  compound_id = compound_summary$compound_id, compound_name = compound_summary$compound_name,
  resolved_pubchem_cid = compound_summary$resolved_pubchem_cid,
  full_inchi_key = compound_summary$full_inchi_key,
  public_vendor_information_present = compound_summary$public_vendor_information_present,
  supplier_count = compound_summary$supplier_count,
  identity_review_status = ifelse(compound_summary$public_vendor_information_present, "manual_product_identity_review_required", "no_public_supplier_record_to_review"),
  review_reason = ifelse(compound_summary$public_vendor_information_present,
    "Supplier names were retrieved from PubChem references; exact product identity, form, purity, stock, price and lead time require manual verification.",
    "No public PubChem supplier reference was found."),
  manual_review_completed = FALSE, reviewer = NA_character_, review_date = NA_character_,
  review_notes = NA_character_, stringsAsFactors = FALSE
)

metrics <- data.frame(
  metric = c("compound_count", "pubchem_cid_resolved_count", "compound_count_with_supplier_names", "compound_count_without_supplier_names", "commercial_product_record_count", "unique_supplier_count", "manual_identity_review_required_count", "manual_identity_reviews_completed"),
  value = c(nrow(manifest), sum(!is.na(resolved)), sum(compound_summary$public_vendor_information_present), sum(!compound_summary$public_vendor_information_present), nrow(products), nrow(supplier_directory), sum(identity_review$identity_review_status == "manual_product_identity_review_required"), 0),
  stringsAsFactors = FALSE
)

cat("Writing A18B supplier-aware outputs...\n")
safe_write_csv(compound_summary, outputs[["summary"]])
safe_write_csv(products, outputs[["products"]])
safe_write_csv(supplier_directory, outputs[["suppliers"]])
safe_write_csv(identity_review, outputs[["review"]])
safe_write_csv(metrics, outputs[["metrics"]])

report <- paste0(
  "<!doctype html><html><head><meta charset='utf-8'><title>RESKO A18B Supplier Names</title>",
  "<style>body{font-family:Arial,sans-serif;max-width:1300px;margin:auto;padding:25px}table{border-collapse:collapse;width:100%;font-size:12px}th,td{border:1px solid #ddd;padding:6px;text-align:left}th{background:#eee}</style></head><body>",
  "<h1>RESKO A18B PubChem Supplier Names</h1><p>Supplier references require manual product verification.</p>",
  make_html_table(compound_summary[, c("compound_name", "resolved_pubchem_cid", "commercial_availability_status", "supplier_count", "supplier_names", "date_checked")]),
  "<h2>Supplier directory</h2>", make_html_table(supplier_directory), "</body></html>"
)
tmp_report <- tempfile(pattern = "A18B_report_", tmpdir = results_dir, fileext = ".html")
writeLines(report, tmp_report, useBytes = TRUE)
backup_file(outputs[["report"]])
if (!file.rename(tmp_report, outputs[["report"]])) stop("Could not install report.", call. = FALSE)

for (path in outputs) if (!file.exists(path) || file.info(path)$size <= 0L) stop("Output validation failed: ", path, call. = FALSE)
check_summary <- readr::read_csv(outputs[["summary"]], show_col_types = FALSE)
check_metrics <- readr::read_csv(outputs[["metrics"]], show_col_types = FALSE)
if (nrow(check_summary) != 16L) stop("Commercial summary must contain 16 rows.", call. = FALSE)
if (nrow(check_metrics) != 8L) stop("A18B metrics must contain 8 rows.", call. = FALSE)

cat("\nA18B supplier-name enrichment completed.\n")
cat("Compounds processed: ", nrow(manifest), "\n", sep = "")
cat("PubChem CIDs resolved: ", sum(!is.na(resolved)), "\n", sep = "")
cat("Compounds with supplier names: ", sum(compound_summary$public_vendor_information_present), "\n", sep = "")
cat("Commercial supplier records: ", nrow(products), "\n", sep = "")
cat("Unique suppliers: ", nrow(supplier_directory), "\n", sep = "")
cat("Verified outputs: ", length(outputs), "\n", sep = "")
