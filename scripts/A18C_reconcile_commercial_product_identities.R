#!/usr/bin/env Rscript

# ============================================================
# RESKO A18C: Reconcile commercial-product identities
# ============================================================
#
# Run from the RESKO project root:
#
# Rscript --vanilla scripts/A18C_reconcile_commercial_product_identities.R
#
# Inputs:
#
# results/A18A_compound_identifiers.csv
# results/A18B_compound_commercial_summary.csv
# results/A18B_commercial_products.csv
#
# Outputs:
#
# results/A18C_commercial_products_reconciled.csv
# results/A18C_commercial_identity_review.csv
# results/A18C_compound_commercial_summary_reconciled.csv
# results/A18C_summary.csv
# results/A18C_commercial_identity_report.html
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

required_packages <- c(
  "readr"
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

if (
  !dir.exists(results_dir) ||
  !dir.exists(scripts_dir)
) {
  stop(
    "Run A18C from the RESKO project root.",
    call. = FALSE
  )
}

input_identifiers <- file.path(
  results_dir,
  "A18A_compound_identifiers.csv"
)

input_commercial_summary <- file.path(
  results_dir,
  "A18B_compound_commercial_summary.csv"
)

input_products <- file.path(
  results_dir,
  "A18B_commercial_products.csv"
)

outputs <- c(
  reconciled_products = file.path(
    results_dir,
    "A18C_commercial_products_reconciled.csv"
  ),
  identity_review = file.path(
    results_dir,
    "A18C_commercial_identity_review.csv"
  ),
  reconciled_summary = file.path(
    results_dir,
    "A18C_compound_commercial_summary_reconciled.csv"
  ),
  summary = file.path(
    results_dir,
    "A18C_summary.csv"
  ),
  report = file.path(
    results_dir,
    "A18C_commercial_identity_report.html"
  )
)

for (path in c(
  input_identifiers,
  input_commercial_summary,
  input_products
)) {
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

clean_text <- function(x) {
  if (
    is.null(x) ||
    length(x) == 0L
  ) {
    return(NA_character_)
  }

  output <- as.character(x)

  output[
    is.na(output) |
      trimws(output) == ""
  ] <- NA_character_

  output
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

safe_logical <- function(x) {
  values <- tolower(
    clean_text(x)
  )

  ifelse(
    values %in% c(
      "true",
      "t",
      "1",
      "yes"
    ),
    TRUE,
    ifelse(
      values %in% c(
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

html_escape <- function(x) {
  output <- as.character(x)
  output <- gsub("&", "&amp;", output, fixed = TRUE)
  output <- gsub("<", "&lt;", output, fixed = TRUE)
  output <- gsub(">", "&gt;", output, fixed = TRUE)
  output <- gsub("\"", "&quot;", output, fixed = TRUE)
  output
}

backup_file <- function(path) {
  if (!file.exists(path)) {
    return(invisible(NULL))
  }

  stamp <- format(
    Sys.time(),
    "%Y%m%d_%H%M%S"
  )

  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)

  backup_path <- paste0(
    stem,
    "_previous_",
    stamp,
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

safe_write_csv <- function(data, path) {
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
      logical(1)
    )
  ]

  if (length(list_columns) > 0L) {
    stop(
      "Cannot write list column(s): ",
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
    file.info(temporary_path)$size <= 0L
  ) {
    stop(
      "Could not create output: ",
      path,
      call. = FALSE
    )
  }

  backup_file(path)

  moved <- file.rename(
    temporary_path,
    path
  )

  if (!moved) {
    stop(
      "Could not install output: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}

make_html_table <- function(data) {
  table_data <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  if (nrow(table_data) == 0L) {
    return(
      "<p>No records available.</p>"
    )
  }

  for (column_name in names(table_data)) {
    table_data[[column_name]] <- ifelse(
      is.na(table_data[[column_name]]),
      "",
      html_escape(
        table_data[[column_name]]
      )
    )
  }

  header <- paste0(
    "<tr>",
    paste0(
      "<th>",
      html_escape(names(table_data)),
      "</th>",
      collapse = ""
    ),
    "</tr>"
  )

  rows <- character(
    nrow(table_data)
  )

  for (index in seq_len(nrow(table_data))) {
    rows[index] <- paste0(
      "<tr>",
      paste0(
        "<td>",
        unlist(
          table_data[
            index,
            ,
            drop = FALSE
          ]
        ),
        "</td>",
        collapse = ""
      ),
      "</tr>"
    )
  }

  paste0(
    "<table><thead>",
    header,
    "</thead><tbody>",
    paste(rows, collapse = "\n"),
    "</tbody></table>"
  )
}

classify_product_text <- function(
  vendor_heading,
  information_name,
  information_text,
  section_path
) {
  combined_text <- paste(
    clean_text(vendor_heading),
    clean_text(information_name),
    clean_text(information_text),
    clean_text(section_path),
    collapse = " "
  )

  combined_text <- tolower(
    combined_text
  )

  if (
    grepl(
      paste(
        "isotope",
        "isotopic",
        "isotopically",
        "deuterated",
        "radiolabel",
        "radiolabeled",
        "radiolabelled",
        "13c",
        "14c",
        "15n",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "isotopically_labelled"
    )
  }

  if (
    grepl(
      paste(
        "analytical standard",
        "reference standard",
        "certified reference",
        "certified standard",
        "crm",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "analytical_or_reference_standard"
    )
  }

  if (
    grepl(
      paste(
        "mixture",
        "formulation",
        "solution",
        "suspension",
        "tablet",
        "capsule",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "mixture_or_formulation"
    )
  }

  if (
    grepl(
      paste(
        "hydrochloride",
        " hcl",
        "sodium",
        "potassium",
        "mesylate",
        "besylate",
        "tosylate",
        "tartrate",
        "sulfate",
        "phosphate",
        "acetate",
        "citrate",
        "fumarate",
        "maleate",
        "succinate",
        "oxalate",
        "lactate",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "salt_form_requires_review"
    )
  }

  if (
    grepl(
      paste(
        "hydrate",
        "solvate",
        "monohydrate",
        "dihydrate",
        "hemihydrate",
        "ethanolate",
        "methanolate",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "solvate_or_hydrate_requires_review"
    )
  }

  if (
    grepl(
      paste(
        "stereoisomer",
        "racemate",
        "racemic",
        "enantiomer",
        "diastereomer",
        sep = "|"
      ),
      combined_text
    )
  ) {
    return(
      "stereochemical_identity_requires_review"
    )
  }

  "public_listing_identity_requires_review"
}

# ============================================================
# Read inputs
# ============================================================

cat(
  "Reading A18C inputs...\n"
)

identifiers <- readr::read_csv(
  input_identifiers,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

commercial_summary <- readr::read_csv(
  input_commercial_summary,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

products <- readr::read_csv(
  input_products,
  col_types = readr::cols(
    .default = readr::col_character()
  ),
  progress = FALSE
)

if (nrow(identifiers) != 16L) {
  stop(
    "A18A identifier manifest must contain 16 compounds.",
    call. = FALSE
  )
}

if (nrow(commercial_summary) != 16L) {
  stop(
    "A18B commercial summary must contain 16 compounds.",
    call. = FALSE
  )
}

required_product_columns <- c(
  "commercial_product_id",
  "compound_id",
  "compound_name",
  "pubchem_cid",
  "vendor_heading",
  "information_name",
  "information_text",
  "product_url",
  "section_path"
)

missing_product_columns <- setdiff(
  required_product_columns,
  names(products)
)

if (length(missing_product_columns) > 0L) {
  stop(
    "A18B product table is missing column(s): ",
    paste(
      missing_product_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

# ============================================================
# Reconcile product records
# ============================================================

cat(
  "Classifying commercial product records...\n"
)

if (nrow(products) > 0L) {
  identity_classes <- character(
    nrow(products)
  )

  for (index in seq_len(nrow(products))) {
    identity_classes[index] <- classify_product_text(
      vendor_heading =
        products$vendor_heading[index],

      information_name =
        products$information_name[index],

      information_text =
        products$information_text[index],

      section_path =
        products$section_path[index]
    )
  }

  products$identity_match_class <-
    identity_classes

  products$exact_identity_confirmed <-
    FALSE

  products$scientific_parent_form_confirmed <-
    FALSE

  products$manual_review_required <-
    TRUE

  products$manual_verification_status <-
    "not_verified"

  products$manual_review_notes <-
    NA_character_

  identifier_match <- match(
    products$compound_id,
    identifiers$compound_id
  )

  products$scientific_full_inchi_key <-
    identifiers$full_inchi_key[
      identifier_match
    ]

  products$scientific_parent_connectivity_key <-
    identifiers$parent_connectivity_key[
      identifier_match
    ]

  products$scientific_canonical_smiles <-
    identifiers$canonical_smiles[
      identifier_match
    ]

  reconciled_products <- products
} else {
  reconciled_products <- products

  reconciled_products$exact_identity_confirmed <-
    logical()

  reconciled_products$scientific_parent_form_confirmed <-
    logical()

  reconciled_products$manual_review_required <-
    logical()

  reconciled_products$manual_review_notes <-
    character()

  reconciled_products$scientific_full_inchi_key <-
    character()

  reconciled_products$scientific_parent_connectivity_key <-
    character()

  reconciled_products$scientific_canonical_smiles <-
    character()
}

# ============================================================
# Build compound-level reconciliation
# ============================================================

cat(
  "Building compound-level identity review...\n"
)

identity_review <- data.frame(
  compound_id =
    commercial_summary$compound_id,

  compound_name =
    commercial_summary$compound_name,

  resolved_pubchem_cid =
    commercial_summary$resolved_pubchem_cid,

  full_inchi_key =
    commercial_summary$full_inchi_key,

  commercial_availability_status =
    commercial_summary$
      commercial_availability_status,

  public_vendor_information_present =
    safe_logical(
      commercial_summary$
        public_vendor_information_present
    ),

  supplier_count =
    suppressWarnings(
      as.numeric(
        commercial_summary$supplier_count
      )
    ),

  commercial_product_record_count =
    suppressWarnings(
      as.numeric(
        commercial_summary$
          commercial_product_record_count
      )
    ),

  listing_identity_classes =
    NA_character_,

  exact_identity_confirmed =
    FALSE,

  manual_product_review_required =
    FALSE,

  procurement_readiness =
    "no_verified_product",

  reviewer =
    NA_character_,

  review_date =
    NA_character_,

  review_notes =
    NA_character_,

  stringsAsFactors = FALSE
)

if (nrow(reconciled_products) > 0L) {
  for (index in seq_len(nrow(identity_review))) {
    current_id <- identity_review$compound_id[index]

    current_products <- reconciled_products[
      reconciled_products$compound_id ==
        current_id,
      ,
      drop = FALSE
    ]

    if (nrow(current_products) == 0L) {
      next
    }

    identity_review$
      listing_identity_classes[index] <-
      collapse_unique(
        current_products$
          identity_match_class
      )

    identity_review$
      exact_identity_confirmed[index] <-
      any(
        current_products$
          exact_identity_confirmed,
        na.rm = TRUE
      )

    identity_review$
      manual_product_review_required[index] <-
      TRUE

    identity_review$
      procurement_readiness[index] <-
      "manual_product_identity_review_required"
  }
}

reconciled_summary <- commercial_summary

summary_match <- match(
  reconciled_summary$compound_id,
  identity_review$compound_id
)

reconciled_summary$listing_identity_classes <-
  identity_review$listing_identity_classes[
    summary_match
  ]

reconciled_summary$exact_identity_confirmed <-
  identity_review$exact_identity_confirmed[
    summary_match
  ]

reconciled_summary$manual_product_review_required <-
  identity_review$manual_product_review_required[
    summary_match
  ]

reconciled_summary$procurement_readiness <-
  identity_review$procurement_readiness[
    summary_match
  ]

# ============================================================
# Create summary metrics
# ============================================================

cat(
  "Creating A18C summary metrics...\n"
)

count_class <- function(class_name) {
  if (nrow(reconciled_products) == 0L) {
    return(0L)
  }

  sum(
    reconciled_products$
      identity_match_class ==
      class_name,
    na.rm = TRUE
  )
}

summary_table <- data.frame(
  metric = c(
    "compound_count",
    "commercial_product_record_count",
    "compound_count_with_public_vendor_information",
    "compound_count_requiring_manual_product_review",
    "exact_identity_confirmed_count",
    "public_listing_identity_review_count",
    "salt_form_review_count",
    "solvate_or_hydrate_review_count",
    "stereochemical_review_count",
    "mixture_or_formulation_count",
    "analytical_or_reference_standard_count",
    "isotopically_labelled_count"
  ),

  value = c(
    as.character(
      nrow(commercial_summary)
    ),

    as.character(
      nrow(reconciled_products)
    ),

    as.character(
      sum(
        identity_review$
          public_vendor_information_present,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        identity_review$
          manual_product_review_required,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        identity_review$
          exact_identity_confirmed,
        na.rm = TRUE
      )
    ),

    as.character(
      count_class(
        "public_listing_identity_requires_review"
      )
    ),

    as.character(
      count_class(
        "salt_form_requires_review"
      )
    ),

    as.character(
      count_class(
        "solvate_or_hydrate_requires_review"
      )
    ),

    as.character(
      count_class(
        "stereochemical_identity_requires_review"
      )
    ),

    as.character(
      count_class(
        "mixture_or_formulation"
      )
    ),

    as.character(
      count_class(
        "analytical_or_reference_standard"
      )
    ),

    as.character(
      count_class(
        "isotopically_labelled"
      )
    )
  ),

  stringsAsFactors = FALSE
)

# ============================================================
# Write outputs
# ============================================================

cat(
  "Writing A18C outputs...\n"
)

safe_write_csv(
  reconciled_products,
  outputs[["reconciled_products"]]
)

safe_write_csv(
  identity_review,
  outputs[["identity_review"]]
)

safe_write_csv(
  reconciled_summary,
  outputs[["reconciled_summary"]]
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

# ============================================================
# Create report
# ============================================================

cat(
  "Creating A18C report...\n"
)

report_columns <- c(
  "compound_name",
  "commercial_availability_status",
  "supplier_count",
  "commercial_product_record_count",
  "listing_identity_classes",
  "exact_identity_confirmed",
  "manual_product_review_required",
  "procurement_readiness"
)

report_data <- reconciled_summary[
  ,
  report_columns,
  drop = FALSE
]

report_html <- paste0(
  "<!doctype html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<title>RESKO A18C Commercial Identity Reconciliation</title>",
  "<style>",
  "body{font-family:Arial,sans-serif;max-width:1250px;",
  "margin:0 auto;padding:28px;line-height:1.5;color:#222}",
  "table{border-collapse:collapse;width:100%;font-size:13px}",
  "th,td{border:1px solid #ddd;padding:7px;text-align:left}",
  "th{background:#f3f4f6}",
  ".note{background:#fff4e5;border-left:5px solid #e67e22;",
  "padding:14px;margin:16px 0}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>RESKO A18C Commercial Identity Reconciliation</h1>",
  "<div class='note'>",
  "No commercial product was treated as an exact scientific-compound ",
  "match without manual structure and product-form verification.",
  "</div>",
  make_html_table(report_data),
  "<h2>Summary</h2>",
  make_html_table(summary_table),
  "</body>",
  "</html>"
)

temporary_report <- tempfile(
  pattern = "A18C_report_",
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
  file.info(temporary_report)$size <= 0L
) {
  stop(
    "The A18C report was not created.",
    call. = FALSE
  )
}

backup_file(
  outputs[["report"]]
)

if (
  !file.rename(
    temporary_report,
    outputs[["report"]]
  )
) {
  stop(
    "Could not install the A18C report.",
    call. = FALSE
  )
}

# ============================================================
# Validate outputs
# ============================================================

cat(
  "Validating A18C outputs...\n"
)

for (path in outputs) {
  if (
    !file.exists(path) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Output validation failed: ",
      path,
      call. = FALSE
    )
  }
}

identity_check <- readr::read_csv(
  outputs[["identity_review"]],
  show_col_types = FALSE
)

summary_check <- readr::read_csv(
  outputs[["summary"]],
  show_col_types = FALSE
)

reconciled_summary_check <- readr::read_csv(
  outputs[["reconciled_summary"]],
  show_col_types = FALSE
)

if (nrow(identity_check) != 16L) {
  stop(
    "A18C identity review must contain 16 compounds.",
    call. = FALSE
  )
}

if (nrow(reconciled_summary_check) != 16L) {
  stop(
    "A18C reconciled summary must contain 16 compounds.",
    call. = FALSE
  )
}

if (nrow(summary_check) != 12L) {
  stop(
    "A18C summary must contain 12 rows.",
    call. = FALSE
  )
}

cat(
  "All A18C output validations passed.\n"
)

# ============================================================
# Completion
# ============================================================

cat(
  "\nA18C commercial identity reconciliation completed.\n"
)

cat(
  "Compounds reviewed: ",
  nrow(identity_review),
  "\n",
  sep = ""
)

cat(
  "Commercial product records classified: ",
  nrow(reconciled_products),
  "\n",
  sep = ""
)

cat(
  "Compounds requiring manual product review: ",
  sum(
    identity_review$
      manual_product_review_required,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Exact product identities confirmed automatically: ",
  sum(
    identity_review$
      exact_identity_confirmed,
    na.rm = TRUE
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