#!/usr/bin/env Rscript

# ============================================================
# RESKO A18D: Prepare live-network compound detail data
# ============================================================
#
# Run:
#
# Rscript --vanilla scripts/A18D_prepare_live_network_compound_data.R
#
# Required inputs:
#
# results/A18A_compound_detail_manifest.csv
# results/A18A_compound_evidence_summary.csv
# results/A18C_compound_commercial_summary_reconciled.csv
# results/A18C_commercial_products_reconciled.csv
#
# Optional inputs:
#
# results/A17D_candidate_side_effects.csv
# results/A17D_candidate_indications.csv
# results/A17D_side_effect_similarity.csv
#
# Outputs:
#
# results/A18D_live_compound_details.csv
# results/A18D_live_compound_details.json
# results/A18D_live_commercial_products.json
# results/A18D_live_compound_aliases.csv
# results/A18D_live_compound_aliases.json
# results/A18D_summary.csv
# results/A18D_live_network_data_report.html
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

required_packages <- c(
  "readr",
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
    "Run A18D from the RESKO project root.",
    call. = FALSE
  )
}

input_details <- file.path(
  results_dir,
  "A18A_compound_detail_manifest.csv"
)

input_evidence <- file.path(
  results_dir,
  "A18A_compound_evidence_summary.csv"
)

input_commercial <- file.path(
  results_dir,
  "A18C_compound_commercial_summary_reconciled.csv"
)

input_products <- file.path(
  results_dir,
  "A18C_commercial_products_reconciled.csv"
)

input_side_effects <- file.path(
  results_dir,
  "A17D_candidate_side_effects.csv"
)

input_indications <- file.path(
  results_dir,
  "A17D_candidate_indications.csv"
)

input_similarity <- file.path(
  results_dir,
  "A17D_side_effect_similarity.csv"
)

outputs <- c(
  compound_csv = file.path(
    results_dir,
    "A18D_live_compound_details.csv"
  ),
  compound_json = file.path(
    results_dir,
    "A18D_live_compound_details.json"
  ),
  products_json = file.path(
    results_dir,
    "A18D_live_commercial_products.json"
  ),
  aliases_csv = file.path(
    results_dir,
    "A18D_live_compound_aliases.csv"
  ),
  aliases_json = file.path(
    results_dir,
    "A18D_live_compound_aliases.json"
  ),
  summary = file.path(
    results_dir,
    "A18D_summary.csv"
  ),
  report = file.path(
    results_dir,
    "A18D_live_network_data_report.html"
  )
)

required_inputs <- c(
  input_details,
  input_evidence,
  input_commercial,
  input_products
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
  candidates,
  default = NA_character_
) {
  available <- candidates[
    candidates %in% names(data)
  ]

  if (length(available) == 0L) {
    return(
      rep(
        default,
        nrow(data)
      )
    )
  }

  data[[available[1]]]
}

read_character_csv <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(
      .default = readr::col_character()
    ),
    progress = FALSE
  )
}

read_optional_csv <- function(path) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    return(
      data.frame(
        stringsAsFactors = FALSE
      )
    )
  }

  read_character_csv(path)
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

  readr::write_csv(
    as.data.frame(
      data,
      stringsAsFactors = FALSE
    ),
    temporary_path,
    na = ""
  )

  if (
    !file.exists(temporary_path) ||
    file.info(temporary_path)$size <= 0L
  ) {
    stop(
      "Failed to create output: ",
      path,
      call. = FALSE
    )
  }

  backup_file(path)

  if (!file.rename(temporary_path, path)) {
    stop(
      "Could not install output: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}

safe_write_json <- function(
  object,
  path,
  dataframe = "rows"
) {
  temporary_path <- paste0(
    path,
    ".tmp"
  )

  if (file.exists(temporary_path)) {
    unlink(temporary_path)
  }

  jsonlite::write_json(
    object,
    path = temporary_path,
    pretty = TRUE,
    auto_unbox = TRUE,
    dataframe = dataframe,
    na = "null",
    null = "null"
  )

  if (
    !file.exists(temporary_path) ||
    file.info(temporary_path)$size <= 0L
  ) {
    stop(
      "Failed to create JSON output: ",
      path,
      call. = FALSE
    )
  }

  backup_file(path)

  if (!file.rename(temporary_path, path)) {
    stop(
      "Could not install JSON output: ",
      path,
      call. = FALSE
    )
  }

  invisible(path)
}

html_escape <- function(x) {
  output <- as.character(x)
  output <- gsub("&", "&amp;", output, fixed = TRUE)
  output <- gsub("<", "&lt;", output, fixed = TRUE)
  output <- gsub(">", "&gt;", output, fixed = TRUE)
  output <- gsub("\"", "&quot;", output, fixed = TRUE)
  output
}

make_html_table <- function(data) {
  data <- as.data.frame(
    data,
    stringsAsFactors = FALSE
  )

  if (nrow(data) == 0L) {
    return(
      "<p>No records available.</p>"
    )
  }

  for (column_name in names(data)) {
    data[[column_name]] <- ifelse(
      is.na(data[[column_name]]),
      "",
      html_escape(
        data[[column_name]]
      )
    )
  }

  header <- paste0(
    "<tr>",
    paste0(
      "<th>",
      html_escape(names(data)),
      "</th>",
      collapse = ""
    ),
    "</tr>"
  )

  rows <- character(
    nrow(data)
  )

  for (index in seq_len(nrow(data))) {
    rows[index] <- paste0(
      "<tr>",
      paste0(
        "<td>",
        unlist(
          data[
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

# ============================================================
# Read inputs
# ============================================================

cat(
  "Reading A18D inputs...\n"
)

details <- read_character_csv(
  input_details
)

evidence <- read_character_csv(
  input_evidence
)

commercial <- read_character_csv(
  input_commercial
)

products <- read_character_csv(
  input_products
)

side_effects <- read_optional_csv(
  input_side_effects
)

indications <- read_optional_csv(
  input_indications
)

similarity <- read_optional_csv(
  input_similarity
)

if (nrow(details) != 16L) {
  stop(
    "A18A detail manifest must contain 16 compounds.",
    call. = FALSE
  )
}

if (nrow(evidence) != 16L) {
  stop(
    "A18A evidence summary must contain 16 compounds.",
    call. = FALSE
  )
}

if (nrow(commercial) != 16L) {
  stop(
    "A18C commercial summary must contain 16 compounds.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    details$compound_id
  ) > 0L
) {
  stop(
    "Duplicate compound IDs were detected.",
    call. = FALSE
  )
}

# ============================================================
# Build side-effect summaries
# ============================================================

cat(
  "Summarising side-effect information...\n"
)

side_effect_summary <- data.frame(
  normalised_compound_name = character(),
  side_effect_record_count = integer(),
  side_effect_names = character(),
  stringsAsFactors = FALSE
)

if (nrow(side_effects) > 0L) {
  side_effect_candidate <- clean_text(
    get_column(
      side_effects,
      c(
        "candidate_name",
        "sider_drug_name"
      )
    )
  )

  side_effect_name <- clean_text(
    get_column(
      side_effects,
      c(
        "side_effect_name",
        "meddra_concept_name"
      )
    )
  )

  side_effect_work <- data.frame(
    normalised_compound_name =
      normalise_name(
        side_effect_candidate
      ),

    side_effect_name =
      side_effect_name,

    stringsAsFactors = FALSE
  )

  valid_names <- sort(
    unique(
      side_effect_work$
        normalised_compound_name[
          !is.na(
            side_effect_work$
              normalised_compound_name
          )
        ]
    )
  )

  side_effect_rows <- vector(
    "list",
    length(valid_names)
  )

  for (index in seq_along(valid_names)) {
    current_name <- valid_names[index]

    current_rows <- side_effect_work[
      side_effect_work$
        normalised_compound_name ==
        current_name,
      ,
      drop = FALSE
    ]

    side_effect_rows[[index]] <- data.frame(
      normalised_compound_name =
        current_name,

      side_effect_record_count =
        nrow(current_rows),

      side_effect_names =
        collapse_unique(
          current_rows$side_effect_name
        ),

      stringsAsFactors = FALSE
    )
  }

  side_effect_summary <- do.call(
    rbind,
    side_effect_rows
  )
}

# ============================================================
# Build indication summaries
# ============================================================

cat(
  "Summarising indication information...\n"
)

indication_summary <- data.frame(
  normalised_compound_name = character(),
  indication_record_count = integer(),
  indication_names = character(),
  stringsAsFactors = FALSE
)

if (nrow(indications) > 0L) {
  indication_candidate <- clean_text(
    get_column(
      indications,
      c(
        "candidate_name",
        "sider_drug_name"
      )
    )
  )

  indication_name <- clean_text(
    get_column(
      indications,
      c(
        "indication_name",
        "meddra_concept_name"
      )
    )
  )

  indication_work <- data.frame(
    normalised_compound_name =
      normalise_name(
        indication_candidate
      ),

    indication_name =
      indication_name,

    stringsAsFactors = FALSE
  )

  valid_names <- sort(
    unique(
      indication_work$
        normalised_compound_name[
          !is.na(
            indication_work$
              normalised_compound_name
          )
        ]
    )
  )

  indication_rows <- vector(
    "list",
    length(valid_names)
  )

  for (index in seq_along(valid_names)) {
    current_name <- valid_names[index]

    current_rows <- indication_work[
      indication_work$
        normalised_compound_name ==
        current_name,
      ,
      drop = FALSE
    ]

    indication_rows[[index]] <- data.frame(
      normalised_compound_name =
        current_name,

      indication_record_count =
        nrow(current_rows),

      indication_names =
        collapse_unique(
          current_rows$indication_name
        ),

      stringsAsFactors = FALSE
    )
  }

  indication_summary <- do.call(
    rbind,
    indication_rows
  )
}

# ============================================================
# Build similarity summaries
# ============================================================

cat(
  "Summarising side-effect-profile similarities...\n"
)

similarity_summary <- data.frame(
  normalised_compound_name = character(),
  closest_side_effect_neighbour = character(),
  highest_side_effect_jaccard = numeric(),
  highest_side_effect_cosine = numeric(),
  stringsAsFactors = FALSE
)

if (nrow(similarity) > 0L) {
  drug_1 <- clean_text(
    get_column(
      similarity,
      c("drug_1")
    )
  )

  drug_2 <- clean_text(
    get_column(
      similarity,
      c("drug_2")
    )
  )

  jaccard <- safe_numeric(
    get_column(
      similarity,
      c("jaccard_similarity")
    )
  )

  cosine <- safe_numeric(
    get_column(
      similarity,
      c("tfidf_cosine_similarity")
    )
  )

  long_similarity <- rbind(
    data.frame(
      compound_name = drug_1,
      neighbour_name = drug_2,
      jaccard = jaccard,
      cosine = cosine,
      stringsAsFactors = FALSE
    ),
    data.frame(
      compound_name = drug_2,
      neighbour_name = drug_1,
      jaccard = jaccard,
      cosine = cosine,
      stringsAsFactors = FALSE
    )
  )

  long_similarity$normalised_compound_name <-
    normalise_name(
      long_similarity$compound_name
    )

  valid_names <- sort(
    unique(
      long_similarity$
        normalised_compound_name[
          !is.na(
            long_similarity$
              normalised_compound_name
          )
        ]
    )
  )

  similarity_rows <- vector(
    "list",
    length(valid_names)
  )

  for (index in seq_along(valid_names)) {
    current_name <- valid_names[index]

    current_rows <- long_similarity[
      long_similarity$
        normalised_compound_name ==
        current_name,
      ,
      drop = FALSE
    ]

    current_rows <- current_rows[
      order(
        -current_rows$cosine,
        -current_rows$jaccard,
        na.last = TRUE
      ),
      ,
      drop = FALSE
    ]

    best_row <- current_rows[
      1,
      ,
      drop = FALSE
    ]

    similarity_rows[[index]] <- data.frame(
      normalised_compound_name =
        current_name,

      closest_side_effect_neighbour =
        best_row$neighbour_name,

      highest_side_effect_jaccard =
        best_row$jaccard,

      highest_side_effect_cosine =
        best_row$cosine,

      stringsAsFactors = FALSE
    )
  }

  similarity_summary <- do.call(
    rbind,
    similarity_rows
  )
}

# ============================================================
# Combine compound details
# ============================================================

cat(
  "Combining compound, evidence, and commercial fields...\n"
)

detail_match <- match(
  details$compound_id,
  evidence$compound_id
)

commercial_match <- match(
  details$compound_id,
  commercial$compound_id
)

side_effect_match <- match(
  details$normalised_compound_name,
  side_effect_summary$
    normalised_compound_name
)

indication_match <- match(
  details$normalised_compound_name,
  indication_summary$
    normalised_compound_name
)

similarity_match <- match(
  details$normalised_compound_name,
  similarity_summary$
    normalised_compound_name
)

compound_details <- data.frame(
  compound_id =
    details$compound_id,

  node_lookup_id =
    details$compound_id,

  compound_name =
    details$compound_name,

  node_type =
    "compound",

  compound_class =
    details$compound_class,

  compound_origin =
    details$compound_origin,

  chembl_id =
    details$chembl_id,

  pubchem_cid =
    commercial$resolved_pubchem_cid[
      commercial_match
    ],

  full_inchi_key =
    details$full_inchi_key,

  parent_connectivity_key =
    details$parent_connectivity_key,

  canonical_smiles =
    details$canonical_smiles,

  analysis_smiles =
    details$analysis_smiles,

  nearest_query_name =
    details$nearest_query_name,

  nearest_query_tanimoto =
    safe_numeric(
      details$nearest_query_tanimoto
    ),

  similarity_band =
    details$similarity_band,

  biological_classification =
    evidence$biological_classification[
      detail_match
    ],

  progression_status =
    evidence$progression_status[
      detail_match
    ],

  scientific_evidence_status =
    evidence$scientific_evidence_status[
      detail_match
    ],

  matched_network_proteins =
    evidence$matched_network_proteins[
      detail_match
    ],

  independent_assay_count =
    safe_numeric(
      evidence$independent_assay_count[
        detail_match
      ]
    ),

  independent_document_count =
    safe_numeric(
      evidence$independent_document_count[
        detail_match
      ]
    ),

  maximum_clinical_phase =
    safe_numeric(
      details$max_phase
    ),

  first_approval_year =
    safe_numeric(
      details$first_approval
    ),

  sider_representation =
    details$sider_representation,

  side_effect_record_count =
    side_effect_summary$
      side_effect_record_count[
        side_effect_match
      ],

  side_effect_names =
    side_effect_summary$
      side_effect_names[
        side_effect_match
      ],

  indication_record_count =
    indication_summary$
      indication_record_count[
        indication_match
      ],

  indication_names =
    indication_summary$
      indication_names[
        indication_match
      ],

  closest_side_effect_neighbour =
    similarity_summary$
      closest_side_effect_neighbour[
        similarity_match
      ],

  highest_side_effect_jaccard =
    similarity_summary$
      highest_side_effect_jaccard[
        similarity_match
      ],

  highest_side_effect_cosine =
    similarity_summary$
      highest_side_effect_cosine[
        similarity_match
      ],

  commercial_availability_status =
    commercial$
      commercial_availability_status[
        commercial_match
      ],

  public_vendor_information_present =
    safe_logical(
      commercial$
        public_vendor_information_present[
          commercial_match
        ]
    ),

  supplier_count =
    safe_numeric(
      commercial$supplier_count[
        commercial_match
      ]
    ),

  commercial_product_record_count =
    safe_numeric(
      commercial$
        commercial_product_record_count[
          commercial_match
        ]
    ),

  supplier_names =
    commercial$supplier_names[
      commercial_match
    ],

  public_product_urls =
    commercial$product_urls[
      commercial_match
    ],

  listing_identity_classes =
    commercial$listing_identity_classes[
      commercial_match
    ],

  exact_product_identity_confirmed =
    safe_logical(
      commercial$
        exact_identity_confirmed[
          commercial_match
        ]
    ),

  manual_product_review_required =
    safe_logical(
      commercial$
        manual_product_review_required[
          commercial_match
        ]
    ),

  procurement_readiness =
    commercial$procurement_readiness[
      commercial_match
    ],

  commercial_data_checked =
    commercial$date_checked[
      commercial_match
    ],

  detail_panel_status =
    "ready",

  stringsAsFactors = FALSE
)

compound_details$side_effect_record_count[
  is.na(compound_details$side_effect_record_count)
] <- 0

compound_details$indication_record_count[
  is.na(compound_details$indication_record_count)
] <- 0

compound_details$supplier_count[
  is.na(compound_details$supplier_count)
] <- 0

compound_details$commercial_product_record_count[
  is.na(
    compound_details$
      commercial_product_record_count
  )
] <- 0

if (nrow(compound_details) != 16L) {
  stop(
    "Live compound-detail table must contain 16 rows.",
    call. = FALSE
  )
}

# ============================================================
# Build alias table
# ============================================================

cat(
  "Building compound-node alias records...\n"
)

alias_rows <- list()
alias_index <- 1L

add_alias <- function(
  compound_id,
  alias,
  alias_type
) {
  alias <- clean_text(alias)

  if (
    length(alias) == 0L ||
    is.na(alias)
  ) {
    return(invisible(NULL))
  }

  alias_rows[[alias_index]] <<- data.frame(
    compound_id = compound_id,
    alias = alias,
    normalised_alias = normalise_name(alias),
    alias_type = alias_type,
    stringsAsFactors = FALSE
  )

  alias_index <<- alias_index + 1L

  invisible(NULL)
}

for (index in seq_len(nrow(compound_details))) {
  add_alias(
    compound_details$compound_id[index],
    compound_details$compound_id[index],
    "compound_id"
  )

  add_alias(
    compound_details$compound_id[index],
    compound_details$compound_name[index],
    "compound_name"
  )

  add_alias(
    compound_details$compound_id[index],
    compound_details$chembl_id[index],
    "chembl_id"
  )

  add_alias(
    compound_details$compound_id[index],
    as.character(
      compound_details$pubchem_cid[index]
    ),
    "pubchem_cid"
  )

  add_alias(
    compound_details$compound_id[index],
    compound_details$full_inchi_key[index],
    "full_inchi_key"
  )
}

compound_aliases <- do.call(
  rbind,
  alias_rows
)

compound_aliases <- unique(
  compound_aliases
)

# ============================================================
# Build product JSON structure
# ============================================================

cat(
  "Preparing commercial-product JSON records...\n"
)

if (nrow(products) > 0L) {
  products$manual_review_required <-
    safe_logical(
      products$manual_review_required
    )

  products$exact_identity_confirmed <-
    safe_logical(
      products$exact_identity_confirmed
    )
}

# ============================================================
# Build summary
# ============================================================

summary_table <- data.frame(
  metric = c(
    "compound_detail_count",
    "alias_record_count",
    "commercial_product_record_count",
    "compounds_with_public_vendor_information",
    "compounds_requiring_manual_product_review",
    "compounds_with_exact_product_identity_confirmed",
    "compounds_with_sider_information",
    "compounds_with_side_effect_records",
    "compounds_with_indication_records",
    "compound_json_record_count",
    "product_json_record_count",
    "detail_panel_ready_count"
  ),

  value = c(
    as.character(
      nrow(compound_details)
    ),

    as.character(
      nrow(compound_aliases)
    ),

    as.character(
      nrow(products)
    ),

    as.character(
      sum(
        compound_details$
          public_vendor_information_present,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        compound_details$
          manual_product_review_required,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        compound_details$
          exact_product_identity_confirmed,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        compound_details$
          sider_representation ==
          "represented_in_sider_4.1",
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        compound_details$
          side_effect_record_count > 0,
        na.rm = TRUE
      )
    ),

    as.character(
      sum(
        compound_details$
          indication_record_count > 0,
        na.rm = TRUE
      )
    ),

    as.character(
      nrow(compound_details)
    ),

    as.character(
      nrow(products)
    ),

    as.character(
      sum(
        compound_details$
          detail_panel_status ==
          "ready"
      )
    )
  ),

  stringsAsFactors = FALSE
)

# ============================================================
# Write live-network files
# ============================================================

cat(
  "Writing A18D live-network files...\n"
)

safe_write_csv(
  compound_details,
  outputs[["compound_csv"]]
)

safe_write_json(
  compound_details,
  outputs[["compound_json"]]
)

safe_write_json(
  products,
  outputs[["products_json"]]
)

safe_write_csv(
  compound_aliases,
  outputs[["aliases_csv"]]
)

safe_write_json(
  compound_aliases,
  outputs[["aliases_json"]]
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

# ============================================================
# Create report
# ============================================================

report_columns <- c(
  "compound_name",
  "compound_origin",
  "biological_classification",
  "maximum_clinical_phase",
  "sider_representation",
  "side_effect_record_count",
  "commercial_availability_status",
  "supplier_count",
  "commercial_product_record_count",
  "procurement_readiness"
)

report_data <- compound_details[
  ,
  report_columns,
  drop = FALSE
]

report_html <- paste0(
  "<!doctype html>",
  "<html lang='en'>",
  "<head>",
  "<meta charset='utf-8'>",
  "<title>RESKO A18D Live-Network Compound Data</title>",
  "<style>",
  "body{font-family:Arial,sans-serif;max-width:1250px;",
  "margin:0 auto;padding:28px;line-height:1.5;color:#222}",
  "table{border-collapse:collapse;width:100%;font-size:13px}",
  "th,td{border:1px solid #ddd;padding:7px;text-align:left}",
  "th{background:#f3f4f6}",
  ".note{background:#eaf4ff;border-left:5px solid #2176ae;",
  "padding:14px;margin:16px 0}",
  "</style>",
  "</head>",
  "<body>",
  "<h1>RESKO A18D Live-Network Compound Detail Data</h1>",
  "<div class='note'>",
  "The generated JSON records are intended to populate a compound-detail ",
  "panel after a compound node is selected in the live RESKO network.",
  "</div>",
  make_html_table(report_data),
  "<h2>Summary</h2>",
  make_html_table(summary_table),
  "</body>",
  "</html>"
)

temporary_report <- tempfile(
  pattern = "A18D_report_",
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
    "The A18D report was not created.",
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
    "Could not install the A18D report.",
    call. = FALSE
  )
}

# ============================================================
# Validate outputs
# ============================================================

cat(
  "Validating A18D outputs...\n"
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

compound_csv_check <- readr::read_csv(
  outputs[["compound_csv"]],
  show_col_types = FALSE
)

alias_check <- readr::read_csv(
  outputs[["aliases_csv"]],
  show_col_types = FALSE
)

summary_check <- readr::read_csv(
  outputs[["summary"]],
  show_col_types = FALSE
)

compound_json_check <- jsonlite::fromJSON(
  outputs[["compound_json"]]
)

product_json_check <- jsonlite::fromJSON(
  outputs[["products_json"]]
)

if (nrow(compound_csv_check) != 16L) {
  stop(
    "A18D compound details must contain 16 rows.",
    call. = FALSE
  )
}

if (nrow(compound_json_check) != 16L) {
  stop(
    "A18D compound JSON must contain 16 records.",
    call. = FALSE
  )
}

if (nrow(alias_check) == 0L) {
  stop(
    "A18D alias table is empty.",
    call. = FALSE
  )
}

if (nrow(summary_check) != 12L) {
  stop(
    "A18D summary must contain 12 rows.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    compound_csv_check$compound_id
  ) > 0L
) {
  stop(
    "A18D compound details contain duplicate IDs.",
    call. = FALSE
  )
}

cat(
  "All A18D output validations passed.\n"
)

# ============================================================
# Completion
# ============================================================

cat(
  "\nA18D live-network compound data preparation completed.\n"
)

cat(
  "Compound-detail records: ",
  nrow(compound_details),
  "\n",
  sep = ""
)

cat(
  "Compound alias records: ",
  nrow(compound_aliases),
  "\n",
  sep = ""
)

cat(
  "Commercial product records: ",
  nrow(products),
  "\n",
  sep = ""
)

cat(
  "Compounds with public vendor information: ",
  sum(
    compound_details$
      public_vendor_information_present,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Compounds requiring manual product review: ",
  sum(
    compound_details$
      manual_product_review_required,
    na.rm = TRUE
  ),
  "\n",
  sep = ""
)

cat(
  "Detail-panel-ready compounds: ",
  sum(
    compound_details$
      detail_panel_status ==
      "ready"
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