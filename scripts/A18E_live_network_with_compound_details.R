#!/usr/bin/env Rscript

# ============================================================
# RESKO A18E: Live ranked network with compound details
# ============================================================
#
# This version embeds the completed detail HTML into each
# compound node and uses visEvents() to display the content.
#
# Run:
#
# Rscript --vanilla scripts/A18E_live_network_with_compound_details.R
#
# Outputs:
#
# results/eef1a_network_A18E_compound_details.html
# results/A18E_network_compound_mapping.csv
# results/A18E_interface_validation.csv
# results/A18E_summary.csv
#
# ============================================================

options(
  stringsAsFactors = FALSE,
  warn = 1
)

required_packages <- c(
  "readr",
  "dplyr",
  "visNetwork",
  "htmlwidgets",
  "htmltools"
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
# 1. Paths
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

if (
  !dir.exists(results_dir) ||
  !dir.exists(scripts_dir)
) {
  stop(
    "Run A18E from the RESKO project root.",
    call. = FALSE
  )
}

files <- list(
  proteins = file.path(
    results_dir,
    "nodes_proteins.csv"
  ),

  pathways = file.path(
    results_dir,
    "nodes_pathways.csv"
  ),

  go = file.path(
    results_dir,
    "nodes_biological_process.csv"
  ),

  compounds = file.path(
    results_dir,
    "nodes_drugs.csv"
  ),

  ppi = file.path(
    results_dir,
    "edges_interacts_with.csv"
  ),

  protein_pathway = file.path(
    results_dir,
    "edges_protein_pathway.csv"
  ),

  protein_go = file.path(
    results_dir,
    "edges_protein_go.csv"
  ),

  compound_edges = file.path(
    results_dir,
    "edges_compound_protein_corrected.csv"
  ),

  scores = file.path(
    results_dir,
    "A15_candidate_score_components.csv"
  ),

  direct = file.path(
    results_dir,
    "A15_direct_eef1a1_ranking.csv"
  ),

  network = file.path(
    results_dir,
    "A15_translation_network_ranking.csv"
  ),

  compound_details = file.path(
    results_dir,
    "A18D_live_compound_details.csv"
  ),

  aliases = file.path(
    results_dir,
    "A18D_live_compound_aliases.csv"
  ),

  products = file.path(
    results_dir,
    "A18C_commercial_products_reconciled.csv"
  )
)

outputs <- c(
  network_html = file.path(
    results_dir,
    "eef1a_network_A18E_compound_details.html"
  ),

  mapping = file.path(
    results_dir,
    "A18E_network_compound_mapping.csv"
  ),

  validation = file.path(
    results_dir,
    "A18E_interface_validation.csv"
  ),

  summary = file.path(
    results_dir,
    "A18E_summary.csv"
  )
)

missing_files <- unlist(files)[
  !file.exists(
    unlist(files)
  )
]

if (length(missing_files) > 0L) {
  stop(
    "Missing input(s):\n",
    paste(
      missing_files,
      collapse = "\n"
    ),
    call. = FALSE
  )
}

# ============================================================
# 2. Helpers
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

normalise_alias <- function(x) {
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

escape_html <- function(x) {
  output <- as.character(x)

  output <- gsub(
    "&",
    "&amp;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    "<",
    "&lt;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    ">",
    "&gt;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    "\"",
    "&quot;",
    output,
    fixed = TRUE
  )

  output <- gsub(
    "'",
    "&#039;",
    output,
    fixed = TRUE
  )

  output
}

display_text <- function(
  value,
  fallback = "Not available"
) {
  value <- clean_text(value)

  if (
    length(value) == 0L ||
    all(is.na(value))
  ) {
    return(fallback)
  }

  escape_html(
    value[!is.na(value)][1]
  )
}

display_number <- function(
  value,
  digits = 3L,
  fallback = "Not available"
) {
  numeric_value <- safe_numeric(value)

  if (
    length(numeric_value) == 0L ||
    is.na(numeric_value[1])
  ) {
    return(fallback)
  }

  format(
    round(
      numeric_value[1],
      digits
    ),
    trim = TRUE,
    scientific = FALSE
  )
}

html_row <- function(
  label,
  value
) {
  paste0(
    "<div class='detail-row'>",
    "<div class='detail-label'>",
    escape_html(label),
    "</div>",
    "<div class='detail-value'>",
    display_text(value),
    "</div>",
    "</div>"
  )
}

html_card <- function(
  title,
  content
) {
  paste0(
    "<div class='detail-card'>",
    "<h3>",
    escape_html(title),
    "</h3>",
    content,
    "</div>"
  )
}

read_data <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    progress = FALSE
  )
}

read_character_data <- function(path) {
  readr::read_csv(
    path,
    col_types = readr::cols(
      .default = readr::col_character()
    ),
    progress = FALSE
  )
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
      "Could not preserve existing output: ",
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

# ============================================================
# 3. Read network and A18 data
# ============================================================

cat(
  "Reading A18E inputs...\n"
)

proteins <- read_data(
  files$proteins
)

pathways <- read_data(
  files$pathways
)

go_processes <- read_data(
  files$go
)

compounds <- read_data(
  files$compounds
)

ppi <- read_data(
  files$ppi
)

protein_pathway <- read_data(
  files$protein_pathway
)

protein_go <- read_data(
  files$protein_go
)

compound_edges_data <- read_data(
  files$compound_edges
)

scores <- read_data(
  files$scores
)

direct_ranking <- read_data(
  files$direct
)

network_ranking <- read_data(
  files$network
)

compound_details <- read_character_data(
  files$compound_details
)

aliases <- read_character_data(
  files$aliases
)

products <- read_character_data(
  files$products
)

if (nrow(compound_details) != 16L) {
  stop(
    "A18D compound details must contain 16 records.",
    call. = FALSE
  )
}

# ============================================================
# 4. Prepare ranked compounds
# ============================================================

cat(
  "Preparing ranked compound data...\n"
)

compound_rank_data <- compounds |>
  dplyr::left_join(
    scores,
    by = "molecule_chembl_id"
  ) |>
  dplyr::left_join(
    direct_ranking |>
      dplyr::select(
        "molecule_chembl_id",
        "direct_rank"
      ),
    by = "molecule_chembl_id"
  ) |>
  dplyr::left_join(
    network_ranking |>
      dplyr::select(
        "molecule_chembl_id",
        "network_rank"
      ),
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    name =
      dplyr::coalesce(
        .data$name,
        .data$molecule_chembl_id
      ),

    overall_priority_score =
      dplyr::coalesce(
        safe_numeric(
          .data$overall_priority_score
        ),
        0
      ),

    compound_size =
      18 +
      1.6 *
      .data$overall_priority_score,

    rank_label =
      paste0(
        "D",
        dplyr::coalesce(
          as.character(.data$direct_rank),
          "-"
        ),
        " / N",
        dplyr::coalesce(
          as.character(.data$network_rank),
          "-"
        ),
        "  ",
        .data$name
      )
  )

# ============================================================
# 5. Map graph compound nodes to A18D records
# ============================================================

cat(
  "Mapping network compounds to A18D details...\n"
)

alias_keys <- normalise_alias(
  aliases$alias
)

mapped_compound_ids <- rep(
  NA_character_,
  nrow(compound_rank_data)
)

for (index in seq_len(nrow(compound_rank_data))) {
  candidate_keys <- unique(
    c(
      normalise_alias(
        compound_rank_data$
          molecule_chembl_id[index]
      ),

      normalise_alias(
        compound_rank_data$name[index]
      )
    )
  )

  candidate_keys <- candidate_keys[
    !is.na(candidate_keys)
  ]

  matches <- aliases$compound_id[
    alias_keys %in%
      candidate_keys
  ]

  matches <- unique(
    clean_text(matches)
  )

  matches <- matches[
    !is.na(matches)
  ]

  if (length(matches) > 0L) {
    mapped_compound_ids[index] <-
      matches[1]
  }
}

compound_rank_data$a18d_compound_id <-
  mapped_compound_ids

mapping_table <- data.frame(
  graph_molecule_id =
    compound_rank_data$molecule_chembl_id,

  graph_display_name =
    compound_rank_data$name,

  graph_node_id =
    paste0(
      "compound::",
      compound_rank_data$molecule_chembl_id
    ),

  a18d_compound_id =
    compound_rank_data$a18d_compound_id,

  mapping_status =
    ifelse(
      is.na(
        compound_rank_data$a18d_compound_id
      ),
      "not_mapped",
      "mapped"
    ),

  stringsAsFactors = FALSE
)

# ============================================================
# 6. Create product HTML
# ============================================================

cat(
  "Preparing commercial-product detail HTML...\n"
)

build_product_html <- function(compound_id) {
  current_products <- products[
    products$compound_id ==
      compound_id,
    ,
    drop = FALSE
  ]

  if (nrow(current_products) == 0L) {
    return(
      paste0(
        "<p class='product-empty'>",
        "No public commercial-product record was available.",
        "</p>"
      )
    )
  }

  product_html <- character(
    nrow(current_products)
  )

  for (index in seq_len(nrow(current_products))) {
    heading <- display_text(
      current_products$vendor_heading[index],
      "Public commercial listing"
    )

    description <- display_text(
      dplyr::coalesce(
        clean_text(
          current_products$
            information_text[index]
        ),
        clean_text(
          current_products$
            information_name[index]
        )
      ),
      "No public listing description"
    )

    identity_class <- display_text(
      current_products$
        identity_match_class[index],
      "Manual identity review required"
    )

    verification <- display_text(
      current_products$
        manual_verification_status[index],
      "Not verified"
    )

    product_url <- clean_text(
      current_products$product_url[index]
    )

    product_link <- ""

    if (
      length(product_url) > 0L &&
      !is.na(product_url)
    ) {
      first_url <- trimws(
        strsplit(
          product_url,
          ";",
          fixed = TRUE
        )[[1]][1]
      )

      product_link <- paste0(
        "<br><a href='",
        first_url,
        "' target='_blank' rel='noopener noreferrer'>",
        "Open public product information",
        "</a>"
      )
    }

    product_html[index] <- paste0(
      "<div class='product-card'>",
      "<strong>",
      heading,
      "</strong><br>",
      description,
      "<br><small>Identity classification: ",
      identity_class,
      "</small><br>",
      "<small>Verification: ",
      verification,
      "</small>",
      product_link,
      "</div>"
    )
  }

  paste(
    product_html,
    collapse = ""
  )
}

# ============================================================
# 7. Create compound detail HTML
# ============================================================

cat(
  "Embedding detail HTML into compound nodes...\n"
)

build_detail_html <- function(
  compound_id
) {
  record <- compound_details[
    compound_details$compound_id ==
      compound_id,
    ,
    drop = FALSE
  ]

  if (nrow(record) == 0L) {
    return(
      paste0(
        "<h2>Compound details unavailable</h2>",
        "<p>No A18D record matched this compound node.</p>"
      )
    )
  }

  record <- record[
    1,
    ,
    drop = FALSE
  ]

  identity_card <- html_card(
    "Compound identity",

    paste0(
      html_row(
        "Compound ID",
        record$compound_id
      ),

      html_row(
        "Preferred name",
        record$compound_name
      ),

      html_row(
        "Compound class",
        record$compound_class
      ),

      html_row(
        "ChEMBL ID",
        record$chembl_id
      ),

      html_row(
        "PubChem CID",
        record$pubchem_cid
      ),

      html_row(
        "Full InChIKey",
        record$full_inchi_key
      ),

      html_row(
        "Parent key",
        record$parent_connectivity_key
      ),

      html_row(
        "Canonical SMILES",
        record$canonical_smiles
      )
    )
  )

  evidence_card <- html_card(
    "Scientific evidence",

    paste0(
      html_row(
        "Classification",
        record$biological_classification
      ),

      html_row(
        "Evidence status",
        record$scientific_evidence_status
      ),

      html_row(
        "Progression status",
        record$progression_status
      ),

      html_row(
        "Matched proteins",
        record$matched_network_proteins
      ),

      html_row(
        "Independent assays",
        record$independent_assay_count
      ),

      html_row(
        "Independent documents",
        record$independent_document_count
      )
    )
  )

  similarity_card <- html_card(
    "Chemical and phenotypic neighbourhood",

    paste0(
      html_row(
        "Nearest RESKO query",
        record$nearest_query_name
      ),

      html_row(
        "Tanimoto similarity",
        record$nearest_query_tanimoto
      ),

      html_row(
        "Similarity band",
        record$similarity_band
      ),

      html_row(
        "Closest side-effect neighbour",
        record$closest_side_effect_neighbour
      ),

      html_row(
        "Side-effect Jaccard",
        record$highest_side_effect_jaccard
      ),

      html_row(
        "Side-effect cosine",
        record$highest_side_effect_cosine
      )
    )
  )

  clinical_card <- html_card(
    "Clinical and safety information",

    paste0(
      html_row(
        "Maximum clinical phase",
        record$maximum_clinical_phase
      ),

      html_row(
        "First approval year",
        record$first_approval_year
      ),

      html_row(
        "SIDER representation",
        record$sider_representation
      ),

      html_row(
        "Side-effect records",
        record$side_effect_record_count
      ),

      html_row(
        "Side-effect terms",
        record$side_effect_names
      ),

      html_row(
        "Indication records",
        record$indication_record_count
      ),

      html_row(
        "Indications",
        record$indication_names
      )
    )
  )

  commercial_card <- html_card(
    "Commercial availability",

    paste0(
      html_row(
        "Availability status",
        record$commercial_availability_status
      ),

      html_row(
        "Supplier count",
        record$supplier_count
      ),

      html_row(
        "Commercial records",
        record$commercial_product_record_count
      ),

      html_row(
        "Supplier names",
        record$supplier_names
      ),

      html_row(
        "Identity classes",
        record$listing_identity_classes
      ),

      html_row(
        "Exact identity confirmed",
        record$exact_product_identity_confirmed
      ),

      html_row(
        "Manual review required",
        record$manual_product_review_required
      ),

      html_row(
        "Procurement readiness",
        record$procurement_readiness
      ),

      html_row(
        "Commercial data checked",
        record$commercial_data_checked
      )
    )
  )

  product_card <- html_card(
    "Public commercial-product records",
    build_product_html(compound_id)
  )

  sider_note <- ""

  if (
    clean_text(
      record$sider_representation
    ) ==
    "not_represented_in_sider_4.1"
  ) {
    sider_note <- paste0(
      "<div class='detail-note'>",
      "<strong>SIDER interpretation:</strong> ",
      "Not represented in SIDER 4.1 does not mean that the ",
      "compound has no side effects.",
      "</div>"
    )
  }

  warning <- paste0(
    "<div class='detail-warning'>",
    "<strong>Procurement warning:</strong> ",
    "Public vendor information does not confirm exact identity, ",
    "stereochemistry, chemical form, purity, stock, price, pack size, ",
    "lead time, or suitability for biological testing.",
    "</div>"
  )

  paste0(
    "<div class='compound-detail-header'>",
    "<h2>",
    display_text(
      record$compound_name,
      compound_id
    ),
    "</h2>",
    "<div class='detail-subtitle'>",
    display_text(compound_id),
    "</div>",
    "</div>",

    "<div class='detail-grid'>",
    identity_card,
    evidence_card,
    similarity_card,
    clinical_card,
    commercial_card,
    product_card,
    "</div>",

    warning,
    sider_note
  )
}

detail_html <- character(
  nrow(compound_rank_data)
)

for (
  index in seq_len(
    nrow(compound_rank_data)
  )
) {
  compound_id <-
    compound_rank_data$
      a18d_compound_id[index]

  if (
    is.na(compound_id) ||
    compound_id == ""
  ) {
    detail_html[index] <- paste0(
      "<h2>Compound mapping unavailable</h2>",
      "<p>The selected graph compound could not be mapped to ",
      "an A18D compound-detail record.</p>"
    )
  } else {
    detail_html[index] <- build_detail_html(
      compound_id
    )
  }
}

compound_rank_data$detail_html <-
  detail_html

# ============================================================
# 8. Calculate biological node metrics
# ============================================================

cat(
  "Calculating biological-node metrics...\n"
)

ppi$score <- safe_numeric(
  ppi$score
)

ppi$score[
  is.na(ppi$score)
] <- 0

ppi_degree <- table(
  c(
    ppi$source,
    ppi$target
  )
)

ppi_degree <- data.frame(
  protein =
    names(ppi_degree),

  ppi_degree =
    as.numeric(ppi_degree),

  stringsAsFactors = FALSE
)

pathway_degree <- protein_pathway |>
  dplyr::count(
    .data$protein,
    name = "pathway_degree"
  )

go_degree <- protein_go |>
  dplyr::count(
    .data$protein,
    name = "go_degree"
  )

compound_degree <- compound_edges_data |>
  dplyr::count(
    .data$queried_protein,
    name = "compound_degree"
  ) |>
  dplyr::rename(
    protein =
      "queried_protein"
  )

protein_metrics <- proteins |>
  dplyr::left_join(
    ppi_degree,
    by = "protein"
  ) |>
  dplyr::left_join(
    pathway_degree,
    by = "protein"
  ) |>
  dplyr::left_join(
    go_degree,
    by = "protein"
  ) |>
  dplyr::left_join(
    compound_degree,
    by = "protein"
  ) |>
  dplyr::mutate(
    dplyr::across(
      c(
        "ppi_degree",
        "pathway_degree",
        "go_degree",
        "compound_degree"
      ),
      function(value) {
        dplyr::coalesce(
          value,
          0
        )
      }
    ),

    node_size =
      18 +
      3 *
      sqrt(
        .data$ppi_degree +
        .data$pathway_degree +
        .data$go_degree +
        .data$compound_degree
      )
  )

pathways <- pathways |>
  dplyr::left_join(
    protein_pathway |>
      dplyr::count(
        .data$pathway,
        name = "protein_count"
      ),
    by = "pathway"
  ) |>
  dplyr::mutate(
    protein_count =
      dplyr::coalesce(
        .data$protein_count,
        0L
      ),

    node_size =
      10 +
      2.3 *
      sqrt(
        .data$protein_count
      )
  )

go_processes <- go_processes |>
  dplyr::left_join(
    protein_go |>
      dplyr::count(
        .data$biological_process,
        name = "protein_count"
      ),
    by = c(
      "go_term" =
        "biological_process"
    )
  ) |>
  dplyr::mutate(
    protein_count =
      dplyr::coalesce(
        .data$protein_count,
        0L
      ),

    node_size =
      9 +
      2.1 *
      sqrt(
        .data$protein_count
      )
  )

# ============================================================
# 9. Build nodes
# ============================================================

cat(
  "Building network nodes...\n"
)

protein_nodes <- protein_metrics |>
  dplyr::transmute(
    id =
      paste0(
        "protein::",
        .data$protein
      ),

    label =
      .data$protein,

    group =
      "Protein",

    shape =
      "dot",

    size =
      .data$node_size,

    title =
      paste0(
        "<b>Protein:</b> ",
        .data$protein,
        "<br><b>Compound links:</b> ",
        .data$compound_degree
      ),

    detail_html =
      paste0(
        "<h2>Protein selected</h2>",
        "<p><strong>Protein:</strong> ",
        escape_html(.data$protein),
        "</p>",
        "<p>Protein details remain available in the graph tooltip.</p>"
      )
  )

pathway_nodes <- pathways |>
  dplyr::transmute(
    id =
      paste0(
        "pathway::",
        .data$pathway
      ),

    label =
      .data$pathway,

    group =
      "Pathway",

    shape =
      "dot",

    size =
      .data$node_size,

    title =
      paste0(
        "<b>Reactome pathway:</b> ",
        .data$pathway
      ),

    detail_html =
      paste0(
        "<h2>Pathway selected</h2>",
        "<p><strong>Pathway:</strong> ",
        escape_html(.data$pathway),
        "</p>"
      )
  )

go_nodes <- go_processes |>
  dplyr::transmute(
    id =
      paste0(
        "go::",
        .data$go_term
      ),

    label =
      .data$go_term,

    group =
      "BiologicalProcess",

    shape =
      "dot",

    size =
      .data$node_size,

    title =
      paste0(
        "<b>GO biological process:</b> ",
        .data$go_term
      ),

    detail_html =
      paste0(
        "<h2>Biological process selected</h2>",
        "<p><strong>GO process:</strong> ",
        escape_html(.data$go_term),
        "</p>"
      )
  )

compound_nodes <- compound_rank_data |>
  dplyr::transmute(
    id =
      paste0(
        "compound::",
        .data$molecule_chembl_id
      ),

    label =
      .data$rank_label,

    group =
      "Compound",

    shape =
      "diamond",

    size =
      .data$compound_size,

    borderWidth =
      3,

    title =
      paste0(
        "<b>Compound:</b> ",
        .data$name,
        "<br><b>ChEMBL ID:</b> ",
        .data$molecule_chembl_id,
        "<br><b>Direct rank:</b> ",
        .data$direct_rank,
        "<br><b>Network rank:</b> ",
        .data$network_rank,
        "<br><br><b>Click to show full details below.</b>"
      ),

    detail_html =
      .data$detail_html
  )

nodes <- dplyr::bind_rows(
  protein_nodes,
  pathway_nodes,
  go_nodes,
  compound_nodes
)

# ============================================================
# 10. Build edges
# ============================================================

cat(
  "Building network edges...\n"
)

ppi_edges <- ppi |>
  dplyr::transmute(
    from =
      paste0(
        "protein::",
        .data$source
      ),

    to =
      paste0(
        "protein::",
        .data$target
      ),

    color =
      "#9A9A9A",

    width =
      0.9 +
      2.7 *
      .data$score,

    arrows =
      ""
  )

pathway_edges <- protein_pathway |>
  dplyr::transmute(
    from =
      paste0(
        "protein::",
        .data$protein
      ),

    to =
      paste0(
        "pathway::",
        .data$pathway
      ),

    color =
      "#40916C",

    width =
      1.15,

    arrows =
      ""
  )

go_edges <- protein_go |>
  dplyr::transmute(
    from =
      paste0(
        "protein::",
        .data$protein
      ),

    to =
      paste0(
        "go::",
        .data$biological_process
      ),

    color =
      "#5B8FD1",

    width =
      1.05,

    arrows =
      ""
  )

compound_edges <- compound_edges_data |>
  dplyr::mutate(
    edge_color =
      dplyr::case_when(
        .data$relationship ==
          "HAS_INHIBITORY_ACTIVITY_AGAINST" ~
          "#C43C39",

        .data$relationship ==
          "BINDS_TO" ~
          "#8E63B0",

        TRUE ~
          "#C77729"
      )
  ) |>
  dplyr::transmute(
    from =
      paste0(
        "compound::",
        .data$molecule_chembl_id
      ),

    to =
      paste0(
        "protein::",
        .data$queried_protein
      ),

    color =
      .data$edge_color,

    width =
      2.4,

    arrows =
      "",

    title =
      paste0(
        "<b>",
        .data$relationship,
        "</b><br>",
        .data$display_name,
        " - ",
        .data$queried_protein
      )
  )

edges <- dplyr::bind_rows(
  ppi_edges,
  pathway_edges,
  go_edges,
  compound_edges
)

# ============================================================
# 11. Build network and click handler
# ============================================================

cat(
  "Building visNetwork widget...\n"
)

click_javascript <- paste0(
  "function(params){",
  "if(!params.nodes || params.nodes.length===0){return;}",
  "var nodeId=params.nodes[0];",
  "var node=this.body.data.nodes.get(nodeId);",
  "var panel=document.getElementById('resko-compound-detail-panel');",
  "if(!panel){return;}",
  "if(node && node.detail_html){",
  "panel.innerHTML=node.detail_html;",
  "panel.scrollIntoView({behavior:'smooth',block:'start'});",
  "}else{",
  "panel.innerHTML='<h2>Details unavailable</h2>",
  "<p>No detail record was attached to this node.</p>';",
  "}",
  "}"
)

stabilisation_javascript <- paste0(
  "function(){",
  "this.setOptions({physics:{enabled:false}});",
  "this.fit();",
  "}"
)

network <- visNetwork::visNetwork(
  nodes,
  edges,
  width = "100%",
  height = "900px",
  main =
    "Ranked eEF1A Compound Network with Compound Details"
) |>
  visNetwork::visGroups(
    groupname = "Compound",
    shape = "diamond",
    color = list(
      background = "#E58B3A",
      border = "#713900",
      highlight = list(
        background = "#5B9CFF",
        border = "#1769D2"
      )
    ),
    font = list(
      face = "Calibri",
      size = 14
    )
  ) |>
  visNetwork::visGroups(
    groupname = "Protein",
    shape = "dot",
    color = list(
      background = "#D95F59",
      border = "#8E2F2A"
    )
  ) |>
  visNetwork::visGroups(
    groupname = "Pathway",
    shape = "dot",
    color = list(
      background = "#4E9F6E",
      border = "#25633A"
    )
  ) |>
  visNetwork::visGroups(
    groupname = "BiologicalProcess",
    shape = "dot",
    color = list(
      background = "#5B8FD1",
      border = "#2B5797"
    )
  ) |>
  visNetwork::visNodes(
    font = list(
      face = "Calibri",
      strokeWidth = 3,
      strokeColor = "#FFFFFF"
    )
  ) |>
  visNetwork::visEdges(
    arrows = "",
    smooth = FALSE,
    selectionWidth = 3,
    hoverWidth = 2
  ) |>
  visNetwork::visOptions(
    highlightNearest = list(
      enabled = TRUE,
      degree = 1,
      hover = FALSE,
      algorithm = "all"
    ),

    nodesIdSelection = list(
      enabled = TRUE,
      useLabels = TRUE,
      main = "Find a node"
    ),

    selectedBy = list(
      variable = "group",
      main = "Filter by node type",
      multiple = TRUE
    )
  ) |>
  visNetwork::visInteraction(
    hover = TRUE,
    navigationButtons = TRUE,
    keyboard = TRUE,
    hideEdgesOnDrag = TRUE
  ) |>
  visNetwork::visPhysics(
    solver = "forceAtlas2Based",

    forceAtlas2Based = list(
      gravitationalConstant = -110,
      centralGravity = 0.006,
      springLength = 265,
      springConstant = 0.033,
      damping = 0.65,
      avoidOverlap = 0.8
    ),

    stabilization = list(
      enabled = TRUE,
      iterations = 3500,
      fit = TRUE
    )
  ) |>
  visNetwork::visLayout(
    randomSeed = 42,
    improvedLayout = TRUE
  ) |>
  visNetwork::visEvents(
    click =
      click_javascript,

    stabilizationIterationsDone =
      stabilisation_javascript
  )

# ============================================================
# 12. Detail panel and styling
# ============================================================

detail_panel <- htmltools::tags$div(
  id = "resko-compound-detail-panel",
  class = "compound-detail-panel",

  htmltools::tags$h2(
    "Selected compound details"
  ),

  htmltools::tags$p(
    "Click an orange compound node to display compound identity, ",
    "scientific evidence, side effects, and commercial availability."
  )
)

page_style <- htmltools::tags$style(
  htmltools::HTML(
    paste0(
      "body{font-family:Calibri,Arial,sans-serif;",
      "margin:0;padding:14px;color:#222;}",

      ".vis-network{border:1px solid #ddd;",
      "border-radius:6px;background:#fff;}",

      ".compound-detail-panel{margin-top:18px;padding:18px;",
      "border:1px solid #cdd5df;border-radius:8px;",
      "background:#f7f9fc;min-height:170px;}",

      ".compound-detail-header h2{margin-bottom:4px;",
      "color:#713900;}",

      ".detail-subtitle{color:#596675;",
      "font-size:13px;margin-bottom:14px;}",

      ".detail-grid{display:grid;",
      "grid-template-columns:repeat(auto-fit,minmax(290px,1fr));",
      "gap:14px;}",

      ".detail-card{background:#fff;",
      "border:1px solid #d8dee6;",
      "padding:13px;border-radius:6px;",
      "overflow-wrap:anywhere;}",

      ".detail-card h3{margin:0 0 10px 0;",
      "color:#35475a;}",

      ".detail-row{display:grid;",
      "grid-template-columns:145px 1fr;",
      "gap:8px;padding:5px 0;",
      "border-bottom:1px solid #edf0f4;}",

      ".detail-row:last-child{border-bottom:none;}",

      ".detail-label{font-weight:600;color:#566473;}",

      ".detail-value{overflow-wrap:anywhere;}",

      ".detail-warning{margin-top:14px;padding:12px;",
      "background:#fff4e5;",
      "border-left:5px solid #e67e22;}",

      ".detail-note{margin-top:14px;padding:12px;",
      "background:#eaf4ff;",
      "border-left:5px solid #2176ae;}",

      ".product-card{margin-top:8px;padding:9px;",
      "border:1px solid #d5dbe3;",
      "border-radius:5px;background:#fff;}",

      ".product-card a{color:#125fa5;",
      "text-decoration:none;}",

      ".product-card a:hover{text-decoration:underline;}",

      ".product-empty{color:#647383;font-style:italic;}",

      "@media(max-width:700px){",
      ".detail-row{grid-template-columns:1fr;}",
      "}"
    )
  )
)

network <- htmlwidgets::appendContent(
  network,
  detail_panel
)

network <- htmlwidgets::prependContent(
  network,
  page_style
)

# ============================================================
# 13. Save network
# ============================================================

cat(
  "Saving A18E live network...\n"
)

backup_file(
  outputs[["network_html"]]
)

htmlwidgets::saveWidget(
  network,
  outputs[["network_html"]],
  selfcontained = TRUE
)

if (
  !file.exists(
    outputs[["network_html"]]
  ) ||
  file.info(
    outputs[["network_html"]]
  )$size <= 0L
) {
  stop(
    "The A18E live-network HTML was not created.",
    call. = FALSE
  )
}

# ============================================================
# 14. Build validation outputs
# ============================================================

mapped_count <- sum(
  mapping_table$mapping_status ==
    "mapped"
)

unmapped_count <- sum(
  mapping_table$mapping_status ==
    "not_mapped"
)

detail_html_count <- sum(
  !is.na(
    compound_nodes$detail_html
  ) &
  compound_nodes$detail_html != ""
)

validation_table <- data.frame(
  validation_item = c(
    "ranked_compound_node_count",
    "ranked_compounds_mapped",
    "ranked_compounds_unmapped",
    "compound_nodes_with_detail_html",
    "compound_detail_record_count",
    "commercial_product_record_count",
    "detail_panel_present",
    "visEvents_click_handler_present",
    "HTML_output_created"
  ),

  observed_value = c(
    as.character(
      nrow(compound_nodes)
    ),

    as.character(
      mapped_count
    ),

    as.character(
      unmapped_count
    ),

    as.character(
      detail_html_count
    ),

    as.character(
      nrow(compound_details)
    ),

    as.character(
      nrow(products)
    ),

    "TRUE",
    "TRUE",
    "TRUE"
  ),

  validation_status = c(
    ifelse(
      nrow(compound_nodes) > 0L,
      "passed",
      "failed"
    ),

    ifelse(
      mapped_count > 0L,
      "passed",
      "failed"
    ),

    ifelse(
      unmapped_count == 0L,
      "passed",
      "review_required"
    ),

    ifelse(
      detail_html_count ==
        nrow(compound_nodes),
      "passed",
      "failed"
    ),

    ifelse(
      nrow(compound_details) == 16L,
      "passed",
      "failed"
    ),

    "passed",
    "passed",
    "passed",
    "passed"
  ),

  stringsAsFactors = FALSE
)

summary_table <- data.frame(
  metric = c(
    "network_node_count",
    "network_edge_count",
    "ranked_compound_node_count",
    "ranked_compounds_mapped",
    "ranked_compounds_unmapped",
    "compound_nodes_with_embedded_details",
    "compound_detail_record_count",
    "commercial_product_record_count",
    "compounds_with_public_vendor_information",
    "detail_panel_enabled",
    "output_html_created"
  ),

  value = c(
    as.character(
      nrow(nodes)
    ),

    as.character(
      nrow(edges)
    ),

    as.character(
      nrow(compound_nodes)
    ),

    as.character(
      mapped_count
    ),

    as.character(
      unmapped_count
    ),

    as.character(
      detail_html_count
    ),

    as.character(
      nrow(compound_details)
    ),

    as.character(
      nrow(products)
    ),

    as.character(
      sum(
        safe_logical(
          compound_details$
            public_vendor_information_present
        ),
        na.rm = TRUE
      )
    ),

    "TRUE",
    "TRUE"
  ),

  stringsAsFactors = FALSE
)

safe_write_csv(
  mapping_table,
  outputs[["mapping"]]
)

safe_write_csv(
  validation_table,
  outputs[["validation"]]
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

# ============================================================
# 15. Validate final HTML
# ============================================================

cat(
  "Validating A18E outputs...\n"
)

for (path in outputs) {
  if (
    !file.exists(path) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "A18E output validation failed: ",
      path,
      call. = FALSE
    )
  }
}

html_text <- paste(
  readLines(
    outputs[["network_html"]],
    warn = FALSE
  ),
  collapse = "\n"
)

required_markers <- c(
  "resko-compound-detail-panel",
  "detail_html",
  "this.body.data.nodes.get",
  "Selected compound details",
  "Procurement warning"
)

missing_markers <- required_markers[
  !vapply(
    required_markers,
    function(marker) {
      grepl(
        marker,
        html_text,
        fixed = TRUE
      )
    },
    FUN.VALUE = logical(1)
  )
]

if (length(missing_markers) > 0L) {
  stop(
    "A18E HTML is missing marker(s): ",
    paste(
      missing_markers,
      collapse = ", "
    ),
    call. = FALSE
  )
}

cat(
  "All A18E output validations passed.\n"
)

# ============================================================
# 16. Completion
# ============================================================

cat(
  "\nA18E live network created successfully.\n"
)

cat(
  "Network nodes: ",
  nrow(nodes),
  "\n",
  sep = ""
)

cat(
  "Network edges: ",
  nrow(edges),
  "\n",
  sep = ""
)

cat(
  "Ranked compound nodes: ",
  nrow(compound_nodes),
  "\n",
  sep = ""
)

cat(
  "Ranked compounds mapped to A18D: ",
  mapped_count,
  "\n",
  sep = ""
)

cat(
  "Compound nodes with embedded details: ",
  detail_html_count,
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
  "Output: ",
  outputs[["network_html"]],
  "\n",
  sep = ""
)

utils::browseURL(
  normalizePath(
    outputs[["network_html"]],
    winslash = "/",
    mustWork = TRUE
  )
)