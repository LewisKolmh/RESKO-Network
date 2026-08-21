#!/usr/bin/env Rscript

# ============================================================
# RESKO A18A: Build compound detail and supplier manifests
# ============================================================
#
# Run from the RESKO project root:
#
# Rscript scripts/A18A_build_compound_detail_manifest.R
#
# Required inputs:
#
# results/A16_structures_combined.csv
# results/A17C_sider_similarity_candidates.csv
# results/A17E_candidate_classification_corrected.csv
#
# Optional inputs:
#
# results/A17D_candidate_annotations.csv
# results/A17D_candidate_side_effects.csv
# results/A17D_candidate_indications.csv
# results/A17D_side_effect_similarity.csv
# results/A17E_provenance_corrected_evidence.csv
#
# Outputs:
#
# results/A18A_compound_detail_manifest.csv
# results/A18A_compound_identifiers.csv
# results/A18A_compound_evidence_summary.csv
# results/A18A_supplier_lookup_manifest.csv
# results/A18A_data_completeness_review.csv
# results/A18A_summary.csv
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
  "tibble"
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
# 2. Define paths
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

if (!dir.exists(results_dir)) {
  stop(
    "The results directory was not found. ",
    "Run A18A from the RESKO project root.",
    call. = FALSE
  )
}

if (!dir.exists(scripts_dir)) {
  stop(
    "The scripts directory was not found. ",
    "Run A18A from the RESKO project root.",
    call. = FALSE
  )
}

input_a16 <- file.path(
  results_dir,
  "A16_structures_combined.csv"
)

input_a17c <- file.path(
  results_dir,
  "A17C_sider_similarity_candidates.csv"
)

input_a17e_classification <- file.path(
  results_dir,
  "A17E_candidate_classification_corrected.csv"
)

input_a17d_annotations <- file.path(
  results_dir,
  "A17D_candidate_annotations.csv"
)

input_a17d_side_effects <- file.path(
  results_dir,
  "A17D_candidate_side_effects.csv"
)

input_a17d_indications <- file.path(
  results_dir,
  "A17D_candidate_indications.csv"
)

input_a17d_side_effect_similarity <- file.path(
  results_dir,
  "A17D_side_effect_similarity.csv"
)

input_a17e_provenance <- file.path(
  results_dir,
  "A17E_provenance_corrected_evidence.csv"
)

outputs <- c(
  detail_manifest = file.path(
    results_dir,
    "A18A_compound_detail_manifest.csv"
  ),
  identifiers = file.path(
    results_dir,
    "A18A_compound_identifiers.csv"
  ),
  evidence = file.path(
    results_dir,
    "A18A_compound_evidence_summary.csv"
  ),
  supplier_lookup = file.path(
    results_dir,
    "A18A_supplier_lookup_manifest.csv"
  ),
  completeness = file.path(
    results_dir,
    "A18A_data_completeness_review.csv"
  ),
  summary = file.path(
    results_dir,
    "A18A_summary.csv"
  )
)

required_inputs <- c(
  input_a16,
  input_a17c,
  input_a17e_classification
)

for (path in required_inputs) {
  if (
    !file.exists(path) ||
    is.na(file.info(path)$size) ||
    file.info(path)$size <= 0L
  ) {
    stop(
      "Required input is missing or empty: ",
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

connectivity_key <- function(inchi_key) {
  key <- clean_text(inchi_key)

  ifelse(
    !is.na(key) &
      nchar(key) >= 14L,
    substr(
      key,
      1L,
      14L
    ),
    NA_character_
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

get_column <- function(
  data,
  possible_names,
  default = NA_character_
) {
  available <- possible_names[
    possible_names %in% names(data)
  ]

  if (length(available) == 0L) {
    return(
      rep(
        default,
        nrow(data)
      )
    )
  }

  data[[
    available[1]
  ]]
}

read_required_csv <- function(path) {
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
      tibble::tibble()
    )
  }

  readr::read_csv(
    path,
    col_types = readr::cols(
      .default = readr::col_character()
    ),
    progress = FALSE
  )
}

backup_existing_file <- function(path) {
  if (!file.exists(path)) {
    return(
      invisible(NULL)
    )
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
    from = path,
    to = backup_path,
    overwrite = FALSE
  )

  if (!copied) {
    stop(
      "Could not preserve previous output: ",
      path,
      call. = FALSE
    )
  }

  if (
    !file.exists(backup_path) ||
    is.na(file.info(backup_path)$size) ||
    file.info(backup_path)$size <= 0L
  ) {
    stop(
      "Backup verification failed: ",
      backup_path,
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
      "Refusing to write list-column(s): ",
      paste(
        list_columns,
        collapse = ", "
      ),
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

# ============================================================
# 4. Read input data
# ============================================================

cat(
  "Reading A18A inputs...\n"
)

a16 <- read_required_csv(
  input_a16
)

a17c <- read_required_csv(
  input_a17c
)

a17e_classification <- read_required_csv(
  input_a17e_classification
)

a17d_annotations <- read_optional_csv(
  input_a17d_annotations
)

a17d_side_effects <- read_optional_csv(
  input_a17d_side_effects
)

a17d_indications <- read_optional_csv(
  input_a17d_indications
)

a17d_side_effect_similarity <- read_optional_csv(
  input_a17d_side_effect_similarity
)

a17e_provenance <- read_optional_csv(
  input_a17e_provenance
)

cat(
  "Required and optional inputs loaded.\n"
)

# ============================================================
# 5. Validate essential input columns
# ============================================================

cat(
  "Validating input columns and row counts...\n"
)

required_a16_columns <- c(
  "compound_id",
  "compound_name",
  "compound_class",
  "pubchem_cid",
  "canonical_smiles",
  "analysis_smiles",
  "inchi_key"
)

required_a17c_columns <- c(
  "sider_stitch_flat_id",
  "sider_drug_name",
  "sider_pubchem_cid",
  "sider_inchi_key",
  "nearest_query_name",
  "nearest_query_tanimoto",
  "similarity_band",
  "unique_side_effect_count",
  "atc_codes"
)

required_a17e_columns <- c(
  "candidate_name",
  "corrected_biological_classification",
  "corrected_progression_status"
)

missing_a16_columns <- setdiff(
  required_a16_columns,
  names(a16)
)

missing_a17c_columns <- setdiff(
  required_a17c_columns,
  names(a17c)
)

missing_a17e_columns <- setdiff(
  required_a17e_columns,
  names(a17e_classification)
)

if (length(missing_a16_columns) > 0L) {
  stop(
    "A16 input is missing column(s): ",
    paste(
      missing_a16_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(missing_a17c_columns) > 0L) {
  stop(
    "A17C input is missing column(s): ",
    paste(
      missing_a17c_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (length(missing_a17e_columns) > 0L) {
  stop(
    "A17E classification is missing column(s): ",
    paste(
      missing_a17e_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (nrow(a16) != 10L) {
  stop(
    "Expected 10 A16 compounds, found ",
    nrow(a16),
    ".",
    call. = FALSE
  )
}

if (nrow(a17c) != 6L) {
  stop(
    "Expected 6 A17C comparators, found ",
    nrow(a17c),
    ".",
    call. = FALSE
  )
}

if (nrow(a17e_classification) != 6L) {
  stop(
    "Expected 6 A17E classifications, found ",
    nrow(a17e_classification),
    ".",
    call. = FALSE
  )
}

cat(
  "Input validation completed.\n"
)

# ============================================================
# 6. Prepare the ten A16 compounds
# ============================================================

cat(
  "Preparing the 10 A16 compound records...\n"
)

a16_compound_id <- clean_text(
  a16$compound_id
)

a16_compound_name <- clean_text(
  a16$compound_name
)

a16_compound_class <- clean_text(
  a16$compound_class
)

a16_inchi_key <- clean_text(
  a16$inchi_key
)

a16_compounds <- tibble::tibble(
  compound_id =
    a16_compound_id,

  compound_name =
    a16_compound_name,

  normalised_compound_name =
    normalise_name(
      a16_compound_name
    ),

  compound_class =
    a16_compound_class,

  compound_origin =
    "A16_current_or_reference_compound",

  chembl_id =
    ifelse(
      grepl(
        "^CHEMBL[0-9]+$",
        a16_compound_id
      ),
      a16_compound_id,
      NA_character_
    ),

  pubchem_cid =
    safe_numeric(
      a16$pubchem_cid
    ),

  sider_stitch_flat_id =
    NA_character_,

  full_inchi_key =
    a16_inchi_key,

  parent_connectivity_key =
    connectivity_key(
      a16_inchi_key
    ),

  canonical_smiles =
    clean_text(
      a16$canonical_smiles
    ),

  analysis_smiles =
    clean_text(
      a16$analysis_smiles
    ),

  atc_codes =
    NA_character_,

  nearest_query_name =
    a16_compound_name,

  nearest_query_tanimoto =
    1,

  similarity_band =
    "self_or_reference_query",

  biological_classification =
    ifelse(
      grepl(
        "reference",
        tolower(
          dplyr::coalesce(
            a16_compound_class,
            ""
          )
        )
      ),
      "Established eEF1A reference ligand",
      "Existing RESKO candidate"
    ),

  progression_status =
    "existing_resko_compound",

  sider_representation =
    "not_represented_in_sider_4.1",

  unique_side_effect_count =
    0,

  indication_count =
    0,

  indication_names =
    NA_character_,

  independent_assay_count =
    NA_real_,

  independent_document_count =
    NA_real_,

  matched_network_proteins =
    NA_character_,

  max_phase =
    NA_real_,

  first_approval =
    NA_real_,

  data_source =
    "A16"
)

if (nrow(a16_compounds) != 10L) {
  stop(
    "A16 preparation did not produce 10 rows.",
    call. = FALSE
  )
}

cat(
  "A16 compound records completed: 10 rows.\n"
)

# ============================================================
# 7. Prepare A17 classifications
# ============================================================

cat(
  "Preparing A17 classification data...\n"
)

classification_independent_assays <- safe_numeric(
  get_column(
    a17e_classification,
    c("independent_assay_count"),
    default = "0"
  )
)

classification_independent_documents <- safe_numeric(
  get_column(
    a17e_classification,
    c("independent_document_count"),
    default = "0"
  )
)

classification_network_proteins <- clean_text(
  get_column(
    a17e_classification,
    c("matched_network_proteins"),
    default = NA_character_
  )
)

classification_prepared <- tibble::tibble(
  normalised_compound_name =
    normalise_name(
      a17e_classification$candidate_name
    ),

  biological_classification =
    clean_text(
      a17e_classification$
        corrected_biological_classification
    ),

  progression_status =
    clean_text(
      a17e_classification$
        corrected_progression_status
    ),

  independent_assay_count =
    classification_independent_assays,

  independent_document_count =
    classification_independent_documents,

  matched_network_proteins =
    classification_network_proteins
)

classification_prepared <- classification_prepared |>
  dplyr::distinct(
    .data$normalised_compound_name,
    .keep_all = TRUE
  )

if (nrow(classification_prepared) != 6L) {
  stop(
    "A17 classification preparation did not produce 6 rows.",
    call. = FALSE
  )
}

cat(
  "A17 classification data completed: 6 rows.\n"
)

# ============================================================
# 8. Prepare optional A17 annotations
# ============================================================

cat(
  "Preparing optional A17 annotation data...\n"
)

if (
  nrow(a17d_annotations) > 0L &&
  "candidate_name" %in%
    names(a17d_annotations)
) {
  annotation_chembl_id <- clean_text(
    get_column(
      a17d_annotations,
      c("molecule_chembl_id"),
      default = NA_character_
    )
  )

  annotation_max_phase <- safe_numeric(
    get_column(
      a17d_annotations,
      c("max_phase"),
      default = NA_character_
    )
  )

  annotation_first_approval <- safe_numeric(
    get_column(
      a17d_annotations,
      c("first_approval"),
      default = NA_character_
    )
  )

  annotation_prepared <- tibble::tibble(
    normalised_compound_name =
      normalise_name(
        a17d_annotations$candidate_name
      ),

    chembl_id =
      annotation_chembl_id,

    max_phase =
      annotation_max_phase,

    first_approval =
      annotation_first_approval
  ) |>
    dplyr::distinct(
      .data$normalised_compound_name,
      .keep_all = TRUE
    )
} else {
  annotation_prepared <- tibble::tibble(
    normalised_compound_name =
      character(),
    chembl_id =
      character(),
    max_phase =
      numeric(),
    first_approval =
      numeric()
  )
}

cat(
  "Optional A17 annotation data completed: ",
  nrow(annotation_prepared),
  " rows.\n",
  sep = ""
)

# ============================================================
# 9. Prepare optional indication summary
# ============================================================

cat(
  "Preparing optional indication summary...\n"
)

if (nrow(a17d_indications) > 0L) {
  indication_candidate_name <- clean_text(
    get_column(
      a17d_indications,
      c(
        "candidate_name",
        "sider_drug_name"
      ),
      default = NA_character_
    )
  )

  indication_identifier <- clean_text(
    get_column(
      a17d_indications,
      c(
        "meddra_umls_id",
        "indication_umls_id"
      ),
      default = NA_character_
    )
  )

  indication_name <- clean_text(
    get_column(
      a17d_indications,
      c(
        "indication_name",
        "meddra_concept_name"
      ),
      default = NA_character_
    )
  )

  indication_rows <- tibble::tibble(
    normalised_compound_name =
      normalise_name(
        indication_candidate_name
      ),

    indication_identifier =
      indication_identifier,

    indication_name =
      indication_name
  )

  indication_rows <- indication_rows |>
    dplyr::filter(
      !is.na(
        .data$normalised_compound_name
      )
    )

  indication_split <- split(
    indication_rows,
    indication_rows$normalised_compound_name
  )

  indication_summary_list <- lapply(
    names(indication_split),
    function(current_name) {
      current_data <- indication_split[[
        current_name
      ]]

      unique_identifiers <- unique(
        clean_text(
          current_data$indication_identifier
        )
      )

      unique_identifiers <- unique_identifiers[
        !is.na(unique_identifiers)
      ]

      tibble::tibble(
        normalised_compound_name =
          current_name,

        indication_count =
          length(unique_identifiers),

        indication_names =
          collapse_unique(
            current_data$indication_name
          )
      )
    }
  )

  indication_summary <- dplyr::bind_rows(
    indication_summary_list
  )
} else {
  indication_summary <- tibble::tibble(
    normalised_compound_name =
      character(),
    indication_count =
      numeric(),
    indication_names =
      character()
  )
}

cat(
  "Optional indication summary completed: ",
  nrow(indication_summary),
  " rows.\n",
  sep = ""
)

# ============================================================
# 10. Prepare the six A17 comparator compounds
# ============================================================

cat(
  "Preparing the 6 A17 comparator records...\n"
)

a17_names <- clean_text(
  a17c$sider_drug_name
)

a17_normalised_names <- normalise_name(
  a17_names
)

a17_inchi_keys <- clean_text(
  a17c$sider_inchi_key
)

a17_canonical_smiles <- clean_text(
  get_column(
    a17c,
    c(
      "canonical_smiles",
      "sider_canonical_smiles"
    ),
    default = NA_character_
  )
)

a17_analysis_smiles <- clean_text(
  get_column(
    a17c,
    c(
      "analysis_smiles",
      "canonical_smiles",
      "sider_canonical_smiles"
    ),
    default = NA_character_
  )
)

a17_base <- tibble::tibble(
  compound_id =
    paste0(
      "A17_",
      toupper(
        a17_normalised_names
      )
    ),

  compound_name =
    a17_names,

  normalised_compound_name =
    a17_normalised_names,

  compound_class =
    ifelse(
      a17_normalised_names %in%
        c(
          "nilotinib",
          "imatinib",
          "ponatinib"
        ),
      "kinase_inhibitor_comparator",
      ifelse(
        a17_normalised_names %in%
          c(
            "alprazolam",
            "triazolam",
            "temazepam"
          ),
        "benzodiazepine_related_comparator",
        "A17_chemical_comparator"
      )
    ),

  compound_origin =
    "A17_SIDER_chemical_neighbour",

  pubchem_cid =
    safe_numeric(
      a17c$sider_pubchem_cid
    ),

  sider_stitch_flat_id =
    clean_text(
      a17c$sider_stitch_flat_id
    ),

  full_inchi_key =
    a17_inchi_keys,

  parent_connectivity_key =
    connectivity_key(
      a17_inchi_keys
    ),

  canonical_smiles =
    a17_canonical_smiles,

  analysis_smiles =
    a17_analysis_smiles,

  atc_codes =
    clean_text(
      a17c$atc_codes
    ),

  nearest_query_name =
    clean_text(
      a17c$nearest_query_name
    ),

  nearest_query_tanimoto =
    safe_numeric(
      a17c$nearest_query_tanimoto
    ),

  similarity_band =
    clean_text(
      a17c$similarity_band
    ),

  sider_representation =
    "represented_in_sider_4.1",

  unique_side_effect_count =
    safe_numeric(
      a17c$unique_side_effect_count
    ),

  data_source =
    "A17"
)

a17_compounds <- a17_base |>
  dplyr::left_join(
    classification_prepared,
    by = "normalised_compound_name"
  ) |>
  dplyr::left_join(
    annotation_prepared,
    by = "normalised_compound_name"
  ) |>
  dplyr::left_join(
    indication_summary,
    by = "normalised_compound_name"
  ) |>
  dplyr::mutate(
    biological_classification =
      dplyr::coalesce(
        .data$biological_classification,
        "Chemical and side-effect comparator only"
      ),

    progression_status =
      dplyr::coalesce(
        .data$progression_status,
        "retain_as_chemical_comparator"
      ),

    independent_assay_count =
      dplyr::coalesce(
        .data$independent_assay_count,
        0
      ),

    independent_document_count =
      dplyr::coalesce(
        .data$independent_document_count,
        0
      ),

    indication_count =
      dplyr::coalesce(
        .data$indication_count,
        0
      )
  )

if (nrow(a17_compounds) != 6L) {
  stop(
    "A17 preparation did not produce 6 rows.",
    call. = FALSE
  )
}

cat(
  "A17 comparator records completed: 6 rows.\n"
)

# ============================================================
# 11. Build unified compound-detail manifest
# ============================================================

cat(
  "Building the unified compound-detail manifest...\n"
)

compound_detail_manifest <- dplyr::bind_rows(
  a16_compounds,
  a17_compounds
)

compound_detail_manifest <-
  compound_detail_manifest |>
  dplyr::mutate(
    commercial_availability_status =
      "not_yet_checked",

    exact_commercial_product_present =
      as.logical(NA),

    supplier_count =
      NA_real_,

    commercial_product_count =
      NA_real_,

    commercial_identity_status =
      "not_yet_checked",

    supplier_last_checked =
      NA_character_,

    supplier_manual_review_required =
      TRUE,

    commercial_data_source =
      NA_character_,

    live_network_detail_ready =
      TRUE
  ) |>
  dplyr::arrange(
    .data$compound_origin,
    .data$compound_name
  )

if (nrow(compound_detail_manifest) != 16L) {
  stop(
    "Expected 16 compound-detail records, found ",
    nrow(compound_detail_manifest),
    ".",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    compound_detail_manifest$compound_id
  ) > 0L
) {
  stop(
    "Duplicate compound identifiers were detected.",
    call. = FALSE
  )
}

cat(
  "Unified compound-detail manifest completed: 16 rows.\n"
)

# ============================================================
# 12. Build identifier table
# ============================================================

cat(
  "Building the compound identifier table...\n"
)

required_identifier_columns <- c(
  "compound_id",
  "compound_name",
  "normalised_compound_name",
  "chembl_id",
  "pubchem_cid",
  "sider_stitch_flat_id",
  "full_inchi_key",
  "parent_connectivity_key",
  "canonical_smiles",
  "analysis_smiles"
)

missing_identifier_columns <- setdiff(
  required_identifier_columns,
  names(compound_detail_manifest)
)

if (length(missing_identifier_columns) > 0L) {
  stop(
    "Compound manifest is missing identifier column(s): ",
    paste(
      missing_identifier_columns,
      collapse = ", "
    ),
    call. = FALSE
  )
}

compound_identifiers <-
  compound_detail_manifest[
    ,
    required_identifier_columns,
    drop = FALSE
  ]

compound_identifiers <- tibble::as_tibble(
  compound_identifiers
)

compound_identifiers$
  preferred_supplier_lookup_identifier <-
  ifelse(
    !is.na(
      compound_identifiers$full_inchi_key
    ),
    "full_inchi_key",
    ifelse(
      !is.na(
        compound_identifiers$
          parent_connectivity_key
      ),
      "parent_connectivity_key",
      ifelse(
        !is.na(
          compound_identifiers$pubchem_cid
        ),
        "pubchem_cid",
        ifelse(
          !is.na(
            compound_identifiers$
              canonical_smiles
          ),
          "canonical_smiles",
          "compound_name"
        )
      )
    )
  )

compound_identifiers$
  identifier_review_required <-
  is.na(
    compound_identifiers$full_inchi_key
  )

if (nrow(compound_identifiers) != 16L) {
  stop(
    "Expected 16 identifier rows, found ",
    nrow(compound_identifiers),
    ".",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    compound_identifiers$compound_id
  ) > 0L
) {
  stop(
    "Duplicate compound IDs were found in the identifier table.",
    call. = FALSE
  )
}

cat(
  "Compound identifier table completed: 16 rows.\n"
)

# ============================================================
# 13. Build scientific evidence summary
# ============================================================

cat(
  "Building the scientific evidence summary...\n"
)

evidence_columns <- c(
  "compound_id",
  "compound_name",
  "compound_class",
  "compound_origin",
  "nearest_query_name",
  "nearest_query_tanimoto",
  "similarity_band",
  "biological_classification",
  "progression_status",
  "matched_network_proteins",
  "independent_assay_count",
  "independent_document_count",
  "sider_representation",
  "unique_side_effect_count",
  "indication_count",
  "indication_names",
  "max_phase",
  "first_approval"
)

compound_evidence_summary <-
  compound_detail_manifest[
    ,
    evidence_columns,
    drop = FALSE
  ]

compound_evidence_summary <-
  tibble::as_tibble(
    compound_evidence_summary
  )

biological_text <- dplyr::coalesce(
  compound_evidence_summary$
    biological_classification,
  ""
)

compound_evidence_summary$
  direct_eef1a_evidence <-
  grepl(
    "Direct eEF1A",
    biological_text,
    ignore.case = TRUE
  )

compound_evidence_summary$
  network_evidence <-
  grepl(
    "network",
    biological_text,
    ignore.case = TRUE
  )

compound_evidence_summary$
  chemical_similarity_only <-
  grepl(
    "Chemical",
    biological_text,
    ignore.case = TRUE
  )

compound_evidence_summary$
  scientific_evidence_status <-
  ifelse(
    compound_evidence_summary$
      direct_eef1a_evidence,
    "direct_eef1a_evidence",
    ifelse(
      compound_evidence_summary$
        network_evidence,
      "network_evidence",
      ifelse(
        compound_evidence_summary$
          chemical_similarity_only,
        "chemical_similarity_only",
        ifelse(
          compound_evidence_summary$
            compound_origin ==
            "A16_current_or_reference_compound",
          "existing_resko_evidence",
          "evidence_not_classified"
        )
      )
    )
  )

if (
  nrow(
    compound_evidence_summary
  ) != 16L
) {
  stop(
    "Evidence summary must contain 16 rows.",
    call. = FALSE
  )
}

cat(
  "Scientific evidence summary completed: 16 rows.\n"
)

# ============================================================
# 14. Build supplier lookup manifest
# ============================================================

cat(
  "Building the supplier lookup manifest...\n"
)

lookup_identifier_type <- ifelse(
  !is.na(
    compound_detail_manifest$
      full_inchi_key
  ),
  "full_inchi_key",
  ifelse(
    !is.na(
      compound_detail_manifest$
        parent_connectivity_key
    ),
    "parent_connectivity_key",
    ifelse(
      !is.na(
        compound_detail_manifest$
          pubchem_cid
      ),
      "pubchem_cid",
      ifelse(
        !is.na(
          compound_detail_manifest$
            canonical_smiles
        ),
        "canonical_smiles",
        "compound_name"
      )
    )
  )
)

supplier_lookup_manifest <- tibble::tibble(
  supplier_lookup_id =
    paste0(
      "SUPLOOKUP_",
      compound_detail_manifest$
        compound_id
    ),

  compound_id =
    compound_detail_manifest$
      compound_id,

  compound_name =
    compound_detail_manifest$
      compound_name,

  compound_origin =
    compound_detail_manifest$
      compound_origin,

  chembl_id =
    compound_detail_manifest$
      chembl_id,

  pubchem_cid =
    compound_detail_manifest$
      pubchem_cid,

  full_inchi_key =
    compound_detail_manifest$
      full_inchi_key,

  parent_connectivity_key =
    compound_detail_manifest$
      parent_connectivity_key,

  canonical_smiles =
    compound_detail_manifest$
      canonical_smiles,

  analysis_smiles =
    compound_detail_manifest$
      analysis_smiles,

  cas_number =
    NA_character_,

  maximum_clinical_phase =
    compound_detail_manifest$
      max_phase,

  scientific_priority_status =
    compound_detail_manifest$
      progression_status,

  lookup_identifier_type =
    lookup_identifier_type,

  pubchem_vendor_lookup_status =
    "not_started",

  molport_lookup_status =
    "not_started",

  emolecules_lookup_status =
    "not_started",

  institutional_supplier_lookup_status =
    "not_started",

  commercial_listing_present =
    as.logical(NA),

  exact_product_present =
    as.logical(NA),

  supplier_count =
    NA_real_,

  commercial_product_count =
    NA_real_,

  identity_match_class =
    "not_yet_assessed",

  manual_identity_review_required =
    TRUE,

  last_checked =
    NA_character_,

  lookup_notes =
    NA_character_
)

if (
  nrow(
    supplier_lookup_manifest
  ) != 16L
) {
  stop(
    "Supplier lookup manifest must contain 16 rows.",
    call. = FALSE
  )
}

cat(
  "Supplier lookup manifest completed: 16 rows.\n"
)

# ============================================================
# 15. Assess data completeness
# ============================================================

cat(
  "Assessing compound-data completeness...\n"
)

has_chembl_id <- !is.na(
  compound_detail_manifest$chembl_id
)

has_pubchem_cid <- !is.na(
  compound_detail_manifest$pubchem_cid
)

has_full_inchi_key <- !is.na(
  compound_detail_manifest$
    full_inchi_key
)

has_parent_connectivity_key <- !is.na(
  compound_detail_manifest$
    parent_connectivity_key
)

has_canonical_smiles <- !is.na(
  compound_detail_manifest$
    canonical_smiles
)

has_analysis_smiles <- !is.na(
  compound_detail_manifest$
    analysis_smiles
)

has_biological_classification <- !is.na(
  compound_detail_manifest$
    biological_classification
)

has_similarity_information <- !is.na(
  compound_detail_manifest$
    nearest_query_tanimoto
)

has_sider_information <-
  compound_detail_manifest$
    sider_representation ==
    "represented_in_sider_4.1"

has_clinical_phase <- !is.na(
  compound_detail_manifest$max_phase
)

supplier_lookup_ready <-
  has_full_inchi_key |
  has_pubchem_cid |
  has_canonical_smiles

missing_critical_identifier <-
  !has_full_inchi_key &
  !has_pubchem_cid &
  !has_canonical_smiles

completeness_status <- ifelse(
  missing_critical_identifier,
  "critical_identifier_review_required",
  ifelse(
    !has_full_inchi_key,
    "inchi_key_review_required",
    ifelse(
      !has_chembl_id,
      "chembl_mapping_optional_or_required",
      "core_identity_complete"
    )
  )
)

completeness_review <- tibble::tibble(
  compound_id =
    compound_detail_manifest$
      compound_id,

  compound_name =
    compound_detail_manifest$
      compound_name,

  compound_origin =
    compound_detail_manifest$
      compound_origin,

  has_chembl_id =
    has_chembl_id,

  has_pubchem_cid =
    has_pubchem_cid,

  has_full_inchi_key =
    has_full_inchi_key,

  has_parent_connectivity_key =
    has_parent_connectivity_key,

  has_canonical_smiles =
    has_canonical_smiles,

  has_analysis_smiles =
    has_analysis_smiles,

  has_biological_classification =
    has_biological_classification,

  has_similarity_information =
    has_similarity_information,

  has_sider_information =
    has_sider_information,

  has_clinical_phase =
    has_clinical_phase,

  supplier_lookup_ready =
    supplier_lookup_ready,

  missing_critical_identifier =
    missing_critical_identifier,

  completeness_status =
    completeness_status
)

if (
  nrow(
    completeness_review
  ) != 16L
) {
  stop(
    "Completeness review must contain 16 rows.",
    call. = FALSE
  )
}

cat(
  "Data-completeness review completed: 16 rows.\n"
)

# ============================================================
# 16. Create summary table
# ============================================================

cat(
  "Creating the A18A summary table...\n"
)

summary_table <- tibble::tibble(
  metric = c(
    "compound_count",
    "a16_compound_count",
    "a17_comparator_count",
    "compounds_with_chembl_id",
    "compounds_with_pubchem_cid",
    "compounds_with_full_inchi_key",
    "compounds_with_parent_connectivity_key",
    "compounds_with_analysis_smiles",
    "compounds_represented_in_sider",
    "compounds_with_clinical_phase",
    "supplier_lookup_ready_count",
    "critical_identifier_review_count",
    "live_network_detail_ready_count",
    "supplier_checks_completed"
  ),

  value = c(
    as.character(
      nrow(
        compound_detail_manifest
      )
    ),

    as.character(
      sum(
        compound_detail_manifest$
          compound_origin ==
          "A16_current_or_reference_compound"
      )
    ),

    as.character(
      sum(
        compound_detail_manifest$
          compound_origin ==
          "A17_SIDER_chemical_neighbour"
      )
    ),

    as.character(
      sum(
        has_chembl_id
      )
    ),

    as.character(
      sum(
        has_pubchem_cid
      )
    ),

    as.character(
      sum(
        has_full_inchi_key
      )
    ),

    as.character(
      sum(
        has_parent_connectivity_key
      )
    ),

    as.character(
      sum(
        has_analysis_smiles
      )
    ),

    as.character(
      sum(
        has_sider_information
      )
    ),

    as.character(
      sum(
        has_clinical_phase
      )
    ),

    as.character(
      sum(
        supplier_lookup_ready
      )
    ),

    as.character(
      sum(
        missing_critical_identifier
      )
    ),

    as.character(
      sum(
        compound_detail_manifest$
          live_network_detail_ready
      )
    ),

    "0"
  )
)

if (nrow(summary_table) != 14L) {
  stop(
    "A18A summary must contain 14 rows.",
    call. = FALSE
  )
}

cat(
  "A18A summary table completed: 14 rows.\n"
)

# ============================================================
# 17. Write outputs
# ============================================================

cat(
  "Writing A18A outputs...\n"
)

safe_write_csv(
  compound_detail_manifest,
  outputs[["detail_manifest"]]
)

cat(
  "  Wrote compound-detail manifest.\n"
)

safe_write_csv(
  compound_identifiers,
  outputs[["identifiers"]]
)

cat(
  "  Wrote compound-identifier manifest.\n"
)

safe_write_csv(
  compound_evidence_summary,
  outputs[["evidence"]]
)

cat(
  "  Wrote compound-evidence summary.\n"
)

safe_write_csv(
  supplier_lookup_manifest,
  outputs[["supplier_lookup"]]
)

cat(
  "  Wrote supplier-lookup manifest.\n"
)

safe_write_csv(
  completeness_review,
  outputs[["completeness"]]
)

cat(
  "  Wrote data-completeness review.\n"
)

safe_write_csv(
  summary_table,
  outputs[["summary"]]
)

cat(
  "  Wrote A18A summary.\n"
)

# ============================================================
# 18. Validate written outputs
# ============================================================

cat(
  "Validating written A18A outputs...\n"
)

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

detail_check <- readr::read_csv(
  outputs[["detail_manifest"]],
  show_col_types = FALSE,
  progress = FALSE
)

identifier_check <- readr::read_csv(
  outputs[["identifiers"]],
  show_col_types = FALSE,
  progress = FALSE
)

evidence_check <- readr::read_csv(
  outputs[["evidence"]],
  show_col_types = FALSE,
  progress = FALSE
)

supplier_check <- readr::read_csv(
  outputs[["supplier_lookup"]],
  show_col_types = FALSE,
  progress = FALSE
)

completeness_check <- readr::read_csv(
  outputs[["completeness"]],
  show_col_types = FALSE,
  progress = FALSE
)

summary_check <- readr::read_csv(
  outputs[["summary"]],
  show_col_types = FALSE,
  progress = FALSE
)

expected_row_counts <- c(
  detail_manifest = 16L,
  identifiers = 16L,
  evidence = 16L,
  supplier_lookup = 16L,
  completeness = 16L,
  summary = 14L
)

observed_row_counts <- c(
  detail_manifest =
    nrow(detail_check),

  identifiers =
    nrow(identifier_check),

  evidence =
    nrow(evidence_check),

  supplier_lookup =
    nrow(supplier_check),

  completeness =
    nrow(completeness_check),

  summary =
    nrow(summary_check)
)

incorrect_tables <- names(
  expected_row_counts
)[
  expected_row_counts !=
    observed_row_counts
]

if (length(incorrect_tables) > 0L) {
  stop(
    "Written row-count mismatch for: ",
    paste(
      incorrect_tables,
      collapse = ", "
    ),
    call. = FALSE
  )
}

if (
  anyDuplicated(
    detail_check$compound_id
  ) > 0L
) {
  stop(
    "Written detail manifest contains duplicate compound IDs.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    identifier_check$compound_id
  ) > 0L
) {
  stop(
    "Written identifier manifest contains duplicate compound IDs.",
    call. = FALSE
  )
}

if (
  anyDuplicated(
    supplier_check$supplier_lookup_id
  ) > 0L
) {
  stop(
    "Written supplier manifest contains duplicate lookup IDs.",
    call. = FALSE
  )
}

cat(
  "All A18A output validations passed.\n"
)

# ============================================================
# 19. Completion summary
# ============================================================

cat(
  "\nA18A compound-detail manifest completed.\n"
)

cat(
  "Compounds integrated: ",
  nrow(compound_detail_manifest),
  "\n",
  sep = ""
)

cat(
  "A16 compounds: ",
  sum(
    compound_detail_manifest$
      compound_origin ==
      "A16_current_or_reference_compound"
  ),
  "\n",
  sep = ""
)

cat(
  "A17 comparators: ",
  sum(
    compound_detail_manifest$
      compound_origin ==
      "A17_SIDER_chemical_neighbour"
  ),
  "\n",
  sep = ""
)

cat(
  "Compounds with ChEMBL IDs: ",
  sum(
    has_chembl_id
  ),
  "\n",
  sep = ""
)

cat(
  "Compounds with PubChem CIDs: ",
  sum(
    has_pubchem_cid
  ),
  "\n",
  sep = ""
)

cat(
  "Compounds with full InChIKeys: ",
  sum(
    has_full_inchi_key
  ),
  "\n",
  sep = ""
)

cat(
  "Compounds ready for supplier lookup: ",
  sum(
    supplier_lookup_ready
  ),
  "\n",
  sep = ""
)

cat(
  "Critical identifier reviews: ",
  sum(
    missing_critical_identifier
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