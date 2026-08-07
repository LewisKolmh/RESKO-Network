# ============================================================
# A11 ANNOTATE CANDIDATE ChEMBL MOLECULES
#
# Input:
#   results/A10_candidate_molecule_ids.csv
#
# Outputs:
#   results/A11_candidate_molecule_annotations.csv
#   results/nodes_drugs.csv
#   results/A11_molecule_retrieval_summary.csv
# ============================================================

# ============================================================
# 1. CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "httr",
  "jsonlite",
  "tibble"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them with:\n",
      "install.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

# ============================================================
# 2. FILE PATHS
# ============================================================

input_file <-
  "results/A10_candidate_molecule_ids.csv"

annotation_file <-
  "results/A11_candidate_molecule_annotations.csv"

nodes_file <-
  "results/nodes_drugs.csv"

retrieval_summary_file <-
  "results/A11_molecule_retrieval_summary.csv"

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

if (!file.exists(input_file)) {
  stop(
    paste0(
      "Input file not found: ",
      input_file,
      "\n\nCurrent working directory:\n",
      getwd()
    )
  )
}

# ============================================================
# 3. LOAD CANDIDATE MOLECULE IDS
# ============================================================

candidates <- readr::read_csv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)

required_columns <- c(
  "molecule_chembl_id",
  "proteins",
  "activity_types",
  "strongest_record_nM",
  "total_measurements"
)

missing_columns <- setdiff(
  required_columns,
  names(candidates)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing columns in candidate file: ",
      paste(missing_columns, collapse = ", ")
    )
  )
}

candidates <- candidates |>
  dplyr::filter(
    !is.na(molecule_chembl_id),
    molecule_chembl_id != ""
  ) |>
  dplyr::distinct(
    molecule_chembl_id,
    .keep_all = TRUE
  )

if (nrow(candidates) == 0) {
  stop("No candidate molecule IDs were found.")
}

# ============================================================
# 4. HELPER FUNCTIONS
# ============================================================

null_to_na <- function(value) {
  if (is.null(value) || length(value) == 0) {
    return(NA)
  }

  if (length(value) > 1) {
    value <- value[[1]]
  }

  if (is.null(value) || length(value) == 0) {
    return(NA)
  }

  value
}

character_or_na <- function(value) {
  value <- null_to_na(value)

  if (length(value) == 1 && is.na(value)) {
    return(NA_character_)
  }

  as.character(value)
}

numeric_or_na <- function(value) {
  value <- null_to_na(value)

  if (length(value) == 1 && is.na(value)) {
    return(NA_real_)
  }

  suppressWarnings(as.numeric(value))
}

logical_or_na <- function(value) {
  value <- null_to_na(value)

  if (length(value) == 1 && is.na(value)) {
    return(NA)
  }

  as.logical(value)
}

safe_nested <- function(object, keys) {
  current <- object

  for (key in keys) {
    if (is.null(current) || is.null(current[[key]])) {
      return(NULL)
    }

    current <- current[[key]]
  }

  current
}

request_molecule <- function(
  molecule_id,
  maximum_retries = 3
) {
  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/molecule/",
    molecule_id,
    ".json"
  )

  for (attempt in seq_len(maximum_retries)) {
    response <- tryCatch(
      httr::GET(
        url,
        httr::timeout(60),
        httr::user_agent(
          "RESKO-eEF1A academic research"
        )
      ),
      error = function(error) NULL
    )

    if (
      !is.null(response) &&
      httr::status_code(response) == 200
    ) {
      text <- httr::content(
        response,
        as = "text",
        encoding = "UTF-8"
      )

      return(
        jsonlite::fromJSON(
          text,
          simplifyVector = FALSE
        )
      )
    }

    Sys.sleep(attempt * 2)
  }

  NULL
}

# ============================================================
# 5. RETRIEVE AND PARSE MOLECULE RECORDS
# ============================================================

annotation_rows <- list()
retrieval_rows <- list()

for (index in seq_len(nrow(candidates))) {
  molecule_id <- candidates$molecule_chembl_id[index]

  message(
    "Retrieving molecule ",
    index,
    " of ",
    nrow(candidates),
    ": ",
    molecule_id
  )

  molecule <- request_molecule(molecule_id)

  retrieval_rows[[index]] <- tibble::tibble(
    molecule_chembl_id = molecule_id,
    retrieved = !is.null(molecule)
  )

  if (is.null(molecule)) {
    annotation_rows[[index]] <- tibble::tibble(
      molecule_chembl_id = molecule_id,
      pref_name = NA_character_,
      display_name = molecule_id,
      molecule_type = NA_character_,
      max_phase = NA_real_,
      first_approval = NA_real_,
      therapeutic_flag = NA,
      natural_product = NA,
      prodrug = NA,
      oral = NA,
      parenteral = NA,
      topical = NA,
      parent_chembl_id = NA_character_,
      canonical_smiles = NA_character_,
      standard_inchi = NA_character_,
      standard_inchi_key = NA_character_,
      full_molecular_weight = NA_real_,
      molecular_weight_freebase = NA_real_,
      alogp = NA_real_,
      polar_surface_area = NA_real_,
      hydrogen_bond_acceptors = NA_real_,
      hydrogen_bond_donors = NA_real_,
      rule_of_five_violations = NA_real_
    )

    next
  }

  pref_name <- character_or_na(molecule$pref_name)

  display_name <- if (
    !is.na(pref_name) &&
    nzchar(pref_name)
  ) {
    pref_name
  } else {
    molecule_id
  }

  annotation_rows[[index]] <- tibble::tibble(
    molecule_chembl_id = molecule_id,
    pref_name = pref_name,
    display_name = display_name,
    molecule_type = character_or_na(molecule$molecule_type),
    max_phase = numeric_or_na(molecule$max_phase),
    first_approval = numeric_or_na(molecule$first_approval),
    therapeutic_flag = logical_or_na(molecule$therapeutic_flag),
    natural_product = logical_or_na(molecule$natural_product),
    prodrug = logical_or_na(molecule$prodrug),
    oral = logical_or_na(molecule$oral),
    parenteral = logical_or_na(molecule$parenteral),
    topical = logical_or_na(molecule$topical),
    parent_chembl_id = character_or_na(
      safe_nested(
        molecule,
        c("molecule_hierarchy", "parent_chembl_id")
      )
    ),
    canonical_smiles = character_or_na(
      safe_nested(
        molecule,
        c("molecule_structures", "canonical_smiles")
      )
    ),
    standard_inchi = character_or_na(
      safe_nested(
        molecule,
        c("molecule_structures", "standard_inchi")
      )
    ),
    standard_inchi_key = character_or_na(
      safe_nested(
        molecule,
        c("molecule_structures", "standard_inchi_key")
      )
    ),
    full_molecular_weight = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "full_mwt")
      )
    ),
    molecular_weight_freebase = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "mw_freebase")
      )
    ),
    alogp = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "alogp")
      )
    ),
    polar_surface_area = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "psa")
      )
    ),
    hydrogen_bond_acceptors = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "hba")
      )
    ),
    hydrogen_bond_donors = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "hbd")
      )
    ),
    rule_of_five_violations = numeric_or_na(
      safe_nested(
        molecule,
        c("molecule_properties", "ro5_violations")
      )
    )
  )

  Sys.sleep(0.25)
}

molecule_annotations <- dplyr::bind_rows(
  annotation_rows
)

retrieval_summary <- dplyr::bind_rows(
  retrieval_rows
)

# ============================================================
# 6. MERGE ANNOTATIONS WITH A10 EVIDENCE
# ============================================================

candidate_annotations <- candidates |>
  dplyr::left_join(
    molecule_annotations,
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    display_name = dplyr::coalesce(
      display_name,
      molecule_chembl_id
    ),
    development_status = dplyr::case_when(
      is.na(max_phase) ~ "Unknown",
      max_phase >= 4 ~ "Approved / launched",
      max_phase == 3 ~ "Phase 3",
      max_phase == 2 ~ "Phase 2",
      max_phase == 1 ~ "Phase 1",
      max_phase == 0 ~ "Preclinical / discovery",
      TRUE ~ "Other"
    )
  ) |>
  dplyr::arrange(
    strongest_record_nM,
    molecule_chembl_id
  )

# ============================================================
# 7. CREATE DRUG / COMPOUND NODES FOR THE GRAPH
# ============================================================

nodes_drugs <- candidate_annotations |>
  dplyr::transmute(
    molecule_chembl_id,
    name = display_name,
    type = "Compound",
    molecule_type,
    max_phase,
    development_status,
    first_approval,
    therapeutic_flag,
    natural_product,
    parent_chembl_id,
    canonical_smiles,
    standard_inchi_key,
    full_molecular_weight,
    alogp,
    polar_surface_area,
    hydrogen_bond_acceptors,
    hydrogen_bond_donors,
    rule_of_five_violations,
    proteins,
    activity_types,
    strongest_record_nM,
    total_measurements
  )

# ============================================================
# 8. WRITE AND VERIFY OUTPUTS
# ============================================================

readr::write_csv(
  candidate_annotations,
  annotation_file,
  na = ""
)

readr::write_csv(
  nodes_drugs,
  nodes_file,
  na = ""
)

readr::write_csv(
  retrieval_summary,
  retrieval_summary_file,
  na = ""
)

output_files <- c(
  annotation_file,
  nodes_file,
  retrieval_summary_file
)

file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(
    file.info(output_files)$size
  )
)

if (!all(file_check$exists)) {
  stop(
    "At least one A11 output file was not created."
  )
}

if (
  any(
    is.na(file_check$size_bytes) |
    file_check$size_bytes <= 0
  )
) {
  stop(
    "At least one A11 output file is empty or invalid."
  )
}

# ============================================================
# 9. PRINT SUMMARY
# ============================================================

cat("\n")
cat("A11 candidate molecule annotation completed.\n")
cat("--------------------------------------------\n")
cat("Candidate molecule IDs: ", nrow(candidates), "\n")
cat("Molecules retrieved:     ", sum(retrieval_summary$retrieved), "\n")
cat("Compound nodes created:  ", nrow(nodes_drugs), "\n\n")

cat("Development status:\n")
print(
  candidate_annotations |>
    dplyr::count(
      development_status,
      sort = TRUE
    ),
  n = Inf
)

cat("\nAnnotated candidates:\n")
print(
  candidate_annotations |>
    dplyr::select(
      molecule_chembl_id,
      display_name,
      proteins,
      activity_types,
      strongest_record_nM,
      molecule_type,
      max_phase,
      development_status
    ),
  n = Inf,
  width = Inf
)

cat("\nOutput verification:\n")
print(
  file_check,
  n = Inf
)
