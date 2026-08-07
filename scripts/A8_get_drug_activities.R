# ============================================================
# A9 RETRIEVE COMPLETE ChEMBL ACTIVITY DATA
#
# Retrieves all available activity pages for each mapped
# ChEMBL target and preserves fields required for filtering.
#
# Input:
# results/A8_chembl_target_map.csv
#
# Outputs:
# results/A9_chembl_activities_complete.csv
# results/A9_chembl_retrieval_summary.csv
# ============================================================


# ============================================================
# 1. REQUIRED PACKAGES
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
      "\nInstall them with:\n",
      "install.packages(c(",
      paste0(
        '"',
        missing_packages,
        '"',
        collapse = ", "
      ),
      "))"
    )
  )
}


# ============================================================
# 2. FILE PATHS
# ============================================================

target_file <-
  "results/A8_chembl_target_map.csv"

activity_output_file <-
  "results/A9_chembl_activities_complete.csv"

summary_output_file <-
  "results/A9_chembl_retrieval_summary.csv"


# ============================================================
# 3. LOAD TARGET MAP
# ============================================================

if (!file.exists(target_file)) {
  stop(
    paste0(
      "Target map not found: ",
      target_file
    )
  )
}

targets <- readr::read_csv(
  target_file,
  show_col_types = FALSE
)

required_columns <- c(
  "protein",
  "chembl_target_id"
)

missing_columns <- setdiff(
  required_columns,
  names(targets)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing columns in target map: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

targets <- targets |>
  dplyr::filter(
    !is.na(protein),
    !is.na(chembl_target_id),
    protein != "",
    chembl_target_id != ""
  ) |>
  dplyr::distinct(
    protein,
    chembl_target_id,
    .keep_all = TRUE
  )


# ============================================================
# 4. API SETTINGS
# ============================================================

base_url <-
  "https://www.ebi.ac.uk/chembl/api/data/activity.json"

page_limit <- 1000

pause_between_requests <- 0.25

maximum_retries <- 3


# ============================================================
# 5. SAFE API REQUEST FUNCTION
# ============================================================

request_chembl_page <- function(
  target_id,
  offset,
  limit
) {

  response <- NULL

  for (attempt in seq_len(maximum_retries)) {

    response <- tryCatch(
      httr::GET(
        url = base_url,
        query = list(
          target_chembl_id = target_id,
          limit = limit,
          offset = offset
        ),
        httr::timeout(120)
      ),
      error = function(error) {
        NULL
      }
    )

    if (
      !is.null(response) &&
      httr::status_code(response) == 200
    ) {
      break
    }

    message(
      "Request failed for ",
      target_id,
      " at offset ",
      offset,
      ". Attempt ",
      attempt,
      " of ",
      maximum_retries,
      "."
    )

    Sys.sleep(attempt * 2)
  }

  if (
    is.null(response) ||
    httr::status_code(response) != 200
  ) {
    return(NULL)
  }

  response_text <- httr::content(
    response,
    as = "text",
    encoding = "UTF-8"
  )

  jsonlite::fromJSON(
    response_text,
    flatten = TRUE
  )
}


# ============================================================
# 6. RETRIEVE ALL PAGES FOR EACH TARGET
# ============================================================

all_activities <- list()

retrieval_summary <- list()

activity_index <- 1

summary_index <- 1

for (target_row in seq_len(nrow(targets))) {

  protein_name <-
    targets$protein[target_row]

  target_id <-
    targets$chembl_target_id[target_row]

  message(
    "\nRetrieving activities for ",
    protein_name,
    " [",
    target_id,
    "]"
  )

  target_pages <- list()

  page_index <- 1

  offset <- 0

  target_total_count <- NA_integer_

  repeat {

    api_result <- request_chembl_page(
      target_id = target_id,
      offset = offset,
      limit = page_limit
    )

    if (is.null(api_result)) {

      warning(
        "Unable to retrieve activities for ",
        target_id,
        " at offset ",
        offset,
        "."
      )

      break
    }

    if (
      is.na(target_total_count) &&
      !is.null(api_result$page_meta$total_count)
    ) {
      target_total_count <-
        as.integer(
          api_result$page_meta$total_count
        )
    }

    if (
      is.null(api_result$activities) ||
      nrow(api_result$activities) == 0
    ) {
      break
    }

    page_data <- tibble::as_tibble(
      api_result$activities
    ) |>
      dplyr::mutate(
        queried_protein = protein_name,
        queried_target_chembl_id = target_id,
        .before = 1
      )

    target_pages[[page_index]] <-
      page_data

    records_received <-
      nrow(page_data)

    message(
      "  Offset ",
      offset,
      ": ",
      records_received,
      " records"
    )

    page_index <-
      page_index + 1

    offset <-
      offset + page_limit

    if (
      !is.na(target_total_count) &&
      offset >= target_total_count
    ) {
      break
    }

    Sys.sleep(
      pause_between_requests
    )
  }

  target_activities <- dplyr::bind_rows(
    target_pages
  )

  all_activities[[activity_index]] <-
    target_activities

  activity_index <-
    activity_index + 1

  retrieval_summary[[summary_index]] <-
    tibble::tibble(
      protein = protein_name,
      chembl_target_id = target_id,
      api_total_count = target_total_count,
      records_retrieved = nrow(
        target_activities
      )
    )

  summary_index <-
    summary_index + 1
}


# ============================================================
# 7. COMBINE RESULTS
# ============================================================

activities_complete <- dplyr::bind_rows(
  all_activities
)

retrieval_summary_df <- dplyr::bind_rows(
  retrieval_summary
)

if (nrow(activities_complete) == 0) {
  stop(
    paste0(
      "No ChEMBL activity records were retrieved. ",
      "Check the target identifiers and network connection."
    )
  )
}


# ============================================================
# 8. SAVE COMPLETE ACTIVITY DATA
# ============================================================

readr::write_csv(
  activities_complete,
  activity_output_file
)

readr::write_csv(
  retrieval_summary_df,
  summary_output_file
)


# ============================================================
# 9. PRINT SUMMARY
# ============================================================

cat("\n")
cat("ChEMBL retrieval completed.\n")
cat("---------------------------\n")

cat(
  "Mapped targets queried: ",
  nrow(targets),
  "\n"
)

cat(
  "Activity records:       ",
  nrow(activities_complete),
  "\n"
)

if (
  "molecule_chembl_id" %in%
  names(activities_complete)
) {
  cat(
    "Unique molecules:       ",
    dplyr::n_distinct(
      activities_complete$molecule_chembl_id
    ),
    "\n"
  )
}

cat(
  "Complete output:        ",
  activity_output_file,
  "\n"
)

cat(
  "Retrieval summary:      ",
  summary_output_file,
  "\n\n"
)

print(
  retrieval_summary_df
)