# ============================================================
# A15 RANK eEF1A AND TRANSLATION-NETWORK CANDIDATES
#
# Inputs:
#   results/edges_compound_protein_corrected.csv
#   results/A11_candidate_molecule_annotations.csv
#
# Outputs:
#   results/A15_candidate_score_components.csv
#   results/A15_direct_eef1a1_ranking.csv
#   results/A15_translation_network_ranking.csv
#   results/A15_ranking_summary.csv
#   results/A15_candidate_rankings.html
# ============================================================

required_packages <- c("readr", "dplyr", "tibble", "htmltools")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(paste0(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))"
  ))
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
results_dir <- file.path(project_root, "results")
if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

edge_file <- file.path(results_dir, "edges_compound_protein_corrected.csv")
annotation_file <- file.path(results_dir, "A11_candidate_molecule_annotations.csv")

components_file <- file.path(results_dir, "A15_candidate_score_components.csv")
direct_file <- file.path(results_dir, "A15_direct_eef1a1_ranking.csv")
network_file <- file.path(results_dir, "A15_translation_network_ranking.csv")
summary_file <- file.path(results_dir, "A15_ranking_summary.csv")
report_file <- file.path(results_dir, "A15_candidate_rankings.html")

missing_files <- c(edge_file, annotation_file)[!file.exists(c(edge_file, annotation_file))]
if (length(missing_files) > 0) {
  stop(paste0("Missing inputs:\n", paste(missing_files, collapse = "\n")))
}

edges <- readr::read_csv(edge_file, show_col_types = FALSE, progress = FALSE)
annotations <- readr::read_csv(annotation_file, show_col_types = FALSE, progress = FALSE)

required_edge_columns <- c(
  "molecule_chembl_id", "display_name", "queried_protein", "relationship",
  "minimum_activity_nM", "unique_assay_count", "unique_document_count",
  "confidence_score", "development_status", "duplicate_kd_ed50_collapsed"
)
missing_columns <- setdiff(required_edge_columns, names(edges))
if (length(missing_columns) > 0) {
  stop(paste0("Missing corrected-edge columns: ", paste(missing_columns, collapse = ", ")))
}

optional_edge_columns <- c(
  "corrected_activity_types", "original_activity_types", "raw_record_count",
  "maximum_pchembl_value", "molecule_type", "max_phase"
)
for (column_name in optional_edge_columns) {
  if (!column_name %in% names(edges)) edges[[column_name]] <- NA
}

edges <- edges |>
  dplyr::mutate(
    display_name = dplyr::coalesce(display_name, molecule_chembl_id),
    minimum_activity_nM = suppressWarnings(as.numeric(minimum_activity_nM)),
    unique_assay_count = dplyr::coalesce(suppressWarnings(as.integer(unique_assay_count)), 0L),
    unique_document_count = dplyr::coalesce(suppressWarnings(as.integer(unique_document_count)), 0L),
    confidence_score = suppressWarnings(as.numeric(confidence_score)),
    max_phase = suppressWarnings(as.numeric(max_phase)),
    raw_record_count = dplyr::coalesce(suppressWarnings(as.integer(raw_record_count)), 0L)
  )

# ----------------------------------------------------------------
# Transparent scoring functions
# ----------------------------------------------------------------

directness_score <- function(protein) {
  dplyr::case_when(
    protein == "EEF1A1" ~ 4,
    protein == "EEF1A2" ~ 3,
    protein %in% c("EEF1G", "EEF1B2", "EEF1D") ~ 2,
    protein %in% c("RPS2", "RPS3", "RPS3A", "RPL3", "RPL4", "RPL7", "RPL18A") ~ 1,
    TRUE ~ 0
  )
}

evidence_score <- function(relationship) {
  dplyr::case_when(
    relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ 3,
    relationship == "BINDS_TO" ~ 2,
    relationship == "HAS_ACTIVITY_AGAINST" ~ 1,
    TRUE ~ 0
  )
}

potency_score <- function(value_nM) {
  dplyr::case_when(
    is.na(value_nM) ~ 0,
    value_nM <= 10 ~ 4,
    value_nM <= 100 ~ 3,
    value_nM <= 1000 ~ 2,
    value_nM <= 10000 ~ 1,
    TRUE ~ 0
  )
}

development_score <- function(max_phase, status) {
  dplyr::case_when(
    !is.na(max_phase) & max_phase >= 2 ~ 2,
    !is.na(max_phase) & max_phase == 1 ~ 1,
    grepl("Phase 2|Phase 3|Approved|launched", status, ignore.case = TRUE) ~ 2,
    grepl("Phase 1", status, ignore.case = TRUE) ~ 1,
    TRUE ~ 0
  )
}

edges_scored <- edges |>
  dplyr::mutate(
    target_directness_score = directness_score(queried_protein),
    evidence_strength_score = evidence_score(relationship),
    potency_score = potency_score(minimum_activity_nM),
    confidence_point = dplyr::if_else(!is.na(confidence_score) & confidence_score == 9, 1, 0),
    assay_point = dplyr::if_else(unique_assay_count >= 1, 1, 0),
    document_point = dplyr::if_else(unique_document_count >= 1, 1, 0),
    provenance_score = confidence_point + assay_point + document_point,
    clinical_development_score = development_score(max_phase, development_status),
    eef1_complex_weight = dplyr::case_when(
      queried_protein %in% c("EEF1A1", "EEF1A2") ~ 4,
      queried_protein %in% c("EEF1G", "EEF1B2", "EEF1D") ~ 3,
      queried_protein %in% c("RPS2", "RPS3", "RPS3A", "RPL3", "RPL4", "RPL7", "RPL18A") ~ 1,
      TRUE ~ 0
    )
  )

# Direct ranking: best evidence edge per compound, with direct-target preference.
direct_components <- edges_scored |>
  dplyr::group_by(molecule_chembl_id, display_name) |>
  dplyr::arrange(
    dplyr::desc(target_directness_score),
    dplyr::desc(evidence_strength_score),
    dplyr::desc(potency_score),
    minimum_activity_nM,
    .by_group = TRUE
  ) |>
  dplyr::summarise(
    direct_primary_target = dplyr::first(queried_protein),
    direct_primary_relationship = dplyr::first(relationship),
    direct_minimum_activity_nM = dplyr::first(minimum_activity_nM),
    direct_target_score = dplyr::first(target_directness_score),
    direct_evidence_score = dplyr::first(evidence_strength_score),
    direct_potency_score = dplyr::first(potency_score),
    direct_provenance_score = dplyr::first(provenance_score),
    direct_development_score = max(clinical_development_score, na.rm = TRUE),
    direct_assay_count = dplyr::first(unique_assay_count),
    direct_document_count = dplyr::first(unique_document_count),
    direct_confidence_score = dplyr::first(confidence_score),
    has_direct_eef1a1_evidence = any(queried_protein == "EEF1A1"),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    direct_total_score =
      direct_target_score + direct_evidence_score + direct_potency_score +
      direct_provenance_score + direct_development_score,
    direct_eligibility = dplyr::if_else(
      has_direct_eef1a1_evidence,
      "Direct EEF1A1 evidence",
      "Network-only evidence"
    )
  )

# Translation-network ranking: relevance-weighted coverage, potency, provenance,
# development, plus a modest breadth penalty to avoid rewarding promiscuity alone.
network_components <- edges_scored |>
  dplyr::group_by(molecule_chembl_id, display_name) |>
  dplyr::summarise(
    network_target_count = dplyr::n_distinct(queried_protein),
    eef1_target_count = dplyr::n_distinct(
      queried_protein[queried_protein %in% c("EEF1A1", "EEF1A2", "EEF1G", "EEF1B2", "EEF1D")]
    ),
    ribosomal_target_count = dplyr::n_distinct(
      queried_protein[queried_protein %in% c("RPS2", "RPS3", "RPS3A", "RPL3", "RPL4", "RPL7", "RPL18A")]
    ),
    strongest_network_target = queried_protein[which.max(eef1_complex_weight * 100 + potency_score)][1],
    strongest_network_relationship = relationship[which.max(eef1_complex_weight * 100 + potency_score)][1],
    strongest_network_activity_nM = {
      vals <- minimum_activity_nM[is.finite(minimum_activity_nM)]
      if (length(vals) == 0) NA_real_ else min(vals)
    },
    network_relevance_score = min(4, max(eef1_complex_weight, na.rm = TRUE)),
    network_coverage_score = min(4, eef1_target_count * 2 + min(2, ribosomal_target_count / 2)),
    network_potency_score = max(potency_score, na.rm = TRUE),
    network_provenance_score = min(3, max(provenance_score, na.rm = TRUE)),
    network_development_score = max(clinical_development_score, na.rm = TRUE),
    network_assay_count = sum(unique_assay_count, na.rm = TRUE),
    network_document_count = sum(unique_document_count, na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::mutate(
    # Penalty begins only once a compound engages more than four network targets.
    network_breadth_penalty = pmin(3, pmax(0, network_target_count - 4) * 0.5),
    network_total_score =
      network_relevance_score + network_coverage_score + network_potency_score +
      network_provenance_score + network_development_score - network_breadth_penalty
  )

score_components <- direct_components |>
  dplyr::full_join(
    network_components,
    by = c("molecule_chembl_id", "display_name")
  ) |>
  dplyr::mutate(
    overall_priority_score = direct_total_score + network_total_score,
    priority_category = dplyr::case_when(
      has_direct_eef1a1_evidence & network_total_score >= 10 ~ "Dual priority",
      has_direct_eef1a1_evidence ~ "Direct priority",
      network_total_score >= 10 ~ "Network priority",
      TRUE ~ "Supporting candidate"
    )
  )

direct_ranking <- score_components |>
  dplyr::arrange(
    dplyr::desc(has_direct_eef1a1_evidence),
    dplyr::desc(direct_total_score),
    direct_minimum_activity_nM,
    display_name
  ) |>
  dplyr::mutate(direct_rank = dplyr::row_number()) |>
  dplyr::select(
    direct_rank, molecule_chembl_id, display_name, direct_eligibility,
    direct_primary_target, direct_primary_relationship,
    direct_minimum_activity_nM, direct_target_score, direct_evidence_score,
    direct_potency_score, direct_provenance_score, direct_development_score,
    direct_total_score, priority_category
  )

network_ranking <- score_components |>
  dplyr::arrange(
    dplyr::desc(network_total_score),
    strongest_network_activity_nM,
    display_name
  ) |>
  dplyr::mutate(network_rank = dplyr::row_number()) |>
  dplyr::select(
    network_rank, molecule_chembl_id, display_name,
    strongest_network_target, strongest_network_relationship,
    strongest_network_activity_nM, network_target_count, eef1_target_count,
    ribosomal_target_count, network_relevance_score, network_coverage_score,
    network_potency_score, network_provenance_score, network_development_score,
    network_breadth_penalty, network_total_score, priority_category
  )

ranking_summary <- tibble::tibble(
  metric = c(
    "Candidates ranked",
    "Candidates with direct EEF1A1 evidence",
    "Direct-priority candidates",
    "Network-priority candidates",
    "Dual-priority candidates"
  ),
  count = c(
    nrow(score_components),
    sum(score_components$has_direct_eef1a1_evidence, na.rm = TRUE),
    sum(score_components$priority_category == "Direct priority", na.rm = TRUE),
    sum(score_components$priority_category == "Network priority", na.rm = TRUE),
    sum(score_components$priority_category == "Dual priority", na.rm = TRUE)
  )
)

readr::write_csv(score_components, components_file, na = "")
readr::write_csv(direct_ranking, direct_file, na = "")
readr::write_csv(network_ranking, network_file, na = "")
readr::write_csv(ranking_summary, summary_file, na = "")

# HTML report.
css <- htmltools::tags$style(htmltools::HTML("\nbody{font-family:Calibri,Arial,sans-serif;margin:28px;color:#222;background:#fff;}\nh1,h2{font-weight:600;}\ntable{border-collapse:collapse;width:100%;margin:12px 0 28px 0;font-size:13px;}\nth{background:#eef2f6;text-align:left;padding:8px;border:1px solid #ccd3da;}\ntd{padding:8px;border:1px solid #d9dee3;}\n.note{background:#f7f8fa;border-left:4px solid #5b8fd1;padding:12px;margin:12px 0;}\n"))

table_tag <- function(data) {
  htmltools::tags$table(
    htmltools::tags$thead(
      htmltools::tags$tr(lapply(names(data), htmltools::tags$th))
    ),
    htmltools::tags$tbody(
      lapply(seq_len(nrow(data)), function(i) {
        htmltools::tags$tr(lapply(data[i, , drop = FALSE], function(x) htmltools::tags$td(as.character(x))))
      })
    )
  )
}

report <- htmltools::tags$html(
  htmltools::tags$head(
    htmltools::tags$title("A15 Candidate Rankings"),
    css
  ),
  htmltools::tags$body(
    htmltools::tags$h1("A15 Candidate Rankings"),
    htmltools::tags$p("Two complementary rankings are reported: direct EEF1A1 candidates and broader translation-network candidates."),
    htmltools::tags$div(class = "note", "Binding and inhibitory-activity evidence are ranked separately. Provenance-corrected Kd/ED50 duplicates contribute once."),
    htmltools::tags$h2("Direct EEF1A1 ranking"),
    table_tag(direct_ranking),
    htmltools::tags$h2("Translation-network ranking"),
    table_tag(network_ranking),
    htmltools::tags$h2("Scoring interpretation"),
    htmltools::tags$ul(
      htmltools::tags$li("Direct score: target directness + evidence strength + potency + provenance + development maturity."),
      htmltools::tags$li("Network score: eEF1 relevance + translation-network coverage + potency + provenance + development maturity - breadth penalty."),
      htmltools::tags$li("The direct ranking is primary for compounds intended to engage EEF1A1 itself."),
      htmltools::tags$li("The network ranking identifies broader perturbation of eEF1 and ribosomal machinery.")
    )
  )
)

htmltools::save_html(report, report_file)

output_files <- c(components_file, direct_file, network_file, summary_file, report_file)
file_check <- tibble::tibble(
  file = output_files,
  exists = file.exists(output_files),
  size_bytes = as.numeric(file.info(output_files)$size)
)
if (!all(file_check$exists)) stop("One or more A15 output files were not created.")
if (any(is.na(file_check$size_bytes) | file_check$size_bytes <= 0)) stop("One or more A15 files are empty.")

cat("\nA15 candidate ranking completed successfully.\n")
cat("--------------------------------------------\n")
cat("\nDirect EEF1A1 ranking:\n")
print(direct_ranking, n = Inf, width = Inf)
cat("\nTranslation-network ranking:\n")
print(network_ranking, n = Inf, width = Inf)
cat("\nOutput verification:\n")
print(file_check, n = Inf)
