#' RESKO Faithful: Step 3 - Search Candidate Drugs (McGarry's SE-coverage filter)
#'
#' Input: Common seed side-effects (data/processed/seed_sideeffects.json)
#'        + SIDER4 all-drugs data
#' Output: Candidate drugs with >25% side-effect coverage against the
#'         common seed side-effect set (data/processed/candidate_drugs.csv)
#'
#' McGarry filter: only keep drugs sharing >25% of the common seed
#' side-effects.
#'
#' STRUCTURAL FINDING (reproduced faithfully from the Python original,
#' src/03_search_candidate_drugs.py): with only 2 SE-eligible seeds
#' surviving pruning (Molibresib, Plitidepsin) and ZERO side-effects in
#' common between them (see 02_extract_seed_sideeffects.R), McGarry's
#' >25%-of-common-SE filter is degenerate here: every candidate's overlap
#' with an empty common-SE set is necessarily empty, so coverage_percent is
#' 0% for every drug in SIDER4 and the candidate list is EMPTY by
#' construction, not due to a bug in the scan itself. The Python original
#' actually crashes at this point (KeyError: 'coverage_percent', trying to
#' sort_values on a column that doesn't exist in an empty DataFrame) --
#' this port reproduces the same empty result but handles it gracefully and
#' writes a header-only CSV instead of crashing, and prints an explicit
#' disclosure rather than silently succeeding.
#'
#' This degeneracy is the reason the pipeline's candidate generation moved
#' to the interactome-based approach for Variant A/B
#' (04_score_candidates_jaccard.R / 04b_score_variantB_interactome.R,
#' fed by candidate_drugs_no_se.csv) instead of this SE-coverage filter.
#' This script and its (empty) output are retained for methods completeness
#' and to show, faithfully, why McGarry's literal SE-coverage step could not
#' be used unmodified for this seed set.

suppressPackageStartupMessages({
  library(arrow)
  library(dplyr)
  library(readr)
  library(jsonlite)
})

source(file.path("R", "seed_compounds.R"))

DATA_DIR <- file.path("data", "raw")
PROCESSED_DIR <- file.path("data", "processed")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: SEARCH CANDIDATE DRUGS BY SIDE-EFFECT COVERAGE\n")
cat(strrep("=", 70), "\n")

# Load common seed side-effects
cat("\n[1/3] Loading seed side-effects\n")
seed_data <- fromJSON(file.path(PROCESSED_DIR, "seed_sideeffects.json"))
common_ses <- unique(unlist(seed_data$common_sideeffects))
final_seed_ids <- unique(unlist(seed_data$seed_drugs))

# Exclude candidates by DrugBank ID, not by seed_compounds key -- SIDER4/ADR
# tables are keyed on drugbank_id, and several seeds share no drugbank_id at
# all (NA), so only real, non-NA ids need to be excluded from the scan.
seed_drugbank_ids <- unique(na.omit(sapply(final_seed_ids, function(cid) {
  if (cid %in% names(SEED_COMPOUNDS)) SEED_COMPOUNDS[[cid]]$drugbank_id else NA
})))

cat(sprintf("  Common side-effects: %d\n", length(common_ses)))
cat(sprintf("  Seed drugs (post-pruning): %d -> %s\n", length(final_seed_ids),
            paste(sort(final_seed_ids), collapse = ", ")))
cat(sprintf("  Seed DrugBank IDs excluded from candidate scan: %s\n",
            paste(sort(seed_drugbank_ids), collapse = ", ")))

if (length(common_ses) == 0) {
  cat("\n  *** STRUCTURAL FINDING ***\n")
  cat("  Common side-effect set is EMPTY (0 common SE between the 2 surviving\n")
  cat("  seed drugs). McGarry's >25%-of-common-SE filter is therefore degenerate:\n")
  cat("  every candidate's overlap with an empty set is empty, so coverage_percent\n")
  cat("  = 0% for every drug in SIDER4 -- the candidate list below is empty BY\n")
  cat("  CONSTRUCTION, not because the scan failed. This is why candidate\n")
  cat("  generation for Variant A/B uses the interactome-based approach instead\n")
  cat("  (see candidate_drugs_no_se.csv / 04_score_candidates_jaccard.R).\n")
}

# Load SIDER data
cat("\n[2/3] Scanning SIDER4 for candidate drugs\n")
sider_df <- read_parquet(file.path(DATA_DIR, "sider_all_se.parquet"))

# Group by drug, count side-effects, apply McGarry's >25% coverage filter
candidates <- sider_df %>%
  filter(!(drugbank_id %in% seed_drugbank_ids)) %>%
  group_by(drugbank_id) %>%
  summarise(
    drug_name = dplyr::first(drug_name),
    total_ses = n_distinct(side_effect_name),
    n_matching_ses = n_distinct(side_effect_name[side_effect_name %in% common_ses]),
    .groups = "drop"
  ) %>%
  mutate(coverage_percent = if (length(common_ses) > 0) 100 * n_matching_ses / length(common_ses) else 0) %>%
  filter(coverage_percent > 25) %>%
  arrange(desc(coverage_percent))

cat(sprintf("\n  Candidates found with >25%% coverage: %d\n", nrow(candidates)))
if (nrow(candidates) > 0) {
  cat("\n  Top 10 by coverage:\n")
  top10 <- head(candidates, 10)
  for (i in seq_len(nrow(top10))) {
    r <- top10[i, ]
    cat(sprintf("    %-30s (%-10s): %5.1f%%\n", r$drug_name, r$drugbank_id, r$coverage_percent))
  }
} else {
  cat("  (empty -- see structural finding above)\n")
}

# Save -- write a header-only CSV rather than crashing when empty, so
# downstream steps / documentation can detect "ran, found nothing" instead
# of "never ran".
cat("\n[3/3] Saving candidate list\n")
if (nrow(candidates) == 0) {
  candidates <- tibble(drugbank_id = character(0), drug_name = character(0),
                        n_matching_ses = integer(0), coverage_percent = double(0),
                        total_ses = integer(0))
} else {
  candidates <- candidates %>% select(drugbank_id, drug_name, n_matching_ses, coverage_percent, total_ses)
}
write_csv(candidates, file.path(PROCESSED_DIR, "candidate_drugs.csv"))
cat(sprintf("  Saved %d candidates to data/processed/candidate_drugs.csv\n", nrow(candidates)))

cat("\nCandidate search complete\n")
if (nrow(candidates) == 0) {
  cat("  NOTE: 0 candidates from the SE-coverage filter (expected, see structural\n")
  cat("  finding above). Candidate generation for scoring instead uses the\n")
  cat("  interactome-based pool -- see candidate_drugs_no_se.csv.\n")
}
cat("  Next: Run 04_score_candidates_jaccard.R\n")
