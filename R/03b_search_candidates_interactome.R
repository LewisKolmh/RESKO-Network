#' RESKO Faithful: Step 3b - Build interactome-based candidate pool (SE-independent)
#'
#' Input: data/raw/interactome_drug_target_links.csv (drug-target edges for
#'        every node in the eEF1A interactome)
#'        data/raw/interactome_nodes.csv (first_shell / is_seed flags per node)
#' Output: data/processed/candidate_drugs_no_se.csv
#'
#' WHY THIS SCRIPT EXISTS (see 03_search_candidate_drugs.R): McGarry's
#' literal SE-coverage filter is degenerate for this seed set -- only 2
#' seeds survive SE pruning and they share 0 common side-effects, so the
#' >25%-coverage filter yields 0 candidates by construction. This script
#' builds the candidate pool actually used for Variant A/B scoring
#' (04_score_candidates_jaccard.R / 04b_score_variantB_interactome.R)
#' instead: every drug with a bioactivity edge into the eEF1A interactome
#' (first-shell partners or the broader network), independent of any
#' side-effect data. This is the disclosed methodological substitution
#' described in METHODS_DEVIATIONS.md.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

source(file.path("R", "seed_compounds.R"))

DATA_DIR <- file.path("data", "raw")
PROCESSED_DIR <- file.path("data", "processed")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: BUILD INTERACTOME-BASED CANDIDATE POOL (SE-INDEPENDENT)\n")
cat(strrep("=", 70), "\n")

links <- read_csv(file.path(DATA_DIR, "interactome_drug_target_links.csv"), show_col_types = FALSE)
nodes <- read_csv(file.path(DATA_DIR, "interactome_nodes.csv"), show_col_types = FALSE)

seed_ids <- names(SEED_COMPOUNDS)
seed_names_upper <- str_replace_all(toupper(seed_ids), "_", " ")

mask_seed <- links$drug_id %in% seed_ids | toupper(links$drug_name) %in% seed_names_upper
cand_links <- links[!mask_seed, ]
cat(sprintf("\n  %d candidate drug-target edges remain after excluding %d seed compound(s) (%d edges excluded)\n",
            nrow(cand_links), length(seed_ids), sum(mask_seed)))

# Aggregate per drug: which interactome nodes (targets) it hits, is any of
# them first-shell/seed
agg <- cand_links %>%
  group_by(drug_id, drug_name, drug_type, max_stage) %>%
  summarise(
    n_interactome_targets_hit = n_distinct(target),
    targets_hit = paste(sort(unique(target)), collapse = ";"),
    moa = dplyr::first(moa),
    moa_n_targets = dplyr::first(moa_n_targets),
    .groups = "drop"
  )

# Flag whether any hit target is first-shell / a seed target (all targets in
# this link table ARE interactome nodes by construction)
node_first_shell <- nodes$symbol[nodes$first_shell]
node_seed <- nodes$symbol[nodes$is_seed]

any_in <- function(targets_str, pool) {
  any(str_split(targets_str, ";")[[1]] %in% pool)
}

agg <- agg %>%
  rowwise() %>%
  mutate(
    hits_first_shell = any_in(targets_hit, node_first_shell),
    hits_eef1a_seed = any_in(targets_hit, node_seed)
  ) %>%
  ungroup() %>%
  arrange(desc(n_interactome_targets_hit))

write_csv(agg, file.path(PROCESSED_DIR, "candidate_drugs_no_se.csv"))

cat(sprintf("\n  Saved %d candidates to data/processed/candidate_drugs_no_se.csv\n", nrow(agg)))
cat(sprintf("  Candidates hitting an eEF1A seed protein directly: %d\n", sum(agg$hits_eef1a_seed)))
cat(sprintf("  Candidates hitting only first-shell partners: %d\n", sum(agg$hits_first_shell & !agg$hits_eef1a_seed)))
cat("\n  Top 10 by number of interactome targets hit:\n")
top10 <- head(agg, 10) %>% select(drug_id, drug_name, max_stage, n_interactome_targets_hit, hits_eef1a_seed)
for (i in seq_len(nrow(top10))) {
  r <- top10[i, ]
  cat(sprintf("    %-16s %-25s stage=%-10s targets=%-4d seed_hit=%s\n",
              r$drug_id, r$drug_name, r$max_stage, r$n_interactome_targets_hit, r$hits_eef1a_seed))
}

cat("\nInteractome-based candidate pool build complete\n")
cat("  Next: Run 04_score_candidates_jaccard.R and 04b_score_variantB_interactome.R\n")
