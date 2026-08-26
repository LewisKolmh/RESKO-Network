#' RESKO Faithful: Step 4 - Score Candidates with McGarry's Jaccard Composite
#'
#' McGarry's three-component scoring:
#' 1. inds_score - therapeutic breadth (indications normalized)
#' 2. se_score - side-effect similarity coverage (%)
#' 3. on_score - on-target promiscuity (targets normalized)
#'
#' Combined via Matrix Jaccard Similarity:
#'   Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))
#'
#' This matches McGarry's reposition_validate_score.R exactly.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

PROCESSED_DIR <- file.path("data", "processed")
RAW_DIR <- file.path("data", "raw")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: SCORE CANDIDATES WITH MCGARRY'S JACCARD COMPOSITE\n")
cat(strrep("=", 70), "\n")

cat(paste(
  "",
  "NOTE ON se_score: real pharmacovigilance data (SIDER4 1430 drugs + FAERS,",
  "synonym-pooled, plus a ClinicalTrials.gov structured-AE check) was",
  "exhausted and confirmed the eEF1A-binder seed set has ZERO common",
  "side-effects across its SE-eligible seeds (Molibresib, Plitidepsin have",
  "data; Didemnin_B and Metarrestin have none in any source, for dated/",
  "trial-status reasons, not a retrieval failure). McGarry's se_score",
  "component is therefore fixed at 0 (undefined, not silently dropped) for",
  "every candidate in this run -- see METHODS_DEVIATIONS.md.", "", sep = "\n"
))

min_max_norm <- function(s) {
  lo <- min(s, na.rm = TRUE); hi <- max(s, na.rm = TRUE)
  if (hi == lo) return(rep(0.0, length(s)))
  (s - lo) / (hi - lo)
}

#' McGarry's Jaccard similarity for an N-component vector.
#' Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))
jaccard_composite <- function(scores_mat) {
  s <- rowSums(scores_mat^2)
  t <- rowSums(scores_mat)
  denom <- 2 * t - s
  out <- numeric(length(s))
  nz <- denom != 0
  out[nz] <- s[nz] / denom[nz]
  out
}

# ============================================================
# VARIANT A: indication-breadth + on-target-promiscuity only
# ============================================================
cat("\n[1/4] Loading SE-free candidates (Variant A input)\n")
cand <- read_csv(file.path(PROCESSED_DIR, "candidate_drugs_no_se.csv"), show_col_types = FALSE)
cat(sprintf("  Loaded %d candidates\n", nrow(cand)))

cat("\n[2/4] Computing real inds_score (ChEMBL indication breadth)\n")
n_ind <- read_csv(file.path(PROCESSED_DIR, "chembl_n_indications.csv"), show_col_types = FALSE)
n_ind_map <- setNames(n_ind$n_indications, n_ind$molecule_chembl_id)
cand$n_indications <- ifelse(cand$drug_id %in% names(n_ind_map), n_ind_map[cand$drug_id], 0)
cand$n_indications[is.na(cand$n_indications)] <- 0
cand$inds_score <- min_max_norm(cand$n_indications)

cat("Computing real on_score (ChEMBL on-target promiscuity)\n")
# moa_n_targets = number of distinct proteins in the drug's ChEMBL-recorded
# mechanism-of-action set -- a direct promiscuity measure, real ChEMBL data,
# not a placeholder.
moa_n_targets_filled <- ifelse(is.na(cand$moa_n_targets), 1, cand$moa_n_targets)
cand$on_score <- min_max_norm(moa_n_targets_filled)

cand$se_score_undefined <- 0.0  # explicitly undefined for ALL variants -- see note above

scores <- as.matrix(cand[, c("inds_score", "se_score_undefined", "on_score")])
cand$composite_score_variantA <- jaccard_composite(scores)
cand <- cand %>% arrange(desc(composite_score_variantA))

cat("\n  Top 15 Variant A candidates by composite score:\n")
top15 <- head(cand, 15)
for (i in seq_len(nrow(top15))) {
  r <- top15[i, ]
  cat(sprintf("    %-30s | score: %.4f | inds: %.3f | on: %.3f | n_ind: %.0f\n",
              substr(as.character(r$drug_name), 1, 30), r$composite_score_variantA,
              r$inds_score, r$on_score, r$n_indications))
}

out_a <- file.path(PROCESSED_DIR, "resko_variantA_no_se_candidates.csv")
write_csv(cand, out_a)
cat(sprintf("\n  Saved Variant A to %s\n", out_a))

# ============================================================
# VARIANT B: adds interactome pathway_score (built in next step)
# ============================================================
cat("\n[3/4] Variant B (interactome-supplemented) is built in a separate step -- see 04b_score_variantB_interactome.R\n")

cat("\n[4/4] Done\n")
cat("\nScoring complete\n")
cat("  Next: Run 04b_score_variantB_interactome.R, then 05_visualize_network_3d.R\n")
