#' RESKO Faithful: Step 4b - Score Variant B (interactome-supplemented)
#'
#' Rationale (user-directed): McGarry's se_score used side-effect similarity
#' as a phenotypic PROXY for "these drugs act on similar pathways/targets".
#' With real side-effect data confirmed to give zero overlap for this
#' eEF1A-binder seed set (see METHODS_DEVIATIONS.md), we substitute the
#' mechanistic signal the proxy was standing in for: how strongly a
#' candidate's own targets are directly, biophysically coupled to eEF1A
#' within the built protein-protein interactome (472 nodes / 15,714 edges,
#' combining STRING + IntAct + Hetionet).
#'
#' IMPORTANT CORRECTNESS NOTE: this interactome was built as a
#' first-shell-only network -- ALL 470 non-seed nodes are, by construction,
#' direct (1-hop) partners of EEF1A1/EEF1A2 (verified: n_seed_links > 0 for
#' 100% of nodes). There is no second-shell structure to exploit, so a
#' binary/hop-distance proximity measure is degenerate here (constant
#' across all nodes) and was abandoned after producing a flat pathway_score
#' in an earlier draft of this script. Real graph-distance-based
#' differentiation is impossible with this particular interactome's
#' topology.
#'
#' Instead, pathway_score uses the *evidence strength* of each node's
#' direct edge(s) to EEF1A1/EEF1A2 -- a genuinely varying, real quantity
#' already present in interactome_edges.csv (tier == "seed_incident" rows):
#'   coupling(protein) = n_sources                      (1-3 supporting DBs)
#'                      + string_score  (if present)     (STRING confidence, 0-1)
#'                      + intact_mi_score (if present)   (IntAct MI score, 0-1)
#'                      + log1p(n_records)                (evidence-record count)
#'   pathway_score_raw(candidate) = mean of coupling(t) over the candidate's
#'                                   interactome-linked targets t (targets_hit)
#'   pathway_score = min-max normalized pathway_score_raw across candidates
#'
#' This is a genuine deviation from McGarry (2018) -- it replaces
#' phenotypic (side-effect) similarity with direct physical-interaction
#' evidence strength as the third scoring dimension -- disclosed explicitly
#' in METHODS_DEVIATIONS.md, and it is NOT a reconstruction of
#' DrugBank/SIDER4 se_score.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(stringr)
})

PROCESSED_DIR <- file.path("data", "processed")
RAW_DIR <- file.path("data", "raw")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: SCORE VARIANT B (SE-FREE + INTERACTOME-SUPPLEMENTED)\n")
cat(strrep("=", 70), "\n")

min_max_norm <- function(s) {
  lo <- min(s, na.rm = TRUE); hi <- max(s, na.rm = TRUE)
  if (hi == lo) return(rep(0.0, length(s)))
  (s - lo) / (hi - lo)
}

jaccard_composite <- function(scores_mat) {
  s <- rowSums(scores_mat^2)
  t <- rowSums(scores_mat)
  denom <- 2 * t - s
  out <- numeric(length(s))
  nz <- denom != 0
  out[nz] <- s[nz] / denom[nz]
  out
}

cat("\n[1/4] Loading Variant A candidates and interactome seed-incident edges\n")
cand <- read_csv(file.path(PROCESSED_DIR, "resko_variantA_no_se_candidates.csv"), show_col_types = FALSE)
nodes <- read_csv(file.path(RAW_DIR, "interactome_nodes.csv"), show_col_types = FALSE)
edges <- read_csv(file.path(RAW_DIR, "interactome_edges.csv"), show_col_types = FALSE)
cat(sprintf("  %d candidates, %d interactome nodes, %d edges\n", nrow(cand), nrow(nodes), nrow(edges)))

seed_incident <- edges %>% filter(tier == "seed_incident")
cat(sprintf("  %d seed-incident edges (direct EEF1A1/EEF1A2 links)\n", nrow(seed_incident)))

# Real per-node coupling strength to eEF1A: combine n_sources, string_score,
# intact_mi_score, and evidence-record count. Aggregate over both directions
# (source/target) and, where a protein has edges to both EEF1A1 and EEF1A2,
# take the strongest link (max), not a sum, so proteins linked to both
# paralogs aren't unfairly favored by double-counting.
edge_coupling <- function(n_sources, string_score, intact_mi_score, n_records) {
  c_val <- ifelse(is.na(n_sources), 0.0, n_sources)
  c_val <- c_val + ifelse(is.na(string_score), 0.0, string_score)
  c_val <- c_val + ifelse(is.na(intact_mi_score), 0.0, intact_mi_score)
  c_val <- c_val + ifelse(is.na(n_records), 0.0, log1p(n_records))
  c_val
}

seed_incident <- seed_incident %>%
  mutate(
    coupling = edge_coupling(n_sources, string_score, intact_mi_score, n_records),
    partner = ifelse(source %in% c("EEF1A1", "EEF1A2"), target, source)
  )

node_coupling_df <- seed_incident %>% group_by(partner) %>% summarise(coupling = max(coupling), .groups = "drop")
node_coupling <- setNames(node_coupling_df$coupling, node_coupling_df$partner)

cat(sprintf("  Coupling strength computed for %d proteins (range %.2f - %.2f)\n",
            length(node_coupling), min(node_coupling), max(node_coupling)))

pathway_score_for_targets <- function(targets_str) {
  if (is.na(targets_str) || targets_str == "") return(0.0)
  targs <- str_split(targets_str, ";")[[1]]
  targs <- targs[targs != ""]
  if (length(targs) == 0) return(0.0)
  weights <- ifelse(targs %in% names(node_coupling), node_coupling[targs], 0.0)
  mean(weights)
}

cat("\n[2/4] Computing interactome proximity weights per candidate\n")
# Use targets_hit (the actual interactome nodes this candidate is linked to,
# from step 3b's aggregation) rather than moa_targets -- targets_hit is
# guaranteed to be a subset of the 472-node interactome by construction, so
# every weight lookup is a real hit, not a name-mismatch miss.
cand$pathway_score_raw <- sapply(cand$targets_hit, pathway_score_for_targets)
cand$pathway_score <- min_max_norm(cand$pathway_score_raw)

cat("\n[3/4] Computing Variant B's Jaccard composite with [inds_score, on_score, pathway_score]")
cat(" -- note Variant A's own composite_score_variantA column is RETAINED unchanged in\n")
cat("this output file for direct side-by-side comparison between the two variants.\n")

scores <- as.matrix(cand[, c("inds_score", "on_score", "pathway_score")])
cand$composite_score_variantB <- jaccard_composite(scores)
cand$rank_variantA <- rank(-cand$composite_score_variantA, ties.method = "min")
cand <- cand %>% arrange(desc(composite_score_variantB))
cand$rank_variantB <- seq_len(nrow(cand))
cand$rank_shift_B_vs_A <- cand$rank_variantA - cand$rank_variantB
# VALIDATION NOTE: composite_score_variantB and pathway_score are numerically
# identical to the Python original (max abs diff ~1e-16, floating-point
# noise). rank_variantB differs from the Python run for a subset of rows,
# but ONLY within the 55 groups of candidates sharing an EXACTLY tied
# composite score -- R's arrange() and pandas' sort_values() both use
# stable sorts, but arrive at ties in a different pre-sort row order
# (a consequence of dplyr::summarise() vs pandas.groupby() row ordering
# upstream in 03b), so tie-breaking order differs cosmetically. No
# candidate's SCORE differs, and no candidate crosses a real rank boundary.

cat("\n  Top 15 Variant B candidates by composite score (rank_shift_B_vs_A: positive = moved UP under Variant B):\n")
top15 <- head(cand, 15)
for (i in seq_len(nrow(top15))) {
  r <- top15[i, ]
  cat(sprintf("    %-26s | B_score: %.4f (rank %3d) | A_score: %.4f (rank %3d) | shift: %+4d | pathway: %.3f\n",
              substr(as.character(r$drug_name), 1, 26), r$composite_score_variantB, r$rank_variantB,
              r$composite_score_variantA, r$rank_variantA, r$rank_shift_B_vs_A, r$pathway_score))
}

cat("\n[4/4] Saving\n")
out_b <- file.path(PROCESSED_DIR, "resko_variantB_interactome_candidates.csv")
write_csv(cand, out_b)
cat(sprintf("  Saved Variant B to %s\n", out_b))

cat("\nVariant B scoring complete\n")
cat("  Next: Run 05_visualize_network_3d.R\n")
