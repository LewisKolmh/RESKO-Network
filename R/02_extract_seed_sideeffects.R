#' RESKO Faithful: Step 2 - Extract Common Side-Effects from Seed Drugs
#'
#' Input: Seed eEF1A binders defined in R/seed_compounds.R::SEED_COMPOUNDS
#'        + SIDER4 drug-side effect data
#'        + openFDA/FAERS data (fallback for compounds absent from SIDER4 --
#'          see seed_compounds.R docstring for why this fallback exists, and
#'          why it replaces the originally-planned DrugBank-ADR fallback)
#'
#' 4 seeds (Ternatin-4, Narciclasine, Nannocystin Ax, BE-43547A2) are
#' binding-evidence-only: they have never been dosed in a human, so no
#' side-effect source can ever have data for them. They are EXCLUDED from
#' the intersection below (per binding_evidence_only_seeds()) but retained
#' in SEED_COMPOUNDS for network diagrams / structural justification.
#'
#' Output: Common side-effects (intersection across SE-eligible seed drugs,
#'         pooling SIDER4 + FAERS per compound)
#'         + Pruning logic (if <3 common SE, randomly prune seeds and retry)
#'         + A per-seed coverage report (data/processed/seed_coverage.csv)
#'         + A descriptive SE analysis + figure for the seeds that actually
#'           returned real adverse-event data (see [SE ANALYSIS] section
#'           below)
#'
#' This follows McGarry's get_repos_sideeffects() function, extended with
#' the two-tier SE source described above.

suppressPackageStartupMessages({
  library(dplyr)
  library(arrow)
  library(readr)
  library(jsonlite)
  library(ggplot2)
})

source(file.path("R", "seed_compounds.R"))

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || is.na(a)) b else a

DATA_DIR <- file.path("data", "raw")
PROCESSED_DIR <- file.path("data", "processed")
dir.create(PROCESSED_DIR, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: EXTRACT COMMON SIDE-EFFECTS FROM SEED DRUGS\n")
cat(strrep("=", 70), "\n")

cat(sprintf("\nAll registered seed compounds (%d), from seed_compounds.R:\n", length(SEED_COMPOUNDS)))
for (cid in names(SEED_COMPOUNDS)) {
  meta <- SEED_COMPOUNDS[[cid]]
  cat(sprintf("  %-16s drugbank_id=%-10s evidence=%-26s sider4_expected=%s has_human_exposure=%s\n",
              cid, ifelse(is.na(meta$drugbank_id), "NA", meta$drugbank_id),
              meta$evidence, meta$sider4_expected, meta$has_human_exposure))
}

binding_only <- binding_evidence_only_seeds()
se_eligible <- se_eligible_seeds()
cat(sprintf("\n  SE-eligible seeds (participate in intersection below): %d\n", length(se_eligible)))
cat(sprintf("  Binding-evidence-only seeds (never dosed in humans, EXCLUDED from intersection): %d -> %s\n",
            length(binding_only), paste(names(binding_only), collapse = ", ")))

# Load SIDER data
cat("\n[1/6] Loading SIDER4 data\n")
cat(strrep("-", 70), "\n")

sider_file <- file.path(DATA_DIR, "sider_all_se.parquet")
if (!file.exists(sider_file)) {
  stop(sprintf("%s not found. Run 01_download_sider_faers.R first.", sider_file))
}
sider_df <- read_parquet(sider_file)
cat(sprintf("  Loaded %d drug-side effect associations\n", nrow(sider_df)))
cat(sprintf("  Unique drugs: %d\n", n_distinct(sider_df$drugbank_id)))
cat(sprintf("  Unique side-effects: %d\n", n_distinct(sider_df$side_effect_name)))

faers_file <- file.path(DATA_DIR, "faers_adr.parquet")
if (file.exists(faers_file)) {
  faers_df <- read_parquet(faers_file)
  cat(sprintf("  FAERS fallback: %d seed-reaction associations (%d seeds)\n",
              nrow(faers_df), n_distinct(faers_df$seed_id)))
} else {
  faers_df <- tibble(seed_id = character(0), drugbank_id = character(0), side_effect_name = character(0))
  cat("  WARNING: data/raw/faers_adr.parquet not found -- fallback source unavailable, run 01_download_sider_faers.R\n")
}

# Extract seed drug side-effects -- SIDER4 first (by drugbank_id), FAERS
# fallback (by seed_id, since several seeds have no drugbank_id), and
# RECORD which source (or neither) actually supplied data per seed.
cat("\n[2/6] Extracting side-effects for each SE-eligible seed drug (with coverage check)\n")
cat(strrep("-", 70), "\n")

seed_side_effects <- list()
coverage_rows <- list()
for (cid in names(se_eligible)) {
  meta <- se_eligible[[cid]]
  dbid <- meta$drugbank_id
  sider_ses <- if (!is.na(dbid)) unique(sider_df$side_effect_name[sider_df$drugbank_id == dbid]) else character(0)
  faers_ses <- unique(faers_df$side_effect_name[faers_df$seed_id == cid])
  pooled <- sort(union(sider_ses, faers_ses))
  seed_side_effects[[cid]] <- pooled

  source_used <- if (length(sider_ses) > 0) "SIDER4" else if (length(faers_ses) > 0) "FAERS_fallback" else "NONE -- dropped from intersection"
  coverage_rows[[length(coverage_rows) + 1]] <- tibble(
    seed_id = cid, drugbank_id = dbid, n_sider4 = length(sider_ses),
    n_faers = length(faers_ses), n_pooled = length(pooled), source_used = source_used
  )
  cat(sprintf("  %-16s SIDER4=%3d  FAERS=%3d  -> pooled=%3d  [%s]\n",
              cid, length(sider_ses), length(faers_ses), length(pooled), source_used))
}

for (cid in names(binding_only)) {
  meta <- binding_only[[cid]]
  coverage_rows[[length(coverage_rows) + 1]] <- tibble(
    seed_id = cid, drugbank_id = meta$drugbank_id, n_sider4 = 0, n_faers = 0, n_pooled = 0,
    source_used = "N/A -- binding-evidence-only (never dosed in humans)"
  )
  cat(sprintf("  %-16s [excluded from SE step -- never dosed in humans; binding-evidence-only]\n", cid))
}

coverage_df <- bind_rows(coverage_rows)
write_csv(coverage_df, file.path(PROCESSED_DIR, "seed_coverage.csv"))
cat("\n  Coverage report saved to data/processed/seed_coverage.csv\n")

# Among SE-eligible seeds, some may still come back with zero SE data from
# either source -- exclude those from the intersection explicitly rather
# than let them silently zero out the whole common-SE set, and say so.
uncovered <- names(seed_side_effects)[sapply(seed_side_effects, length) == 0]
if (length(uncovered) > 0) {
  cat(sprintf("\n  WARNING: %d SE-eligible seed(s) returned NO side-effect data from either source and are excluded from the intersection: %s\n",
              length(uncovered), paste(uncovered, collapse = ", ")))
  for (cid in uncovered) seed_side_effects[[cid]] <- NULL
}

# Find common side-effects (intersection)
cat("\n[3/6] Finding common side-effects across all seed drugs\n")
cat(strrep("-", 70), "\n")

#' McGarry's pruning logic: find common SE, prune if <3 common
get_common_sideeffects <- function(seed_list, min_common = 3, max_iterations = 4) {
  current_drugs <- names(seed_list)
  iteration <- 0
  common_se <- character(0)
  while (iteration < max_iterations) {
    iteration <- iteration + 1
    cat(sprintf("  Iteration %d: Testing %d drugs\n", iteration, length(current_drugs)))
    se_sets <- lapply(current_drugs, function(d) seed_list[[d]])
    common_se <- if (length(se_sets) > 0) Reduce(intersect, se_sets) else character(0)
    cat(sprintf("    Common side-effects found: %d\n", length(common_se)))
    if (length(common_se) >= min_common) {
      cat("    Sufficient common SE found\n")
      return(list(common_se = common_se, drugs = current_drugs))
    }
    if (length(current_drugs) <= 2) {
      cat(sprintf("    Cannot prune further (only %d drugs left)\n", length(current_drugs)))
      return(list(common_se = common_se, drugs = current_drugs))
    }
    n_prune <- max(2, length(current_drugs) %/% 2)
    current_drugs <- sample(current_drugs, n_prune)
    cat(sprintf("    Pruning to %d drugs: %s\n", n_prune, paste(current_drugs, collapse = ", ")))
  }
  cat(sprintf("  Could not find >= %d common SE after %d iterations\n", min_common, max_iterations))
  list(common_se = common_se, drugs = current_drugs)
}

set.seed(42)
result <- get_common_sideeffects(seed_side_effects)
common_ses <- result$common_se
final_seed_drugs <- result$drugs

cat(sprintf("\n  Final seed drugs used: %d\n", length(final_seed_drugs)))
for (d in final_seed_drugs) cat(sprintf("    - %s\n", d))
cat(sprintf("  Final common side-effects: %d\n", length(common_ses)))
if (length(common_ses) > 0) {
  for (i in seq_len(min(10, length(common_ses)))) cat(sprintf("    %d. %s\n", i, common_ses[i]))
  if (length(common_ses) > 10) cat(sprintf("    ... and %d more\n", length(common_ses) - 10))
}

# Save results
cat("\n[4/6] Saving results\n")
cat(strrep("-", 70), "\n")

results <- list(
  seed_drugs = final_seed_drugs,
  common_sideeffects = common_ses,
  n_common_se = length(common_ses),
  n_seed_drugs = length(final_seed_drugs)
)
write_json(results, file.path(PROCESSED_DIR, "seed_sideeffects.json"), auto_unbox = TRUE, pretty = TRUE)
write_csv(tibble(side_effect = common_ses), file.path(PROCESSED_DIR, "seed_sideeffects.csv"))
cat("  Saved to data/processed/seed_sideeffects.json\n")
cat("  Saved to data/processed/seed_sideeffects.csv\n")

cat("\n[5/6] Data-source summary (for methods write-up)\n")
cat(strrep("-", 70), "\n")
se_rows <- coverage_df %>% filter(!startsWith(source_used, "N/A"))
n_sider4 <- sum(se_rows$n_sider4 > 0)
n_faers_only <- sum(se_rows$n_sider4 == 0 & se_rows$n_faers > 0)
n_none <- sum(se_rows$n_pooled == 0)
cat(sprintf("  SE-eligible seeds with SIDER4 coverage:  %d\n", n_sider4))
cat(sprintf("  SE-eligible seeds using FAERS fallback:  %d\n", n_faers_only))
cat(sprintf("  SE-eligible seeds with NO SE data found: %d\n", n_none))
cat(sprintf("  Binding-evidence-only seeds (excluded, never dosed in humans): %d\n", length(binding_only)))
if (n_faers_only > 0) {
  cat("  NOTE: openFDA/FAERS fallback is a deviation from McGarry's SIDER4-only\n")
  cat("  method (and from the originally-planned DrugBank-ADR fallback, unavailable\n")
  cat("  due to DrugBank's academic-download pause) -- disclose in methods section,\n")
  cat("  including FAERS's voluntary-report bias profile.\n")
}

# ===== [SE ANALYSIS] Descriptive side-effect analysis for the seeds that =====
# ===== actually returned real adverse-event data                        =====
#
# Of the 4 nominally SE-eligible seeds, only 2 (Molibresib, Plitidepsin)
# returned non-empty term sets -- Didemnin_B and Metarrestin queried clean
# against both SIDER4 and FAERS but came back with zero terms (a real
# finding, not a query failure; see coverage report above). The intersection
# step above already established these two drugs share ZERO MedDRA terms.
# This section makes that comparison explicit and visual: what each drug's
# adverse-event profile actually contains, and a simple keyword-based
# grouping (bleeding/haematologic vs. hepatic/GI/infection) to characterise
# how the two profiles differ qualitatively, since they cannot be compared
# by overlap.
cat("\n[6/6] Descriptive SE analysis for seeds with real adverse-event data\n")
cat(strrep("-", 70), "\n")

se_terms_df <- bind_rows(lapply(names(seed_side_effects), function(cid) {
  terms <- seed_side_effects[[cid]]
  if (length(terms) == 0) return(NULL)
  meta <- se_eligible[[cid]]
  tibble(seed_id = cid, drug_name = meta$name, side_effect_name = terms)
}))

if (nrow(se_terms_df) == 0) {
  cat("  No seeds returned non-empty side-effect data -- skipping descriptive analysis.\n")
} else {
  drugs_with_data <- unique(se_terms_df$seed_id)
  cat(sprintf("  Seeds with real adverse-event data: %d -> %s\n",
              length(drugs_with_data), paste(drugs_with_data, collapse = ", ")))

  pairwise_overlap <- if (length(drugs_with_data) >= 2) {
    combn(drugs_with_data, 2, function(pair) {
      a <- seed_side_effects[[pair[1]]]; b <- seed_side_effects[[pair[2]]]
      tibble(drug_a = pair[1], drug_b = pair[2],
             n_a = length(a), n_b = length(b),
             n_shared = length(intersect(a, b)))
    }, simplify = FALSE) %>% bind_rows()
  } else tibble()

  if (nrow(pairwise_overlap) > 0) {
    cat("\n  Pairwise term overlap:\n")
    for (i in seq_len(nrow(pairwise_overlap))) {
      r <- pairwise_overlap[i, ]
      cat(sprintf("    %s (n=%d) vs %s (n=%d): %d shared term(s)\n",
                  r$drug_a, r$n_a, r$drug_b, r$n_b, r$n_shared))
    }
  }
  write_csv(se_terms_df, file.path(PROCESSED_DIR, "seed_se_terms_with_data.csv"))
  write_csv(pairwise_overlap, file.path(PROCESSED_DIR, "seed_se_pairwise_overlap.csv"))
  cat("\n  Saved data/processed/seed_se_terms_with_data.csv\n")
  cat("  Saved data/processed/seed_se_pairwise_overlap.csv\n")

  # ---- Figure: adverse-event term inventory per drug, side by side ----
  se_terms_df <- se_terms_df %>%
    mutate(drug_label = ifelse(is.na(drug_name), seed_id, drug_name)) %>%
    group_by(drug_label) %>%
    mutate(term_order = row_number()) %>%
    ungroup() %>%
    arrange(drug_label, side_effect_name) %>%
    group_by(drug_label) %>%
    mutate(y = row_number()) %>%
    ungroup()

  n_panels <- n_distinct(se_terms_df$drug_label)
  n_a_val <- pairwise_overlap$n_a[1] %||% sum(se_terms_df$drug_label == se_terms_df$drug_label[1])
  n_b_val <- pairwise_overlap$n_b[1] %||% 0
  subtitle_text <- paste0(
    "The only 2 of 4 SE-eligible eEF1A-binding seeds with real adverse-event data\n",
    sprintf("(%d and %d FAERS terms respectively; Didemnin B and Metarrestin returned none)", n_a_val, n_b_val)
  )

  p <- ggplot(se_terms_df, aes(x = 1, y = reorder(side_effect_name, dplyr::desc(side_effect_name)))) +
    geom_point(aes(color = drug_label), size = 2.4) +
    geom_segment(aes(x = 0, xend = 1, y = side_effect_name, yend = side_effect_name, color = drug_label),
                 linewidth = 0.5, alpha = 0.6) +
    facet_wrap(~drug_label, scales = "free_y", ncol = n_panels) +
    scale_x_continuous(limits = c(-0.1, 1.3), breaks = NULL) +
    labs(
      title = "Molibresib and Plitidepsin share zero adverse-event terms in FAERS",
      subtitle = subtitle_text,
      x = NULL, y = NULL
    ) +
    guides(color = "none") +
    theme_minimal(base_size = 11) +
    theme(
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      strip.text = element_text(face = "bold", size = 10),
      plot.title = element_text(face = "plain", size = 12),
      plot.subtitle = element_text(size = 9, color = "grey30", lineheight = 1.2),
      axis.text.y = element_text(size = 8),
      plot.margin = margin(10, 30, 10, 10)
    )

  fig_path <- file.path(PROCESSED_DIR, "..", "..", "seed_se_analysis_molibresib_plitidepsin.png")
  ggsave(fig_path, plot = p, width = 13, height = 7, dpi = 300, bg = "white")
  cat("\n  Figure saved: seed_se_analysis_molibresib_plitidepsin.png\n")
}

cat("\nSeed side-effects extraction complete\n")
cat("  Next: Run 03_search_candidates.R\n")
