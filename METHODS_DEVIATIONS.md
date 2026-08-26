# Deviations from McGarry's original RESKO method

This pipeline runs **two variants**, both of which deviate from McGarry
(2018)'s original 3-component Jaccard composite (`inds_score` + `se_score` +
`on_score`) in the same way and then differ from each other in one further
respect. This document is the single source of truth for what changed, why,
and what each variant's output file/column actually means — read this before
comparing scores across variants or against McGarry's published method.

## Why se_score is undefined for both variants

McGarry's `se_score` component measures side-effect (adverse-event) profile
similarity between a seed drug and a candidate, sourced from SIDER4/FAERS/
DrugBank pharmacovigilance data. For this project's seed set (eEF1A-binding
compounds identified from literature), real querying of SIDER4 (1,430 drugs,
309,849 pairs), FAERS (synonym-pooled across generic and brand names), and
ClinicalTrials.gov's structured adverse-event section confirmed:

- Only 4 of 13 seeds have ANY clinical development history at all
  (Molibresib, Plitidepsin, Didemnin_B, Metarrestin) — the other 9 are
  research/binding-assay compounds with no ChEMBL preferred name, synonym,
  or clinical-phase record.
- Of those 4, only 2 (Molibresib, Plitidepsin) return any adverse-event
  terms in SIDER4/FAERS at all.
- Pooling every tested synonym/brand name for those 2 still produces **zero**
  shared MedDRA terms between them.
- Didemnin_B has no ClinicalTrials.gov entry (its trials predate the
  registry, est. 2000). Metarrestin's one registered trial is active with no
  posted results section.

This is a genuine data ceiling for this seed set specifically — not a
retrieval failure, and not fixable by sourcing DrugBank instead of SIDER4/
FAERS (the same structural gap — preclinical/early-trial compounds lacking
real-world exposure data — applies regardless of source). Per user decision,
`se_score` is fixed at `0.0` (stored as the column `se_score_undefined`) for
every candidate in both variants below, rather than silently dropped from
the pipeline or backfilled with a placeholder constant.

## Variant A — `resko_variantA_no_se_candidates.csv`

Composite score (`composite_score_variantA`) = McGarry's Jaccard formula
over `[inds_score, se_score_undefined, on_score]`, where:
- `inds_score` = min-max normalized count of distinct ChEMBL indication
  terms (`mesh_heading`/`efo_term`) for the candidate — real ChEMBL
  `drug_indication` data, not a placeholder.
- `on_score` = min-max normalized `moa_n_targets` — the number of distinct
  proteins in the candidate's ChEMBL-recorded mechanism-of-action set — a
  direct promiscuity measure, also real ChEMBL data.
- `se_score_undefined` = 0.0 for all rows (see above).

Candidate pool: drugs with ChEMBL-confirmed bioactivity against a protein in
the eEF1A first-shell interactome (`interactome_drug_target_links.csv`,
416 drugs before exclusions), independent of any side-effect coverage
filter (the original `03_search_candidate_drugs.py` filtered by SIDER4
coverage, which would yield zero candidates given the se_score finding
above).

## Variant B — `resko_variantB_interactome_candidates.csv`

Adds a third real component, `pathway_score`, replacing `se_score` as the
third Jaccard dimension: `composite_score_variantB` = Jaccard over
`[inds_score, on_score, pathway_score]`.

`pathway_score` is **not** a reconstruction of side-effect similarity. It
measures the direct, evidenced physical-interaction strength between the
candidate's own targets and eEF1A1/eEF1A2 within the built interactome
(472 nodes / 15,714 edges, STRING + IntAct + Hetionet):

```
coupling(protein) = n_sources + string_score + intact_mi_score + log1p(n_records)
pathway_score_raw(candidate) = mean(coupling(t) for t in candidate's interactome-linked targets)
pathway_score = min-max normalize(pathway_score_raw) across candidates
```

An earlier draft used a binary "is this protein a direct eEF1A partner"
weight, which was degenerate: **all 472 nodes in this interactome are
first-shell (1-hop) partners of EEF1A1/EEF1A2 by construction** — there is
no second-shell structure in this network to distinguish "close" from "far"
nodes by hop count. The coupling-strength formula above uses the real,
continuously-varying edge evidence attributes instead, which is why Variant
B produces genuinely different rankings from Variant A (see
`rank_shift_B_vs_A` column) rather than a constant offset.

`resko_variantB_interactome_candidates.csv` retains `composite_score_variantA`
unchanged alongside `composite_score_variantB`, plus `rank_variantA`,
`rank_variantB`, and `rank_shift_B_vs_A` (positive = candidate ranked higher
under Variant B) for direct side-by-side comparison — this is the
recommended file to use for cross-variant analysis; the plain Variant A CSV
lacks the rank-shift context.

## What is genuinely comparable to McGarry's original ranking, and what is not

Neither variant should be described as "RESKO with McGarry's method" without
qualification — both drop the se_score dimension McGarry used and Variant B
substitutes a mechanistically different (though also literature-grounded)
third term. They ARE directly comparable to *each other*, since both reuse
the same `inds_score`/`on_score` computation and the same underlying Jaccard
formula, and the network figure
(`resko_variantA_vs_variantB_network.png`) is built for exactly that
side-by-side reading — same node layout in both panels, colored by that
panel's own score.
