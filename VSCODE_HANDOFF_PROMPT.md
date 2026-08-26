# Handoff prompt: build the interactive HTML network graph for RESKO-Faithful (Variant A & Variant B)

Paste this to the VS Code agent working in this repo (`RESKO-Network`, branch `Faithful`).

---

## Context: what this repo is

RESKO-Faithful re-implements Ken McGarry's RESKO drug-repositioning method
(McGarry et al. 2018, *Knowledge-Based Systems*, three-component Jaccard
composite score: `inds_score` + `se_score` + `on_score`) applied to a new
target: finding drug candidates that may modulate eEF1A1/eEF1A2 (elongation
factor 1-alpha), the protein this project's seed compounds are confirmed
direct binders of. Full background: `README.md`, `IMPLEMENTATION_NOTES.md`.

**Critical: McGarry's original `se_score` (side-effect similarity) could not
be computed for this seed set.** Real querying of SIDER4, FAERS, and
ClinicalTrials.gov confirmed zero shared side-effect terms across the
SE-eligible seeds — a genuine data ceiling, not a retrieval gap. Full
explanation: `METHODS_DEVIATIONS.md` (read this before touching any scoring
code — it is the single source of truth for what each score column means).

Because of that gap, the pipeline produces **two parallel scored variants**
instead of one. Your task is to build the HTML network graph(s) for **both**,
clearly labeled as distinct.

## The two variants — exact differentiating names

| | **Variant A** | **Variant B** |
|---|---|---|
| Full name | RESKO minus side-effects | RESKO minus side-effects, interactome-supplemented |
| Composite score column | `composite_score_variantA` | `composite_score_variantB` |
| Components (Jaccard input) | `inds_score`, `se_score_undefined` (fixed 0.0), `on_score` | `inds_score`, `on_score`, `pathway_score` |
| Script that computes it | `src/04_score_candidates_jaccard.py` | `src/04b_score_variantB_interactome.py` |
| Output CSV | `data/processed/resko_variantA_no_se_candidates.csv` | `data/processed/resko_variantB_interactome_candidates.csv` (also retains `composite_score_variantA`, plus `rank_variantA`, `rank_variantB`, `rank_shift_B_vs_A`) |

**`pathway_score`** (Variant B only) is a real, continuously-varying measure
of direct physical-interaction evidence between a candidate's targets and
eEF1A1/eEF1A2, computed from `interactome_edges.csv`'s `seed_incident` rows:
`coupling(protein) = n_sources + string_score + intact_mi_score +
log1p(n_records)`, averaged over the candidate's linked targets, then
min-max normalized. It is **not** a side-effect proxy — see
`METHODS_DEVIATIONS.md` for why an earlier binary-weight draft was
degenerate and was replaced with this formula.

## What already exists (do not recompute — load and consume)

- `data/raw/interactome_nodes.csv` (472 rows) — protein nodes: `symbol`,
  `degree`, `is_seed`, `first_shell`, `n_seed_links`, `in_string`,
  `in_intact`, `in_hetionet`
- `data/raw/interactome_edges.csv` (15,714 rows) — `source`, `target`,
  `source_dbs`, `n_sources`, `tier` (`seed_incident` vs `partner_partner`),
  `string_score`, `intact_mi_score`, `n_records`
- `data/raw/interactome_drug_target_links.csv` — drug→target links with
  `drug_id`, `drug_name`, `drug_type`, `max_stage`, `moa`, `moa_targets`,
  `moa_n_targets`
- `data/processed/candidate_drugs_no_se.csv` (415 candidates) — per-drug
  aggregation: `targets_hit`, `n_interactome_targets_hit`, `hits_eef1a_seed`,
  `hits_first_shell`
- `data/processed/resko_variantA_no_se_candidates.csv` and
  `resko_variantB_interactome_candidates.csv` — final scored/ranked tables
  (see column list above)
- `src/seed_compounds.py` — the seed registry (`SEED_COMPOUNDS`,
  `EXCLUDED_COMPOUNDS`), 13 seeds, only 4 SE-eligible (Molibresib,
  Plitidepsin, Didemnin_B, Metarrestin)

## What already exists for visualization (starting point, not the deliverable)

`src/05_visualize_network_3d.py` currently renders a **static PNG**
(`resko_variantA_vs_variantB_network.png`, saved as a project artifact) —
two side-by-side matplotlib/networkx panels sharing one node-position
layout (ring layout, eEF1A1/eEF1A2 at center, top-20-per-variant candidates
on the ring, radial labels), colored per-panel by that panel's own composite
score, with Variant B additionally drawing solid edges for
interactome-evidenced links (`pathway_score > 0`) vs. dashed for
unsupported. This script does **not** produce HTML — an earlier
now-superseded draft of this same file used `plotly.graph_objects` for an
interactive 3D graph, but it predates the two-variant split and reads a
filename (`resko_faithful_candidates_scored.csv`) that no longer exists.

## Your task

Build a new script (suggest `src/06_visualize_network_html.py`) that produces
**interactive HTML network graph(s)** — e.g. via Plotly or pyvis — covering
both variants, reusing the reference layout/encoding logic already validated
in the PNG (same node set, same "candidate colored by its own score" idea,
same solid/dashed interactome-support distinction for Variant B) but making
it explorable (hover tooltips with drug name/score/target list, zoom, toggle
between variants or render both in one page). Read `resko_variantA_vs_
variantB_network.png`'s generating script (`src/05_visualize_network_3d.py`)
for the exact node-position/color logic to carry over, and
`METHODS_DEVIATIONS.md` for what each score legitimately means so tooltips
don't misrepresent them (e.g. never label `se_score_undefined` as a real
side-effect score in a tooltip).
