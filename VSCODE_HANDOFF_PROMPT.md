# Handoff prompt: build a full-interactome interactive HTML network graph for RESKO-Faithful (Variant A & Variant B), styled like the Hetnet-eEF1A repo's network figure

Paste this to the VS Code agent working in this repo (`RESKO-Network`, branch `Faithful`).

---

## Status update — read this first

A prior pass already built a **minimal** interactive HTML graph:
`src/06_visualize_network_html.R` (R/visNetwork) produces
`resko_variantA_network_interactive.html` and
`resko_variantB_network_interactive.html`. Those are a small **hub-and-spoke**
graph: eEF1A1/eEF1A2 at the center, only the top-20-per-variant union (24
candidates) on a ring around it, one edge per candidate straight to the
center. They are explorable (hover, zoom, node search) but they do **not**
show the underlying protein interactome at all — no gene/protein nodes,
no protein-protein edges, nothing but "drug → target".

**This task is a different, bigger graph.** The user wants the RESKO
candidate network rendered in the same visual style as the sister
`Hetnet-eEF1A` repo's `results/hetnet_network_graph.png` — see the
description of that figure below — which shows the **actual interactome**
(genes/proteins as nodes, PPI edges from STRING/IntAct/Hetionet, a
force-directed/spring layout, NOT a ring) with significant compounds
overlaid as a distinct node shape/color, connected to the specific protein
target(s) they hit (not all wired to one central hub node). Build this as a
**new Python script**, since no Python HTML/interactive script currently
exists in this repo at all (the existing Python `05_visualize_network_3d.py`
only emits a static PNG — see below).

## Reference design: Hetnet-eEF1A's `hetnet_network_graph.png`

That figure (static PNG, matplotlib/networkx, no runnable script currently
committed for it in that repo either) shows:
- **Node types**: `EEF1A1`/`EEF1A2` (large, dark navy circles) · other
  interactome genes (small grey circles, sized/labeled by degree) ·
  significant compounds (red **squares**, distinct shape from genes)
- **Layout**: force-directed/spring (NOT a ring) — genes cluster by
  connectivity, compounds sit near the specific gene(s) they're evidenced
  to hit
- **Edges**: grey = ordinary interactome PPI edge; red = edge from a
  significant compound to its hit gene(s)
- **Legend**: three-entry categorical legend (EEF1A1/EEF1A2, other
  interactome gene, significant compound)
- Title states the finding in one line (e.g. *"Bonferroni-significant
  compounds sit on paths into the eEF1A neighborhood"*)

Reproduce this shape and encoding, adapted to RESKO's two-variant scoring
(color/highlight compounds by `composite_score_variantA` /
`composite_score_variantB` rather than by a p-value, since RESKO has no
significance test — see `METHODS_DEVIATIONS.md`), and make it an
**interactive HTML** file (Plotly or pyvis), not a static PNG.

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

- `data/raw/interactome_nodes.csv` (472 rows, one header) — protein nodes:
  `symbol`, `uniprot`, `degree`, `is_seed` (2 rows: EEF1A1/EEF1A2),
  `first_shell` (470 of 472 rows — **every non-seed node in this graph is a
  1-hop partner by construction**, see `METHODS_DEVIATIONS.md`),
  `n_seed_links`, `in_string`, `in_intact`, `in_hetionet`, `hiv_rtc_context`
- `data/raw/interactome_edges.csv` (15,714 rows) — `source`, `target`,
  `source_dbs`, `n_sources`, `tier` (`seed_incident`, 615 rows — edges
  touching EEF1A1/EEF1A2 directly — vs `partner_partner`, 15,099 rows —
  edges between two non-seed partners), `in_hetionet`, `string_score`,
  `intact_mi_score`, `n_records`. **This edge table is what a
  Hetnet-style force-directed graph should actually render** — the R port's
  ring graph never touches it.
- `data/raw/interactome_drug_target_links.csv` (228,722 rows) — drug→target
  links: `target`, `drug_id`, `drug_name`, `drug_type`, `max_stage`, `moa`,
  `moa_targets`, `moa_n_targets`, `sider_id`. Use this to draw a compound's
  edge(s) to its *actual* hit gene(s) in `interactome_nodes.csv`, not to a
  single central hub.
- `data/processed/candidate_drugs_no_se.csv` (415 candidates) — per-drug
  aggregation: `drug_id`, `drug_name`, `drug_type`, `max_stage`,
  `n_interactome_targets_hit`, `targets_hit` (semicolon- or list-encoded —
  check the actual delimiter before parsing), `moa`, `moa_n_targets`,
  `hits_first_shell`, `hits_eef1a_seed`
- `data/processed/resko_variantA_no_se_candidates.csv` and
  `resko_variantB_interactome_candidates.csv` — final scored/ranked tables.
  Full column list of the Variant B file (Variant A file is the same minus
  the B-specific columns): `drug_id, drug_name, drug_type, max_stage,
  n_interactome_targets_hit, targets_hit, moa, moa_n_targets,
  hits_first_shell, hits_eef1a_seed, n_indications, inds_score, on_score,
  se_score_undefined, composite_score_variantA, pathway_score_raw,
  pathway_score, composite_score_variantB, rank_variantA, rank_variantB,
  rank_shift_B_vs_A`
- `src/seed_compounds.py` — the seed registry (`SEED_COMPOUNDS`,
  `EXCLUDED_COMPOUNDS`), 13 seeds, only 4 SE-eligible (Molibresib,
  Plitidepsin, Didemnin_B, Metarrestin)

## What already exists for visualization (starting points, not the deliverable)

- `src/05_visualize_network_3d.py` (Python/matplotlib/networkx) →
  `resko_variantA_vs_variantB_network.png` — static, two-panel, **ring**
  layout of only the top-20-per-variant union (24 candidates) wired
  straight to a single eEF1A hub node. Does not touch the interactome
  edge table at all.
- `R/05_visualize_network_3d.R` (ggplot2/patchwork) — R port of the above,
  same ring design, `resko_variantA_vs_variantB_network_R.png`.
- `R/06_visualize_network_html.R` (visNetwork) — interactive HTML version of
  the same ring design: `resko_variantA_network_interactive.html`,
  `resko_variantB_network_interactive.html`. Explorable (hover/zoom/node
  search) but still the reduced 24-node hub-and-spoke shape, not the
  interactome.
- None of the above render `interactome_nodes.csv`/`interactome_edges.csv`
  as a graph. **The Hetnet-style figure this task asks for has never been
  built for RESKO** — you are not upgrading an existing full-interactome
  script, you are writing the first one.

## Your task

Write a new Python script, `src/06_visualize_network_html.py`, that renders
the **full RESKO interactome** (`interactome_nodes.csv` +
`interactome_edges.csv`, force-directed layout — e.g. `networkx.spring_layout`
or `kamada_kawai_layout` for node positions, then draw with Plotly for the
interactive HTML output) in the same visual grammar as Hetnet-eEF1A's
`hetnet_network_graph.png`:

- **EEF1A1/EEF1A2** — large, distinctly colored circular nodes (their
  `is_seed` rows)
- **other interactome genes** — smaller grey circles, sized by `degree`,
  labeled (at least on hover; static text labels for all 472 will be
  unreadable — consider labeling only nodes above a degree threshold, or
  making all labels hover-only)
- **candidate compounds** — a visually distinct node shape/color (square,
  matching the Hetnet figure's convention) for candidates above some
  ranking cutoff (e.g. top 20-30 by `composite_score_variantA` /
  `composite_score_variantB` — don't render all 415, it will be unreadable),
  colored/sized by that variant's composite score instead of Hetnet's
  p-value (RESKO has no significance test — never invent one)
- **edges**: ordinary interactome PPI edges in grey (from
  `interactome_edges.csv`); a compound's edge(s) to its real hit gene(s)
  (from `interactome_drug_target_links.csv`, filtered to genes present in
  `interactome_nodes.csv`) drawn in a distinct color, e.g. red — this
  mirrors Hetnet's highlighted "significant compound → target" edges,
  substituting "top-ranked by composite score" for "Bonferroni-significant"
- Produce **two HTML files**, one per variant (`resko_variantA_full_
  interactome.html`, `resko_variantB_full_interactome.html`), or a single
  page with a variant toggle — either is fine, but the two variants' node
  colorings must not be conflated in one static render
- Hover tooltips: gene symbol + degree for interactome nodes; drug name +
  composite score + rank + hit target(s) for compound nodes. Never label
  `se_score_undefined` as a real side-effect score — see
  `METHODS_DEVIATIONS.md` for what each score legitimately means before
  writing any tooltip text.
- A one-line title stating a finding, in the Hetnet figure's style — e.g.
  naming which candidates sit closest to the eEF1A hub under each variant.
