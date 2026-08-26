# RESKO-Faithful: R port

A complete R re-implementation of the Python pipeline in `../src/`,
covering the current two-variant scoring design (see
`../METHODS_DEVIATIONS.md`). Every step was validated against the
corresponding Python-generated output before being treated as trustworthy;
results below.

## Requirements

R packages (installed into the `r` conda environment used for this port):
`arrow`, `dplyr`, `readr`, `jsonlite`, `stringr`, `tidyr`, `ggplot2`,
`patchwork`, `igraph`, `visNetwork`, `htmlwidgets`, `httr2`.

## Pipeline steps and validation status

Run from the repository root (`resko-faithful/`), in order:

| # | Script | Produces | Validation vs Python |
|---|---|---|---|
| 0 | `R/seed_compounds.R` | `SEED_COMPOUNDS`, `EXCLUDED_COMPOUNDS` (in-memory registry, sourced by later steps) | Exact match: 13 seeds, 4 SE-eligible, 9 binding-evidence-only, 5 excluded |
| 1 | `R/01_download_sider_faers.R` | `data/raw/sider_all_se.parquet`, `data/raw/faers_adr.parquet` | Exact match: 309,849 SIDER4 pairs (1,430 drugs / 6,061 side effects); identical FAERS term counts per seed |
| 2 | `R/02_extract_seed_sideeffects.R` | `data/processed/seed_sideeffects.json`, `seed_coverage.csv`, `seed_se_analysis_molibresib_plitidepsin.png` | Same seed-pruning logic and 0-common-side-effect finding across the 2 data-bearing seeds |
| 3 | `R/03_search_candidate_drugs.R` | `data/processed/candidate_drugs.csv` | Reproduces the same structural finding as Python (0 common SE → 0 candidates by construction); resolves gracefully instead of crashing |
| 3b | `R/03b_search_candidates_interactome.R` | `data/processed/candidate_drugs_no_se.csv` | **Byte-identical** to Python across all 415 rows / 10 columns |
| 4 | `R/04_score_candidates_jaccard.R` | `data/processed/resko_variantA_no_se_candidates.csv` | **Numerically identical** to Python (max abs diff 0.0) on `inds_score`, `on_score`, `composite_score_variantA` |
| 4b | `R/04b_score_variantB_interactome.R` | `data/processed/resko_variantB_interactome_candidates.csv` | **Numerically identical** to Python (diff ~1e-16, float noise) on `pathway_score`, `composite_score_variantB`. `rank_variantB` differs from Python only within 55 exactly-tied-score groups (stable-sort tie-break order, no real score/rank-boundary difference — see comment in the script) |
| 5 | `R/05_visualize_network_3d.R` | `resko_variantA_vs_variantB_network_R.png` | Same shared ring layout, same 24-candidate union, same solid/dashed interactome-support edge encoding, reimplemented with ggplot2 + patchwork (R has no direct networkx/matplotlib equivalent) |
| 6 | `R/06_visualize_network_html.R` | `resko_variantA_network_interactive.html`, `resko_variantB_network_interactive.html` | Same layout/coloring/edge logic as step 5, made explorable via visNetwork (hover tooltips, zoom, node search) — fulfils `../VSCODE_HANDOFF_PROMPT.md`'s interactive-HTML request natively in R |

## Notes on deviations from the Python originals

- Step 03's R port does not crash on the empty-candidate case (the Python
  original raises `KeyError: 'coverage_percent'` trying to sort an empty
  DataFrame) — it writes a header-only CSV and prints the same structural
  disclosure instead.
- Step 05/06 use ggplot2/visNetwork idioms rather than literal
  networkx/matplotlib/pyvis translations, since those Python libraries have
  no R equivalent; the node positions, coloring rule, and edge-encoding
  logic are carried over exactly.
- `se_score_undefined` is fixed at 0.0 throughout (both languages) and is
  deliberately omitted from step 6's tooltips rather than mislabeled as a
  real side-effect score — see `../METHODS_DEVIATIONS.md`.
