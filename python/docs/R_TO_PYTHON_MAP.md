# R script → Python module map

All 42 scripts in `RESKO-Network/scripts/` accounted for. "Status" is one of:

- **ported** — behaviour available in the Python package now
- **planned** — stage exists in the pipeline design, not yet written
- **superseded** — the script's job disappears in the new structure, reason given
- **not ported** — deliberately dropped, reason given

## Stage A2–A3 — SIDER seed coverage

| R script | Python | Status | Notes |
|---|---|---|---|
| `A2_inspect_sider.R` | `sider.load_drug_names`, `sider.load_side_effects` | ported | |
| `A2_seed_lookup.R` | `sider.match_drug_names` | ported | Unanchored `str_detect` (A2:35) replaced with boundary-anchored matching — S3-1 |
| `A3_check_expanded_seed_set.R` | `sider.match_drug_names`, `sider.coverage_report` | ported | Category labels excluded from the coverage denominator — S3-1 |

## Stage A4–A7 — Interaction and annotation network

| R script | Python | Status | Notes |
|---|---|---|---|
| `A4_get_eef1a_string_network.R` | `network.fetch_string_partners` | planned | Config in `config.StringConfig` |
| `A5_create_neo4j_files.R` | — | superseded | Graph lives in `networkx`; no CSV staging layer needed. Export via `network.to_neo4j_csv` if Neo4j is still wanted |
| `A6_reactome_enrichment.R` | `enrichment.reactome` | planned | |
| `A6_create_pathway_edges.R` | `network.add_pathway_edges` | planned | Node/edge CSV pairs collapse into one graph-building call |
| `A6_create_pathway_nodes.R` | `network.add_pathway_edges` | planned | |
| `A6_create_protein_master_list.R` | `network.protein_registry` | planned | |
| `A7_GO_enrichment.R` | `enrichment.gene_ontology` | planned | |
| `A7_create_GO_edges.R` | `network.add_go_edges` | planned | |
| `A7_create_GO_nodes.R` | `network.add_go_edges` | planned | |

## Stage A8–A12 — ChEMBL bioactivity

| R script | Python | Status | Notes |
|---|---|---|---|
| `A8_test_chembl_connection.R` | — | not ported | A connectivity smoke test; `ChemblClient` raises on failure with the URL attached |
| `A8_get_chembl_targets.R` | `chembl.ChemblClient.targets` | planned | |
| `A8_get_drug_activities.R` | `chembl.ChemblClient.activities` | ported | Single canonical `cache_key` — S2-1 |
| `A8_export_protein_list.R` | `network.protein_registry` | planned | |
| `A8_live_network.R` | `viz.live_network` | planned | |
| `A10_filter_chembl_activities.R` | `activities.quality_filter` | planned | |
| `A11_annotate_candidate_molecules.R` | `chembl.ChemblClient.molecules` | planned | |
| `A12_create_compound_protein_edges.R` | `network.add_compound_edges` | planned | |
| `A12B_correct_compound_protein_evidence.R` | `ranking.score_evidence_table` | ported | Relation and assay-format awareness added — S1-1 |
| `A12B_correct_compound_protein_evidence_v3.R` | `ranking.score_evidence_table` | superseded | The `_v3` suffix is version control; git handles it |

Activity flattening (the `A9_chembl_activities_complete.csv` output) is
`chembl.activities_to_frame`, which retains `standard_relation` and
`standard_upper_value` explicitly — those are the fields S1-1 and S2-1 turn on.

## Stage A13–A15 — Provenance and ranking

| R script | Python | Status | Notes |
|---|---|---|---|
| `A13_live_network_with_compounds.R` | `viz.live_network` | planned | |
| `A13B_live_network_corrected_evidence.R` | `viz.live_network` | superseded | Same figure, corrected input — a parameter, not a script |
| `A14_curate_assay_document_provenance.R` | `provenance.annotate` | planned | This stage is the reason the molibresib record was traceable; keep it |
| `A14B_finalise_provenance_outputs.R` | `provenance.annotate` | superseded | |
| `A15_rank_candidates.R` | `ranking.rank_candidates` | **ported** | The S1-1 fix. `ScoringMode.legacy()` reproduces the original |
| `A15B_live_network_ranked.R` | `viz.live_network` | planned | |

## Stage A16–A17 — Chemical similarity and SIDER bridge

| R script | Python | Status | Notes |
|---|---|---|---|
| `A16A_retrieve_reference_structures.R` | `chembl.ChemblClient.molecules` | planned | |
| `A16B_analyse_chemical_similarity.R` | `similarity.screen` | ported | |
| `A16_environment.yml` | `environment.yml` | ported | One environment for the whole project, not per-stage |
| `A17A_prepare_sider_reference_library.R` | `sider.load_drug_names` | ported | |
| `A17B_build_sider_identifier_crosswalk.R` | `sider.build_crosswalk` | planned | InChIKey-based, replacing name matching entirely |
| `A17C_screen_sider_drugs_by_chemical_similarity.R` | `similarity.screen` | **ported** | The S1-2 fix: scaffold columns added. ~1,400 lines → ~150 |
| `A17D_evaluate_sider_candidates.R` | `chembl` + `ranking` | ported | Cache logic to `chembl`, scoring to `ranking` |
| `A17E_correct_sider_candidate_provenance.R` | `provenance.annotate` | planned | |

## Stage A18 — Commercial availability

| R script | Python | Status | Notes |
|---|---|---|---|
| `A18A_build_compound_detail_manifest.R` | `sourcing.build_manifest` | planned | |
| `A18B_enrich_pubchem_commercial_availability.R` | `sourcing.pubchem_vendors` | planned | |
| `A18C_reconcile_commercial_product_identities.R` | `sourcing.reconcile` | planned | |
| `A18D_prepare_live_network_compound_data.R` | `viz.live_network` | planned | |
| `A18E_live_network_with_compound_details.R` | `viz.live_network` | planned | |
| `A18F_supplier_search_preperation_script.R` | `sourcing.export_search_sheet` | planned | Note: filename typo `preperation` in the original |
| `A18G_supplier_results_import_script.R` | `sourcing.import_results` | planned | |

## Not ported, and why

| Item | Reason |
|---|---|
| `A8_test_chembl_connection.R` | Connectivity assertions belong in the client, not a separate script |
| `A5_create_neo4j_files.R` | Neo4j CSV staging is unnecessary with an in-process graph; an exporter can be added if the Neo4j browser is still in use |
| `A15B_script_numbered.txt` | Not a script — a line-numbered copy of `A15B`, committed by accident |

## Structural changes worth knowing

**42 scripts → 5 modules (so far).** The numbered-script layout encodes execution
order in filenames, which means a corrected re-run becomes a new script
(`A12B`, `_v3`, `A13B`). Here a correction is a parameter (`ScoringMode`), so
there is one code path and both variants are reproducible from it.

**`_previous_<timestamp>` outputs disappear.** Those 26 files exist because the
scripts overwrite their outputs. Git already versions outputs; the port writes to
a caller-specified path and leaves history to git.

**Node/edge CSV pairs collapse.** Six `A6`/`A7` scripts produce paired node and
edge tables for the same entity type. In `networkx` that is one
`add_nodes_from` / `add_edges_from` pair inside a single function.

**Line-count comparison.** `A17C` is ~1,400 lines; `similarity.py` is ~150 and
does more (scaffold decomposition, artefact flagging). The difference is mostly
inline column-selection blocks and repeated summary-table construction, which
pandas handles in one expression.
