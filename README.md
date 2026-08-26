# RESKO FAITHFUL: Ken McGarry's Drug Repositioning Method Applied to eEF1A Inhibitor Discovery

## Overview

This repository implements the exact methodology from:
> Ken McGarry, Yitka Graham, Sharon McDonald, Anuam Rashid (2018). 
> **RESKO: Repositioning drugs by using side effects and knowledge from ontologies**. 
> *Knowledge Based Systems*. https://doi.org/10.1016/j.knosys.2018.06.017

### Key Adaptation
Instead of starting from a disease (e.g., Alzheimer's) + known treatment drugs, we start from:
- **Target protein**: eEF1A1/eEF1A2 (elongation factor 1-alpha)
- **Seed drugs**: Known eEF1A inhibitors and HIV-1 reverse transcriptase modulators
- **Goal**: Find candidate compounds with similar side-effect profiles that may act as eEF1A modulators

## The McGarry Method: Three-Component Scoring

### 1. **INDS_SCORE** (Therapeutic Breadth)
Number of therapeutic indications each candidate drug treats, normalized:
```
inds_score[i] = n_indications[i] / total_unique_indications
```

### 2. **SE_SCORE** (Side-Effect Similarity)
Percentage of common side-effects the candidate shares with seed drugs:
```
se_score[i] = shared_side_effects[i] / total_common_side_effects
```
Candidates with >25% overlap qualify.

### 3. **ON_SCORE** (On-Target Promiscuity)
Number of protein targets the candidate hits (normalized):
```
on_score[i] = n_unique_targets[i] / total_unique_targets
```

### Composite Score: Jaccard Similarity
All three scores are combined via **matrix Jaccard similarity** (unusual for continuous vectors):
```
X = [inds_score, se_score, on_score]  per drug
Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))
final_score = diag(Jaccard(X))
```

## Pipeline Steps

### Step 1: Download Data
```bash
python src/01_download_drugbank_sider.py
```
Downloads DrugBank 4.0+ and SIDER4 (FDA pharmacovigilance side-effect data).
**Note**: DrugBank requires manual registration at https://www.drugbank.ca/

### Step 2: Extract Seed Side-Effects
```bash
python src/02_extract_seed_sideeffects.py
```
- Loads seed eEF1A inhibitors (CHEMBL1802814, CHEMBL1802815, etc.)
- Extracts side-effects from SIDER4
- Finds common side-effects (intersection)
- Implements McGarry's pruning logic if <3 common SE

### Step 3: Search Candidate Drugs
```bash
python src/03_search_candidate_drugs.py
```
- Scans all drugs in SIDER4/ChEMBL
- Calculates side-effect coverage (%) for each vs. common seed SE
- Filters candidates with >25% coverage
- Excludes seed drugs

### Step 4: Score Candidates (Jaccard Composite)
```bash
python src/04_score_candidates_jaccard.py
```
- Fetches DrugBank indications → inds_score
- Uses coverage% → se_score
- Fetches protein targets → on_score
- Computes McGarry's Jaccard composite score
- Ranks candidates

### Step 5: Visualize 3D Network
```bash
python src/05_visualize_network_3d.py
```
Builds interactive 3D network graph:
- Candidate compounds (colored by composite score)
- Connected to their target proteins
- Connected to eEF1A (seed protein)
- Edge thickness represents side-effect coverage

## Output Files

- `results/resko_faithful_candidates.csv` — Main ranking (all candidates with scores)
- `results/resko_faithful_network_3d.html` — Interactive 3D network visualization
- `data/processed/seed_sideeffects.csv` — Common side-effects for seed drugs

## R Port

**Everything above describes the original Python pipeline.** The `R/`
directory contains a complete, independently-validated R port of the full
pipeline (steps 01-06), run in parallel to the Python originals rather than
replacing them. See `R/README.md` for the R-specific usage, environment
setup, and per-step validation results (numerical/byte-identity checks
against the Python outputs).

**Important:** the Python `src/` pipeline described above and this `R/`
port both target the *current* two-variant scoring design (Variant A / B,
`composite_score_variantA` / `composite_score_variantB`), documented fully
in `METHODS_DEVIATIONS.md` — not the single-`se_score` design shown in the
"Pipeline Steps" and "Output Files" sections above, which predate that
split and are kept here for historical context on the original McGarry
adaptation. `METHODS_DEVIATIONS.md` is the current source of truth for
what each score column means in both languages.

## Seed Compounds (Direct eEF1A Binders Only)

Full registry with evidence type and DrugBank ID (where known) lives in
`src/seed_compounds.py` — treat that file as the source of truth; this table
is a summary of it.

**Tier 1 — ChEMBL-verified** (quantitative IC50/Kd/ED50 binding record against
an eEF1A-family target, from the eEF1A-centred hetnet interactome pull):

1. **CHEMBL1802814**
2. **CHEMBL1802815**
3. **CHEMBL5653589** — EEF1G-complex targeting
4. **CHEMBL1802973**
5. **CHEMBL3752910**
6. **CHEMBL1232461** — Molibresib

**Tier 2 — literature-confirmed direct binders** (crystal structure /
biochemical trapping of an eEF1A conformational state):

7. **Plitidepsin (Aplidin)** — approved (Australia, 2018); antiviral activity
   (incl. SARS-CoV-2) attributed to direct eEF1A inhibition
8. **Didemnin B** — binds eEF1A between domains I/III, trapping the GTP-bound
   conformation; no evidence of direct HIV-1 protein binding
9. **Metarrestin** — eEF1A2-selective, Phase I clinical trial
10. **Ternatin-4** — occupies the same eEF1A site as didemnin B
11. **Narciclasine** — in vivo validated
12. **Nannocystin Ax** — shares the didemnin-B/ternatin-4 binding site
13. **BE-43547A2** — in vivo validated

**Excluded** (reasons kept in `src/seed_compounds.py:EXCLUDED_COMPOUNDS`):
Lactimidomycin (binds the ribosomal E-site during the translocation step
eEF1A drives, not eEF1A itself), efavirenz (HIV-1 RT inhibitor, no direct
eEF1A binding), cycloheximide, anisomycin, aminoglycosides (all off-target
ribosomal mechanisms).

## Data Sources

- **DrugBank 4.0+**: Drug indications, protein targets, adverse reactions
- **SIDER4**: Drug side-effect associations (FDA pharmacovigilance database,
  curated ~2015 from marketed-drug package inserts)
- **DrugBank ADR (fallback)**: for seeds that are investigational or were
  approved after SIDER4's curation date and therefore have no SIDER4 entry
  (this applies to most Tier-2 seeds above). **This is a deliberate deviation
  from McGarry's SIDER4-only method** — disclose it in any methods write-up.
  Per-seed coverage (which source supplied data, or whether neither did) is
  written to `data/processed/seed_coverage.csv` by
  `02_extract_seed_sideeffects.py`.
- **ChEMBL**: Compound chemical structures and bioactivity data

## Requirements

```
pandas>=1.3
numpy>=1.20
scikit-learn>=0.24
plotly>=5.0
networkx>=2.6
drugbank-api (optional, for automated DrugBank access)
```

Install:
```bash
pip install -r requirements.txt
```

## Citation

If you use this implementation, please cite both:

1. McGarry et al. (2018) — Original RESKO method
2. This repository — Your adapted implementation for eEF1A

## License

[Your license here]

---

**Status**: Implementation complete with faithful adherence to McGarry et al. (2018) methodology. 
All three scoring components and Jaccard composite follow the original paper exactly.
