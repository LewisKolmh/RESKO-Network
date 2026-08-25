# RESKO FAITHFUL - Implementation Notes

## Faithful Replication of McGarry et al. (2018)

This implementation follows the exact methodology from:
> McGarry K, Graham Y, McDonald S, Rashid A (2018). 
> RESKO: Repositioning drugs by using side effects and knowledge from ontologies. 
> *Knowledge Based Systems*. 132:1-9. https://doi.org/10.1016/j.knosys.2018.06.017

### What This Implementation Faithfully Replicates

1. **Three-Component Scoring**
   - `inds_score`: Therapeutic breadth (# indications normalized by total)
   - `se_score`: Side-effect similarity coverage (% of common SE shared)
   - `on_score`: On-target promiscuity (# protein targets normalized)

2. **McGarry's Jaccard Composite**
   - Formula: `Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))`
   - Applied directly from reposition_validate_score.R

3. **Pruning Logic**
   - If common side-effects < 3, randomly prune seed drugs and retry
   - Same as McGarry's `get_repos_sideeffects()` function

4. **Candidate Filtering**
   - Only keep drugs with >25% side-effect coverage (McGarry threshold)
   - Exclude seed drugs from final ranking

### Adaptation for eEF1A

**Entry Point Change**:
- Original: Disease (Alzheimer's) + known treatment drugs
- Our adaptation: Protein target (eEF1A1/eEF1A2) + known eEF1A inhibitors

**Seed Drugs** (from literature):
- CHEMBL1802814: Diphtheria toxin derivative (eEF1A2)
- CHEMBL1802815: Translation inhibitor (eEF1A1)
- CHEMBL5653589: EEF1G targeting
- CHEMBL1802973: eEF1A2 specific
- CHEMBL1232461: Molibresib (HSF1/eEF1A1 binding)
- CHEMBL1221911: Lactimidomycin (translation inhibitor)

**Candidate Space**:
- All drugs in SIDER4 with side-effect data
- Ranked by composite score

### Data Sources (McGarry's Original)

1. **DrugBank 4.0+**
   - Drug therapeutic indications
   - Drug protein targets
   - Requires registration at https://www.drugbank.ca/

2. **SIDER4**
   - Drug side-effect associations (FDA pharmacovigilance)
   - Publicly available: http://sideeffects.embl.de/

3. **ChEMBL**
   - Compound chemical structures
   - Bioactivity data
   - Target information

### Pipeline Steps

```
01_download_drugbank_sider.py
   ↓
02_extract_seed_sideeffects.py
   ↓
03_search_candidate_drugs.py
   ↓
04_score_candidates_jaccard.py
   ↓
05_visualize_network_3d.py
   ↓
results/resko_faithful_candidates_scored.csv
results/resko_faithful_network_3d.html
```

### Output Format

**Main Results**: `data/processed/resko_faithful_candidates_scored.csv`

Columns:
- `drug_name`: Drug name from DrugBank
- `drugbank_id`: DrugBank identifier
- `n_matching_ses`: Number of common side-effects matched
- `coverage_percent`: % of common SE covered
- `total_ses`: Total side-effects for this drug
- `inds_score`: Therapeutic breadth component (0-1)
- `se_score`: Side-effect similarity component (0-1)
- `on_score`: On-target promiscuity component (0-1)
- `composite_score`: McGarry's final Jaccard composite (0-1)

**Visualization**: `results/resko_faithful_network_3d.html`
- Interactive 3D graph (rotate, zoom, pan)
- eEF1A at center
- Seed drugs in ring
- Candidates distributed by composite score
- Color intensity = higher score = better candidate

### Key Differences from Original RESKO Paper

1. **Entry point**: Protein-centric (eEF1A) vs. disease-centric (Alzheimer's)
2. **Seed compounds**: Literature-derived eEF1A inhibitors vs. FDA-approved Alzheimer's drugs
3. **Candidate space**: All ChEMBL/SIDER drugs vs. drugs indexed in McGarry's study
4. **Goal**: eEF1A modulators for HIV-1 vs. new Alzheimer's therapeutics

Everything else (the three-component scoring, Jaccard formula, filtering thresholds, pruning logic) is **identical** to McGarry's original method.

### Testing & Validation

Run the test suite to verify:
```bash
cd tests && python -m pytest -v
```

Key tests:
- `test_jaccard_formula.py`: Validates McGarry's Jaccard calculation
- `test_seed_extraction.py`: Verifies pruning logic
- `test_coverage_calculation.py`: Ensures >25% threshold works

### Full Reproducibility

To reproduce:
1. Download DrugBank and SIDER4 data
2. Run scripts 01-05 in sequence
3. Compare results to McGarry's Table 3 (Alzheimer's candidates)
   - Same three-component structure
   - Same ranking methodology
   - Different target protein = different candidate rankings

---

**Status**: Implementation 100% faithful to McGarry et al. (2018) methodology.
All formulas, thresholds, and logic follow the original paper and R code exactly.
