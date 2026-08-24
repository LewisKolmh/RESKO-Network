# resko-network

Python port of the `RESKO-Network` R pipeline: eEF1A allosteric-inhibitor
screening: relation-aware evidence ranking and chemical-similarity screening
against SIDER, ported from the R pipeline in `scripts/`.

The port is not a translation. Three defects found in the R original are fixed at
source, and each fix is pinned by a regression test so it cannot silently revert.
The diagnosis for each is in `CODE_REVIEW.md`.

| Fix | R original | Here |
|---|---|---|
| **S1-1** Censored measurements | `potency_score()` (A15:101-110) read only the numeric value, so `IC50 >= 10 µM` scored as a weak positive | `ranking.potency_score()` reads `standard_relation`; censored records score 0 |
| **S1-2** Scaffold artefacts | Similarity reported with no scaffold context, so shared-core hits looked like pharmacology | `similarity.screen()` returns `scaffold_tanimoto` and `scaffold_artefact` per pair |
| **S2-1** Cache-key drift | Filenames built inline under two padding conventions across revisions, leaving orphans with a different response schema | `chembl.cache_key()` is the only constructor; `chembl.audit_cache_dir()` detects strays |

## Setup

```bash
conda env create -f environment.yml
conda activate resko-network
pip install -e ".[dev]"
pytest                    # 70 tests (63 unit + 7 doctests), ~4 s
```

Then in VS Code: `mv dot-vscode .vscode`, and select the `resko-network`
interpreter (Cmd-Shift-P → *Python: Select Interpreter*). Tests appear in the
Testing panel automatically.

## Reproducing the ranking correction

The headline result — molibresib falling from rank #1 to #6 — from your existing
R outputs, no refetch needed:

```python
import pandas as pd
from resko_network.ranking import rank_candidates, ScoringMode

evidence = pd.read_csv("results/A12B_compound_protein_evidence_corrected.csv")

legacy    = rank_candidates(evidence, ScoringMode.legacy())   # reproduces R
corrected = rank_candidates(evidence, ScoringMode())          # relation-aware

print(legacy.loc[0, "display_name"])      # MOLIBRESIB
print(corrected.loc[0, "display_name"])   # CHEMBL1802814
```

`ScoringMode.legacy()` exists precisely so the before/after is reproducible
rather than asserted. 4 of 23 evidence records change score; all 4 are
molibresib's, and no other compound is affected.

### How faithful is `legacy()`? (validated, not assumed)

`legacy()` reproduces A15's **relation-blind behaviour** — the defect — not its
exact composite arithmetic. Measured against `A15_translation_network_ranking.csv`:

| Check | Result |
|---|---|
| Top-4 order identical | yes |
| Rank 1 is molibresib | yes |
| Spearman vs A15 | 0.943 |
| Bit-exact composite scores | **no** |

The two disagreements are honest and bounded. A15 assigns each compound a single
potency *bucket* from its strongest record; the port sums potency across records
(`total_potency`), which is a different aggregation. And ranks 5–6 swap because
A15 has a genuine three-way tie at `total_score = 12` (CHEMBL5653589,
CHEMBL3752910, CHEMBL1802973) that it leaves in input order, while the port
breaks ties by `total_potency`.

So: use `legacy()` to demonstrate that relation-blindness puts molibresib first.
Do not cite it as a bit-exact reimplementation of A15. `port_validation.csv`
records these checks.

## Auditing a cache directory

```python
from resko_network.chembl import audit_cache_dir

audit = audit_cache_dir("results/A17D_api_cache")
paged = audit[audit.paged]
print(paged[~paged.canonical_name.astype(bool)])              # orphaned files
print(paged.groupby("canonical_name").n_fields_union.describe())  # schema divergence
```

On the committed cache this reports 36 files — 15 canonical, 15 orphaned, 6
un-paged `molecule_*.json`. Record counts match exactly across every same-offset
pair (12,292 each side), but the orphaned files carry 47 keys on every record
while the canonical ones carry 28–36, omitting `standard_upper_value` entirely.
That is the field a censored measurement's ceiling lives in.

## Screening with scaffold awareness

```python
from resko_network.similarity import prepare, screen

queries = prepare(candidates, "molecule_chembl_id", "canonical_smiles")
library = prepare(sider_drugs, "drug_name", "canonical_smiles")

hits = screen(queries, library)
hits[hits.scaffold_artefact]      # pairs explained by a shared Murcko scaffold
hits[~hits.scaffold_artefact]     # what remains once echoes are removed
```

## Module layout

```
src/resko_network/
├── config.py       every threshold and weight, in one place
├── chembl.py       paged API client, canonical cache keys, cache audit
├── sider.py        SIDER loading, boundary-anchored seed matching, coverage
├── similarity.py   Morgan fingerprints, Murcko scaffolds, artefact flagging
└── ranking.py      relation-aware scoring, candidate ranking
```

`docs/R_TO_PYTHON_MAP.md` maps each of the 42 R scripts to its Python home,
including the ones deliberately not ported and why.

## Design notes

**Config is centralised.** In the R original the 0.30 candidate threshold lived
at A17C:1142 and the potency ladder at A15:101-110, several thousand lines apart.
Changing a parameter should mean editing one line of `config.py`.

**Corrections are switchable.** `ScoringMode.legacy()` and the
`relation_aware` / `demote_proteomics` keyword arguments exist so that "the fix
changed the result" is a claim you can run, not one you have to trust.

**Caches are inputs, not artefacts.** API responses belong in a gitignored cache
directory keyed by one canonical function, with the API version recorded
alongside. They should not be committed — the R repo tracks 42 MB of them.

**Absence is recorded, not dropped.** `sider.match_drug_names()` emits a row for
an unmatched seed rather than omitting it, so a coverage denominator can always
be reconstructed. Category labels like "Didemnin analogues" are flagged and
excluded from that denominator, since they can never match a database row.
