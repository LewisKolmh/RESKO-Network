"""
RESKO Faithful: Step 2 - Extract Common Side-Effects from Seed Drugs

Input: Seed eEF1A binders defined in seed_compounds.SEED_COMPOUNDS
       + SIDER4 drug-side effect data
       + openFDA/FAERS data (fallback for compounds absent from SIDER4 —
         see seed_compounds.py docstring for why this fallback exists, and
         why it replaces the originally-planned DrugBank-ADR fallback)

4 seeds (Ternatin-4, Narciclasine, Nannocystin Ax, BE-43547A2) are
binding-evidence-only: they have never been dosed in a human, so no
side-effect source can ever have data for them. They are EXCLUDED from the
intersection below (per seed_compounds.binding_evidence_only_seeds()) but
retained in SEED_COMPOUNDS for network diagrams / structural justification.

Output: Common side-effects (intersection across SE-eligible seed drugs,
        pooling SIDER4 + FAERS per compound)
        + Pruning logic (if <3 common SE, randomly prune seeds and retry)
        + A per-seed coverage report (data/processed/seed_coverage.csv) so
          missing data is visible rather than silently assumed away.

This follows McGarry's get_repos_sideeffects() function, extended with the
two-tier SE source described above.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json
import sys

sys.path.insert(0, str(Path(__file__).parent))
from seed_compounds import SEED_COMPOUNDS, se_eligible_seeds, binding_evidence_only_seeds

DATA_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

print("="*70)
print("RESKO FAITHFUL: EXTRACT COMMON SIDE-EFFECTS FROM SEED DRUGS")
print("="*70)
print(f"\nAll registered seed compounds ({len(SEED_COMPOUNDS)}), from seed_compounds.py:")
for cid, meta in SEED_COMPOUNDS.items():
    print(f"  {cid:16s} drugbank_id={str(meta['drugbank_id']):10s} "
          f"evidence={meta['evidence']:26s} sider4_expected={meta['sider4_expected']} "
          f"has_human_exposure={meta.get('has_human_exposure')}")

binding_only = binding_evidence_only_seeds()
se_eligible = se_eligible_seeds()
print(f"\n  SE-eligible seeds (participate in intersection below): {len(se_eligible)}")
print(f"  Binding-evidence-only seeds (never dosed in humans, EXCLUDED from "
      f"intersection): {len(binding_only)} -> {list(binding_only)}")

# Load SIDER data
print("\n[1/5] Loading SIDER4 data")
print("-" * 70)

sider_file = DATA_DIR / "sider_all_se.parquet"
if not sider_file.exists():
    print(f"ERROR: {sider_file} not found. Run 01_download_drugbank_sider.py first.")
    exit(1)

sider_df = pd.read_parquet(sider_file)
print(f"  Loaded {len(sider_df)} drug-side effect associations")
print(f"  Unique drugs: {sider_df['drugbank_id'].nunique()}")
print(f"  Unique side-effects: {sider_df['side_effect_name'].nunique()}")

# Load openFDA/FAERS fallback data (produced by 01_download_drugbank_sider.py)
faers_file = DATA_DIR / "faers_adr.parquet"
if faers_file.exists():
    faers_df = pd.read_parquet(faers_file)
    print(f"  FAERS fallback: {len(faers_df)} seed-reaction associations "
          f"({faers_df['seed_id'].nunique() if len(faers_df) else 0} seeds)")
else:
    faers_df = pd.DataFrame(columns=['seed_id', 'drugbank_id', 'side_effect_name'])
    print("  WARNING: data/raw/faers_adr.parquet not found — "
          "fallback source unavailable, run 01_download_drugbank_sider.py")

# Extract seed drug side-effects — SIDER4 first (by drugbank_id), FAERS
# fallback (by seed_id, since several seeds have no drugbank_id), and
# RECORD which source (or neither) actually supplied data per seed.
# Only SE-eligible seeds are processed here; binding-evidence-only seeds are
# reported separately below and never enter the intersection.
print("\n[2/5] Extracting side-effects for each SE-eligible seed drug (with coverage check)")
print("-" * 70)

seed_side_effects = {}
coverage_rows = []
for cid, meta in se_eligible.items():
    dbid = meta["drugbank_id"]
    sider_ses = []
    if dbid is not None:
        sider_ses = sider_df[sider_df['drugbank_id'] == dbid]['side_effect_name'].unique().tolist()
    faers_ses = faers_df[faers_df['seed_id'] == cid]['side_effect_name'].unique().tolist()
    pooled = sorted(set(sider_ses) | set(faers_ses))
    seed_side_effects[cid] = pooled

    if sider_ses:
        source = "SIDER4"
    elif faers_ses:
        source = "FAERS_fallback"
    else:
        source = "NONE — dropped from intersection"
    coverage_rows.append({
        "seed_id": cid, "drugbank_id": dbid, "n_sider4": len(sider_ses),
        "n_faers": len(faers_ses), "n_pooled": len(pooled), "source_used": source,
    })
    print(f"  {cid:16s} SIDER4={len(sider_ses):3d}  FAERS={len(faers_ses):3d}  "
          f"-> pooled={len(pooled):3d}  [{source}]")

# Record the binding-evidence-only seeds in the SAME report, clearly marked,
# so the coverage file documents every registered seed's fate, not just the
# ones that entered the intersection.
for cid, meta in binding_only.items():
    coverage_rows.append({
        "seed_id": cid, "drugbank_id": meta["drugbank_id"], "n_sider4": 0,
        "n_faers": 0, "n_pooled": 0,
        "source_used": "N/A — binding-evidence-only (never dosed in humans)",
    })
    print(f"  {cid:16s} [excluded from SE step — never dosed in humans; "
          f"binding-evidence-only]")

coverage_df = pd.DataFrame(coverage_rows)
coverage_df.to_csv(PROCESSED_DIR / "seed_coverage.csv", index=False)
print(f"\n  ✓ Coverage report saved to data/processed/seed_coverage.csv")

# Among SE-eligible seeds, some may still come back with zero SE data from
# either source (e.g. a FAERS query with no hits) — exclude those from the
# intersection explicitly rather than let them silently zero out the whole
# common-SE set, and say so.
uncovered = [cid for cid, ses in seed_side_effects.items() if len(ses) == 0]
if uncovered:
    print(f"\n  ⚠ {len(uncovered)} SE-eligible seed(s) returned NO side-effect "
          f"data from either source and are excluded from the intersection: "
          f"{uncovered}")
    for cid in uncovered:
        del seed_side_effects[cid]

# Find common side-effects (intersection)
print("\n[3/5] Finding common side-effects across all seed drugs")
print("-" * 70)

def get_common_sideeffects(seed_dict, min_common=3, max_iterations=4):
    """McGarry's pruning logic: find common SE, prune if <3 common"""
    current_drugs = list(seed_dict.keys())
    iteration = 0
    
    while iteration < max_iterations:
        iteration += 1
        print(f"  Iteration {iteration}: Testing {len(current_drugs)} drugs")
        
        # Find intersection
        se_sets = [set(seed_dict[drug_id]) for drug_id in current_drugs]
        common_se = set.intersection(*se_sets) if se_sets else set()
        
        print(f"    Common side-effects found: {len(common_se)}")
        
        if len(common_se) >= min_common:
            print(f"    ✓ Sufficient common SE found!")
            return list(common_se), current_drugs
        
        if len(current_drugs) <= 2:
            print(f"    ✗ Cannot prune further (only {len(current_drugs)} drugs left)")
            # Return what we have
            return list(common_se), current_drugs
        
        # McGarry's pruning: randomly sample half the drugs
        n_prune = max(2, len(current_drugs) // 2)
        current_drugs = np.random.choice(current_drugs, n_prune, replace=False).tolist()
        print(f"    Pruning to {n_prune} drugs: {current_drugs}")
    
    print(f"  ✗ Could not find ≥{min_common} common SE after {max_iterations} iterations")
    return list(common_se), current_drugs

common_ses, final_seed_drugs = get_common_sideeffects(seed_side_effects)

print(f"\n  Final seed drugs used: {len(final_seed_drugs)}")
for drug_id in final_seed_drugs:
    print(f"    - {drug_id}")

print(f"  Final common side-effects: {len(common_ses)}")
for i, se in enumerate(common_ses[:10], 1):
    print(f"    {i}. {se}")
if len(common_ses) > 10:
    print(f"    ... and {len(common_ses) - 10} more")

# Save results
print("\n[4/5] Saving results")
print("-" * 70)

results = {
    "seed_drugs": final_seed_drugs,
    "common_sideeffects": common_ses,
    "n_common_se": len(common_ses),
    "n_seed_drugs": len(final_seed_drugs),
}

with open(PROCESSED_DIR / "seed_sideeffects.json", 'w') as f:
    json.dump(results, f, indent=2)

seed_df = pd.DataFrame({
    "side_effect": common_ses
})
seed_df.to_csv(PROCESSED_DIR / "seed_sideeffects.csv", index=False)

print(f"  ✓ Saved to data/processed/seed_sideeffects.json")
print(f"  ✓ Saved to data/processed/seed_sideeffects.csv")

print("\n[5/5] Data-source summary (for methods write-up)")
print("-" * 70)
se_rows = coverage_df[~coverage_df["source_used"].str.startswith("N/A")]
n_sider4 = (se_rows["n_sider4"] > 0).sum()
n_faers_only = ((se_rows["n_sider4"] == 0) & (se_rows["n_faers"] > 0)).sum()
n_none = (se_rows["n_pooled"] == 0).sum()
print(f"  SE-eligible seeds with SIDER4 coverage:  {n_sider4}")
print(f"  SE-eligible seeds using FAERS fallback:  {n_faers_only}")
print(f"  SE-eligible seeds with NO SE data found: {n_none}")
print(f"  Binding-evidence-only seeds (excluded, never dosed in humans): "
      f"{len(binding_only)}")
if n_faers_only > 0:
    print("  NOTE: openFDA/FAERS fallback is a deviation from McGarry's "
          "SIDER4-only method (and from the originally-planned DrugBank-ADR "
          "fallback, unavailable due to DrugBank's academic-download pause) "
          "— disclose in methods section, including FAERS's voluntary-report "
          "bias profile.")

print("\n✓ Seed side-effects extraction complete")
print("  Next: Run 03_search_candidate_drugs.py")
