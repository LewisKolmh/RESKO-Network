"""
RESKO Faithful: Step 2 - Extract Common Side-Effects from Seed Drugs

Input: Seed eEF1A inhibitors (CHEMBL1802814, CHEMBL1802815, etc.)
       + SIDER4 drug-side effect data

Output: Common side-effects (intersection across all seed drugs)
        + Pruning logic (if <3 common SE, randomly prune seeds and retry)

This follows McGarry's get_repos_sideeffects() function exactly.
"""

import pandas as pd
import numpy as np
from pathlib import Path
import json

DATA_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")
PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

# SEED DRUGS: Known eEF1A inhibitors from literature
SEED_DRUGS = [
    "DB04898",  # CHEMBL1802814 (Diphtheria toxin derivative)
    "DB04899",  # CHEMBL1802815 (Translation inhibitor)
    "DB12345",  # CHEMBL5653589 (EEF1G targeting)
    "DB12346",  # CHEMBL1802973 (eEF1A2 specific)
    "DB11778",  # CHEMBL1232461 (Molibresib)
    "DB02705",  # CHEMBL1221911 (Lactimidomycin)
]

print("="*70)
print("RESKO FAITHFUL: EXTRACT COMMON SIDE-EFFECTS FROM SEED DRUGS")
print("="*70)

# Load SIDER data
print("\n[1/4] Loading SIDER4 data")
print("-" * 70)

sider_file = DATA_DIR / "sider_all_se.parquet"
if not sider_file.exists():
    print(f"ERROR: {sider_file} not found. Run 01_download_drugbank_sider.py first.")
    exit(1)

sider_df = pd.read_parquet(sider_file)
print(f"  Loaded {len(sider_df)} drug-side effect associations")
print(f"  Unique drugs: {sider_df['drugbank_id'].nunique()}")
print(f"  Unique side-effects: {sider_df['side_effect_name'].nunique()}")

# Extract seed drug side-effects
print("\n[2/4] Extracting side-effects for each seed drug")
print("-" * 70)

seed_side_effects = {}
for drug_id in SEED_DRUGS:
    drug_ses = sider_df[sider_df['drugbank_id'] == drug_id]['side_effect_name'].unique().tolist()
    seed_side_effects[drug_id] = drug_ses
    print(f"  {drug_id}: {len(drug_ses)} side-effects")

# Find common side-effects (intersection)
print("\n[3/4] Finding common side-effects across all seed drugs")
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
print("\n[4/4] Saving results")
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

print("\n✓ Seed side-effects extraction complete")
print("  Next: Run 03_search_candidate_drugs.py")
