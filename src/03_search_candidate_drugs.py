"""
RESKO Faithful: Step 3 - Search Candidate Drugs

Input: Common seed side-effects + SIDER4 all-drugs data
Output: Candidate drugs with >25% side-effect coverage

McGarry filter: Only keep drugs sharing >25% of common seed side-effects
"""

import pandas as pd
import json
from pathlib import Path

DATA_DIR = Path("data/raw")
PROCESSED_DIR = Path("data/processed")

print("="*70)
print("RESKO FAITHFUL: SEARCH CANDIDATE DRUGS BY SIDE-EFFECT COVERAGE")
print("="*70)

# Load common seed side-effects
print("\n[1/3] Loading seed side-effects")
with open(PROCESSED_DIR / "seed_sideeffects.json") as f:
    seed_data = json.load(f)
common_ses = set(seed_data["common_sideeffects"])
seed_drugs = set(seed_data["seed_drugs"])

print(f"  Common side-effects: {len(common_ses)}")
print(f"  Seed drugs: {len(seed_drugs)}")

# Load SIDER data
print("\n[2/3] Scanning SIDER4 for candidate drugs")
sider_df = pd.read_parquet(DATA_DIR / "sider_all_se.parquet")

# Group by drug, count side-effects
candidates = []
for drug_id, group in sider_df.groupby('drugbank_id'):
    if drug_id in seed_drugs:
        continue  # Skip seed drugs
    
    drug_ses = set(group['side_effect_name'].unique())
    matching_ses = drug_ses & common_ses
    
    # McGarry filter: >25% coverage
    coverage = (len(matching_ses) / len(common_ses)) * 100 if common_ses else 0
    
    if coverage > 25:
        drug_name = group['drug_name'].iloc[0]
        candidates.append({
            'drugbank_id': drug_id,
            'drug_name': drug_name,
            'n_matching_ses': len(matching_ses),
            'coverage_percent': coverage,
            'total_ses': len(drug_ses),
        })

candidates_df = pd.DataFrame(candidates).sort_values('coverage_percent', ascending=False)

print(f"\n  Candidates found with >25% coverage: {len(candidates_df)}")
if len(candidates_df) > 0:
    print(f"\n  Top 10 by coverage:")
    for idx, row in candidates_df.head(10).iterrows():
        print(f"    {row['drug_name']:30s} ({row['drugbank_id']:10s}): {row['coverage_percent']:5.1f}%")

# Save
print("\n[3/3] Saving candidate list")
candidates_df.to_csv(PROCESSED_DIR / "candidate_drugs.csv", index=False)
print(f"  ✓ Saved {len(candidates_df)} candidates to data/processed/candidate_drugs.csv")

print("\n✓ Candidate search complete")
print("  Next: Run 04_score_candidates_jaccard.py")
