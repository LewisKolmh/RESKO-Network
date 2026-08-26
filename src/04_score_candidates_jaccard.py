"""
RESKO Faithful: Step 4 - Score Candidates with McGarry's Jaccard Composite

McGarry's three-component scoring:
1. inds_score - therapeutic breadth (indications normalized)
2. se_score - side-effect similarity coverage (%)
3. on_score - on-target promiscuity (targets normalized)

Combined via Matrix Jaccard Similarity:
  Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))

This matches McGarry's reposition_validate_score.R exactly.
"""

import pandas as pd
import numpy as np
import json
from pathlib import Path

PROCESSED_DIR = Path("data/processed")
RAW_DIR = Path("data/raw")

print("="*70)
print("RESKO FAITHFUL: SCORE CANDIDATES WITH MCGARRY'S JACCARD COMPOSITE")
print("="*70)
print("""
NOTE ON se_score: real pharmacovigilance data (SIDER4 1430 drugs + FAERS,
synonym-pooled, plus a ClinicalTrials.gov structured-AE check) was
exhausted and confirmed the eEF1A-binder seed set has ZERO common
side-effects across its SE-eligible seeds (Molibresib, Plitidepsin have
data; Didemnin_B and Metarrestin have none in any source, for dated/
trial-status reasons, not a retrieval failure). McGarry's se_score
component is therefore fixed at 0 (undefined, not silently dropped) for
every candidate in this run -- see METHODS_DEVIATIONS.md.
""")

def min_max_norm(s):
    lo, hi = s.min(), s.max()
    if hi == lo:
        return pd.Series(0.0, index=s.index)
    return (s - lo) / (hi - lo)

def jaccard_composite(scores_2d):
    """
    McGarry's Jaccard similarity for an N-component vector.
    Jaccard(X) = sum(X^2) / (2*sum(X) - sum(X^2))
    """
    s = np.sum(scores_2d**2, axis=1)
    t = np.sum(scores_2d, axis=1)
    denom = 2*t - s
    out = np.zeros_like(s)
    nz = denom != 0
    out[nz] = s[nz] / denom[nz]
    return out

# ============================================================
# VARIANT A: indication-breadth + on-target-promiscuity only
# ============================================================
print("\n[1/4] Loading SE-free candidates (Variant A input)")
cand = pd.read_csv(PROCESSED_DIR / "candidate_drugs_no_se.csv")
print(f"  Loaded {len(cand)} candidates")

print("\n[2/4] Computing real inds_score (ChEMBL indication breadth)")
n_ind = pd.read_csv(PROCESSED_DIR / "chembl_n_indications.csv").set_index("molecule_chembl_id")["n_indications"]
cand["n_indications"] = cand["drug_id"].map(n_ind).fillna(0)
cand["inds_score"] = min_max_norm(cand["n_indications"])

print("Computing real on_score (ChEMBL on-target promiscuity)")
# moa_n_targets = number of distinct proteins in the drug's ChEMBL-recorded
# mechanism-of-action set -- a direct promiscuity measure, real ChEMBL data,
# not a placeholder.
cand["on_score"] = min_max_norm(cand["moa_n_targets"].fillna(1))

cand["se_score_undefined"] = 0.0  # explicitly undefined for ALL variants -- see note above

scores = cand[["inds_score", "se_score_undefined", "on_score"]].to_numpy()
cand["composite_score_variantA"] = jaccard_composite(scores)
cand = cand.sort_values("composite_score_variantA", ascending=False)

print(f"\n  Top 15 Variant A candidates by composite score:")
for idx, row in cand.head(15).iterrows():
    print(f"    {str(row['drug_name'])[:30]:30s} | score: {row['composite_score_variantA']:.4f} | "
          f"inds: {row['inds_score']:.3f} | on: {row['on_score']:.3f} | n_ind: {row['n_indications']:.0f}")

out_a = PROCESSED_DIR / "resko_variantA_no_se_candidates.csv"
cand.to_csv(out_a, index=False)
print(f"\n  ✓ Saved Variant A to {out_a}")

# ============================================================
# VARIANT B: adds interactome pathway_score (built in next step)
# ============================================================
print("\n[3/4] Variant B (interactome-supplemented) is built in a separate "
      "step -- see 04b_score_variantB_interactome.py")

print("\n[4/4] Done")
print("\n✓ Scoring complete")
print("  Next: Run 04b_score_variantB_interactome.py, then 05_visualize_network_3d.py")
