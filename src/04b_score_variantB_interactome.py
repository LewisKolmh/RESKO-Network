"""
RESKO Faithful: Step 4b - Score Variant B (interactome-supplemented)

Rationale (user-directed): McGarry's se_score used side-effect similarity as
a phenotypic PROXY for "these drugs act on similar pathways/targets". With
real side-effect data confirmed to give zero overlap for this eEF1A-binder
seed set (see METHODS_DEVIATIONS.md), we substitute the mechanistic signal
the proxy was standing in for: how strongly a candidate's own targets are
directly, biophysically coupled to eEF1A within the built protein-protein
interactome (472 nodes / 15714 edges, combining STRING + IntAct + Hetionet).

IMPORTANT CORRECTNESS NOTE: this interactome was built as a first-shell-only
network -- ALL 470 non-seed nodes are, by construction, direct (1-hop)
partners of EEF1A1/EEF1A2 (verified: n_seed_links > 0 for 100% of nodes).
There is no second-shell structure to exploit, so a binary/hop-distance
proximity measure is degenerate here (constant across all nodes) and was
abandoned after producing a flat pathway_score in an earlier draft of this
script. Real graph-distance-based differentiation is impossible with this
particular interactome's topology.

Instead, pathway_score uses the *evidence strength* of each node's direct
edge(s) to EEF1A1/EEF1A2 -- a genuinely varying, real quantity already
present in interactome_edges.csv (tier == 'seed_incident' rows):
  coupling(protein) = n_sources                      (1-3 supporting DBs)
                     + string_score  (if present)     (STRING confidence, 0-1)
                     + intact_mi_score (if present)   (IntAct MI score, 0-1)
                     + log1p(n_records)                (evidence-record count)
  pathway_score_raw(candidate) = mean of coupling(t) over the candidate's
                                  interactome-linked targets t (targets_hit)
  pathway_score = min-max normalized pathway_score_raw across candidates

This is a genuine deviation from McGarry (2018) -- it replaces phenotypic
(side-effect) similarity with direct physical-interaction evidence strength
as the third scoring dimension -- disclosed explicitly in
METHODS_DEVIATIONS.md, and it is NOT a reconstruction of DrugBank/SIDER4
se_score.
"""

import pandas as pd
import numpy as np
from pathlib import Path

PROCESSED_DIR = Path("data/processed")
RAW_DIR = Path("data/raw")

print("="*70)
print("RESKO FAITHFUL: SCORE VARIANT B (SE-FREE + INTERACTOME-SUPPLEMENTED)")
print("="*70)

def min_max_norm(s):
    lo, hi = s.min(), s.max()
    if hi == lo:
        return pd.Series(0.0, index=s.index)
    return (s - lo) / (hi - lo)

def jaccard_composite(scores_2d):
    s = np.sum(scores_2d**2, axis=1)
    t = np.sum(scores_2d, axis=1)
    denom = 2*t - s
    out = np.zeros_like(s)
    nz = denom != 0
    out[nz] = s[nz] / denom[nz]
    return out

print("\n[1/4] Loading Variant A candidates and interactome seed-incident edges")
cand = pd.read_csv(PROCESSED_DIR / "resko_variantA_no_se_candidates.csv")
nodes = pd.read_csv(RAW_DIR / "interactome_nodes.csv")
edges = pd.read_csv(RAW_DIR / "interactome_edges.csv")
print(f"  {len(cand)} candidates, {len(nodes)} interactome nodes, {len(edges)} edges")

seed_incident = edges[edges["tier"] == "seed_incident"].copy()
print(f"  {len(seed_incident)} seed-incident edges (direct EEF1A1/EEF1A2 links)")

# Real per-node coupling strength to eEF1A: combine n_sources, string_score,
# intact_mi_score, and evidence-record count. Aggregate over both directions
# (source/target) and, where a protein has edges to both EEF1A1 and EEF1A2,
# take the strongest link (max), not a sum, so proteins linked to both
# paralogs aren't unfairly favored by double-counting.
def edge_coupling(row):
    c = row["n_sources"] if pd.notna(row["n_sources"]) else 0.0
    if pd.notna(row["string_score"]):
        c += row["string_score"]
    if pd.notna(row["intact_mi_score"]):
        c += row["intact_mi_score"]
    if pd.notna(row["n_records"]):
        c += np.log1p(row["n_records"])
    return c

seed_incident["coupling"] = seed_incident.apply(edge_coupling, axis=1)
# the non-EEF1A endpoint of each seed-incident edge is the "partner protein"
seed_incident["partner"] = seed_incident.apply(
    lambda r: r["target"] if r["source"] in ("EEF1A1", "EEF1A2") else r["source"], axis=1)
node_coupling = seed_incident.groupby("partner")["coupling"].max().to_dict()

print(f"  Coupling strength computed for {len(node_coupling)} proteins "
      f"(range {min(node_coupling.values()):.2f} - {max(node_coupling.values()):.2f})")

def pathway_score_for_targets(targets_str):
    if pd.isna(targets_str) or not targets_str:
        return 0.0
    targs = [t for t in str(targets_str).split(";") if t]
    if not targs:
        return 0.0
    weights = [node_coupling.get(t, 0.0) for t in targs]
    return float(np.mean(weights))

print("\n[2/4] Computing interactome proximity weights per candidate")
# Use targets_hit (the actual interactome nodes this candidate is linked to,
# from step 2's aggregation) rather than moa_targets -- targets_hit is
# guaranteed to be a subset of the 472-node interactome by construction,
# so every weight lookup is a real hit, not a name-mismatch miss.
cand["pathway_score_raw"] = cand["targets_hit"].apply(pathway_score_for_targets)
cand["pathway_score"] = min_max_norm(cand["pathway_score_raw"])

print("\n[3/4] Computing Variant B's Jaccard composite with "
      "[inds_score, on_score, pathway_score] -- note Variant A's own "
      "composite_score_variantA column is RETAINED unchanged in this output "
      "file for direct side-by-side comparison between the two variants.")
scores = cand[["inds_score", "on_score", "pathway_score"]].to_numpy()
cand["composite_score_variantB"] = jaccard_composite(scores)
cand["rank_variantA"] = cand["composite_score_variantA"].rank(ascending=False, method="min").astype(int)
cand = cand.sort_values("composite_score_variantB", ascending=False)
cand["rank_variantB"] = range(1, len(cand) + 1)
cand["rank_shift_B_vs_A"] = cand["rank_variantA"] - cand["rank_variantB"]

print(f"\n  Top 15 Variant B candidates by composite score "
      f"(rank_shift_B_vs_A: positive = moved UP under Variant B):")
for idx, row in cand.head(15).iterrows():
    print(f"    {str(row['drug_name'])[:26]:26s} | B_score: {row['composite_score_variantB']:.4f} "
          f"(rank {row['rank_variantB']:>3d}) | A_score: {row['composite_score_variantA']:.4f} "
          f"(rank {row['rank_variantA']:>3d}) | shift: {row['rank_shift_B_vs_A']:+4d} | "
          f"pathway: {row['pathway_score']:.3f}")

print("\n[4/4] Saving")
out_b = PROCESSED_DIR / "resko_variantB_interactome_candidates.csv"
cand.to_csv(out_b, index=False)
print(f"  ✓ Saved Variant B to {out_b}")

print("\n✓ Variant B scoring complete")
print("  Next: Run 05_visualize_network_3d.py")
