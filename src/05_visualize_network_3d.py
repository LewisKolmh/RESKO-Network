"""
RESKO Faithful: Step 5 - Render comparison network graphs for Variant A and
Variant B.

Both figures share the SAME node layout (computed once from the union of
top-20 candidates across both variants, seeds, and eEF1A) so a reader can
visually track any candidate's position between the two panels. What
differs between the figures is what is real and different between the
variants:
  - Variant A: edge weight = composite_score_variantA (indication-breadth +
    on-target-promiscuity only; se_score fixed at 0/undefined)
  - Variant B: edge weight = composite_score_variantB (adds pathway_score,
    the interactome direct-coupling-strength term); edges from candidates
    whose pathway_score > 0 (i.e., an interactome-supported link exists)
    are drawn as solid, all others as dashed -- this is the one visual
    encoding with no equivalent in Variant A, because pathway_score does
    not exist there.
"""

import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib as mpl
from pathlib import Path

PROCESSED_DIR = Path("data/processed")

# Minimal inline styling (this script runs standalone via subprocess, not the
# interactive kernel, so the figure-style skill's kernel-loaded helpers are
# not importable here -- these are self-contained equivalents).
def apply_figure_style():
    mpl.rcParams.update({
        "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.frameon": False,
        "savefig.dpi": 300, "figure.dpi": 150,
    })

def panel_letter(ax, letter):
    ax.text(-0.02, 1.05, letter, transform=ax.transAxes, fontsize=13,
            fontweight="bold", va="top", ha="right")

print("="*70)
print("RESKO FAITHFUL: RENDER VARIANT A vs VARIANT B NETWORK GRAPHS")
print("="*70)

variantB = pd.read_csv(PROCESSED_DIR / "resko_variantB_interactome_candidates.csv")
N_TOP = 20
topA = variantB.nlargest(N_TOP, "composite_score_variantA")
topB = variantB.nlargest(N_TOP, "composite_score_variantB")
union_names = sorted(set(topA["drug_name"]) | set(topB["drug_name"]))
print(f"  Top {N_TOP} per variant, union = {len(union_names)} unique candidates")

union_df = variantB[variantB["drug_name"].isin(union_names)].set_index("drug_name")

# Build ONE shared layout: eEF1A at center, candidates evenly spaced on a
# ring (guarantees zero label collision for this star topology -- spring
# layout was tried first and produced two overlapping node positions).
n = len(union_names)
angles = np.linspace(0, 2*np.pi, n, endpoint=False)
pos = {"eEF1A1/eEF1A2": np.array([0.0, 0.0])}
for name, theta in zip(union_names, angles):
    pos[name] = np.array([np.cos(theta), np.sin(theta)]) * 1.6

apply_figure_style()
fig, axes = plt.subplots(1, 2, figsize=(15, 7.5))

def draw_variant(ax, score_col, title, show_interactome_support=False):
    G = nx.Graph()
    G.add_node("eEF1A1/eEF1A2", kind="target")
    for name in union_names:
        row = union_df.loc[name]
        G.add_node(name, kind="candidate", score=row[score_col])
        supported = show_interactome_support and row.get("pathway_score", 0) > 0
        G.add_edge("eEF1A1/eEF1A2", name, weight=row[score_col], supported=supported)

    scores = np.array([G.nodes[n].get("score", 0) for n in union_names])
    norm = mpl.colors.Normalize(vmin=0, vmax=max(scores.max(), 1e-6))
    cmap = mpl.cm.viridis

    # target node
    nx.draw_networkx_nodes(G, pos, nodelist=["eEF1A1/eEF1A2"], ax=ax,
                            node_size=1400, node_color="#d62728", edgecolors="black", linewidths=1.2)
    # candidate nodes, colored by this variant's own score
    node_colors = [cmap(norm(G.nodes[n]["score"])) for n in union_names]
    nx.draw_networkx_nodes(G, pos, nodelist=union_names, ax=ax,
                            node_size=420, node_color=node_colors, edgecolors="black", linewidths=0.6)

    if show_interactome_support:
        solid = [(u, v) for u, v, d in G.edges(data=True) if d["supported"]]
        dashed = [(u, v) for u, v, d in G.edges(data=True) if not d["supported"]]
        nx.draw_networkx_edges(G, pos, edgelist=solid, ax=ax, style="solid", width=1.6, edge_color="#444444")
        nx.draw_networkx_edges(G, pos, edgelist=dashed, ax=ax, style="dashed", width=0.8, edge_color="#aaaaaa")
    else:
        nx.draw_networkx_edges(G, pos, ax=ax, style="solid", width=1.0, edge_color="#999999")

    # Radial label placement just outside each node, angle-aligned so text
    # reads outward from the ring without colliding with neighbors.
    for name in union_names:
        x, y = pos[name]
        theta = np.arctan2(y, x)
        lx, ly = x * 1.14, y * 1.14
        rot = np.degrees(theta)
        ha = "left" if -90 < np.degrees(theta) <= 90 else "right"
        if ha == "right":
            rot += 180
        ax.text(lx, ly, name, fontsize=6, rotation=rot, ha=ha, va="center",
                 rotation_mode="anchor")
    nx.draw_networkx_labels(G, pos, labels={"eEF1A1/eEF1A2": "eEF1A1/\neEF1A2"}, ax=ax, font_size=7, font_weight="bold")

    ax.set_title(title, fontsize=10)
    ax.set_xlim(-2.3, 2.3)
    ax.set_ylim(-2.3, 2.3)
    ax.set_aspect("equal")
    ax.axis("off")
    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.04, pad=0.02)
    cbar.set_label(score_col.replace("composite_score_", "composite score, variant "), fontsize=7)
    cbar.ax.tick_params(labelsize=6)

draw_variant(axes[0], "composite_score_variantA",
             "Variant A: indication-breadth + on-target-promiscuity only\n"
             "(se_score undefined -- no real side-effect data available for this seed set)")
draw_variant(axes[1], "composite_score_variantB",
             "Variant B: Variant A + interactome direct-coupling term\n"
             "(solid edge = candidate target has a direct, evidenced eEF1A interactome link)",
             show_interactome_support=True)

panel_letter(axes[0], "a")
panel_letter(axes[1], "b")

fig.suptitle("RESKO-Faithful candidate ranking: real-data variants without side-effect similarity",
             fontsize=11, y=1.02)
fig.tight_layout()

out_png = "resko_variantA_vs_variantB_network.png"
fig.savefig(out_png, dpi=300, bbox_inches="tight")
print(f"\n  ✓ Saved comparison figure to {out_png}")
print("\n✓ Visualization complete")
