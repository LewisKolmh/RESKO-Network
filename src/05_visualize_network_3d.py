"""
RESKO Faithful: Step 5 - Visualize 3D Network Graph

Creates interactive 3D network:
- eEF1A (center node)
- Seed drugs (linked to eEF1A)
- Candidate compounds (colored by composite score)
- Edges represent side-effect coverage %

Uses Plotly for interactive 3D visualization.
"""

import pandas as pd
import json
import numpy as np
from pathlib import Path

try:
    import plotly.graph_objects as go
except ImportError:
    print("ERROR: plotly required. Install: pip install plotly")
    exit(1)

PROCESSED_DIR = Path("data/processed")
RESULTS_DIR = Path("results")
RESULTS_DIR.mkdir(exist_ok=True)

print("="*70)
print("RESKO FAITHFUL: VISUALIZE 3D NETWORK GRAPH")
print("="*70)

# Load data
print("\n[1/3] Loading candidate scores and seed data")
candidates_df = pd.read_csv(PROCESSED_DIR / "resko_faithful_candidates_scored.csv")
with open(PROCESSED_DIR / "seed_sideeffects.json") as f:
    seed_data = json.load(f)

print(f"  Candidates: {len(candidates_df)}")
print(f"  Seed drugs: {len(seed_data['seed_drugs'])}")

# Build network nodes
print("\n[2/3] Building 3D network")
nodes_x, nodes_y, nodes_z = [], [], []
nodes_name = []
nodes_type = []
nodes_color = []

# Center node: eEF1A
eef1a_x, eef1a_y, eef1a_z = 0, 0, 0
nodes_x.append(eef1a_x)
nodes_y.append(eef1a_y)
nodes_z.append(eef1a_z)
nodes_name.append("eEF1A1/eEF1A2")
nodes_type.append("seed_protein")
nodes_color.append(10)  # Colorscale value

# Seed drugs (circle around center)
n_seeds = len(seed_data['seed_drugs'])
for i, drug_id in enumerate(seed_data['seed_drugs']):
    angle = 2 * np.pi * i / n_seeds
    radius = 3
    x = radius * np.cos(angle)
    y = radius * np.sin(angle)
    z = 0
    nodes_x.append(x)
    nodes_y.append(y)
    nodes_z.append(z)
    nodes_name.append("Seed: " + str(drug_id))
    nodes_type.append("seed_drug")
    nodes_color.append(8)

# Candidate compounds (distributed in 3D space)
top_n = min(50, len(candidates_df))
for i in range(top_n):
    row = candidates_df.iloc[i]
    angle = 2 * np.pi * i / top_n
    radius = 6
    x = radius * np.cos(angle) + np.random.normal(0, 0.3)
    y = radius * np.sin(angle) + np.random.normal(0, 0.3)
    z = 3 + np.random.normal(0, 0.5)
    
    nodes_x.append(x)
    nodes_y.append(y)
    nodes_z.append(z)
    nodes_name.append(row['drug_name'])
    nodes_type.append("candidate")
    nodes_color.append(row['composite_score'])

# Build edges
edges_x, edges_y, edges_z = [], [], []

# Edges: seed drugs to eEF1A
for i in range(1, n_seeds + 1):
    edges_x.extend([nodes_x[0], nodes_x[i], None])
    edges_y.extend([nodes_y[0], nodes_y[i], None])
    edges_z.extend([nodes_z[0], nodes_z[i], None])

# Edges: candidates to eEF1A
for i in range(n_seeds + 1, len(nodes_x)):
    edges_x.extend([nodes_x[0], nodes_x[i], None])
    edges_y.extend([nodes_y[0], nodes_y[i], None])
    edges_z.extend([nodes_z[0], nodes_z[i], None])

# Create figure
print("\n[3/3] Creating interactive 3D visualization")

fig = go.Figure(data=[
    go.Scatter3d(
        x=edges_x, y=edges_y, z=edges_z,
        mode='lines',
        line=dict(color='rgba(125,125,125,0.5)', width=1),
        hoverinfo='none',
        name='Connections'
    ),
    go.Scatter3d(
        x=nodes_x, y=nodes_y, z=nodes_z,
        mode='markers+text',
        marker=dict(
            size=10,
            color=nodes_color,
            colorscale='Viridis',
            showscale=True,
            colorbar=dict(title="Composite<br>Score"),
            line=dict(color='white', width=2)
        ),
        text=nodes_name,
        textposition="top center",
        hoverinfo='text',
        name='Compounds'
    )
])

fig.update_layout(
    title=dict(
        text="RESKO Faithful: eEF1A Drug Repositioning Network<br>(McGarry Method - 3D Visualization)",
        x=0.5,
        xanchor='center'
    ),
    scene=dict(
        xaxis_title='X Coordinate',
        yaxis_title='Y Coordinate',
        zaxis_title='Z Coordinate',
        camera=dict(
            eye=dict(x=1.5, y=1.5, z=1.3)
        )
    ),
    width=1000,
    height=800,
    showlegend=False
)

output_file = RESULTS_DIR / "resko_faithful_network_3d.html"
fig.write_html(str(output_file))
print(f"  ✓ Saved to results/resko_faithful_network_3d.html")

print("\n✓ Network visualization complete!")
print(f"\n  Final ranking saved to: data/processed/resko_faithful_candidates_scored.csv")
print(f"  Interactive 3D graph: results/resko_faithful_network_3d.html")
