# RESKO FAITHFUL: Quick Start Guide

## ⚡ 60-Second Overview

You're implementing **Ken McGarry's RESKO drug repositioning method** for eEF1A inhibitors.

**Original paper**: McGarry et al. (2018) — applied to Alzheimer's disease
**Your adaptation**: Applied to eEF1A1/eEF1A2 protein targets

**Method**: Find drugs with similar side-effects to known eEF1A inhibitors → score by McGarry's Jaccard formula

**Runtime**: ~30 minutes (mostly data downloads)

---

## 🚀 Quick Start (5 Steps)

### Step 1: Setup
```bash
cd resko-faithful
pip install -r requirements.txt
```

### Step 2: Download Data
```bash
python src/01_download_drugbank_sider.py
```
⚠️ Then manually download DrugBank from https://www.drugbank.ca/releases/latest/downloads/
- all-drug-links.csv.zip
- indications.csv.zip  
- all-targets.csv.zip

Place in: `data/raw/`

### Step 3: Extract Seed Side-Effects
```bash
python src/02_extract_seed_sideeffects.py
```

### Step 4: Search & Score
```bash
python src/03_search_candidate_drugs.py
python src/04_score_candidates_jaccard.py
```

### Step 5: Visualize
```bash
python src/05_visualize_network_3d.py
```

---

## 📊 View Results

**CSV ranking** (open in Excel/Sheets):
```
data/processed/resko_faithful_candidates_scored.csv
```

**Interactive 3D graph** (open in browser):
```
results/resko_faithful_network_3d.html
```

---

## 🧬 Seed Compounds (eEF1A Inhibitors)

These 6 known eEF1A inhibitors serve as your entry point:

1. CHEMBL1802814 — Diphtheria toxin (eEF1A2)
2. CHEMBL1802815 — Translation inhibitor (eEF1A1)
3. CHEMBL5653589 — EEF1G targeting
4. CHEMBL1802973 — eEF1A2 specific
5. CHEMBL1232461 — Molibresib (HSF1/eEF1A1)
6. CHEMBL1221911 — Lactimidomycin

All candidates are ranked by side-effect similarity to these seeds.

---

## 📈 Understanding Your Results

**composite_score > 0.4**: High confidence candidates ✓
- Strong side-effect overlap
- Hit multiple translation machinery targets
- Likely eEF1A modulators

**composite_score 0.2-0.4**: Medium confidence ?
- Partial overlap
- Worthy of literature review

**composite_score < 0.2**: Low confidence ✗
- Marginal overlap
- Likely false positive

---

## 📚 Key Formulas

**McGarry's Three Components:**

```
inds_score[i] = n_indications[i] / total_indications
se_score[i] = shared_side_effects[i] / total_common_side_effects  
on_score[i] = n_targets[i] / total_targets
```

**Composite Score (Jaccard):**
```
Jaccard(X) = sum(X²) / (2·sum(X) - sum(X²))
```

---

## 🔗 Comparison to Your Other Work

| Method | Entry | Runtime | Output |
|--------|-------|---------|--------|
| **RESKO** | eEF1A + seed drugs | 30 min | Side-effect ranking |
| **Hetnet** | eEF1A network | 2-3 hrs | Network-validated ranking |
| **Both** | Combined | 3-4 hrs | High-confidence intersection |

Best: Run both, compare top 20 candidates, focus on overlap.

---

## 📖 Full Documentation

- **README.md** — Project overview & pipeline
- **IMPLEMENTATION_NOTES.md** — Methodology verification (McGarry faithful)
- **RESKO_Faithful_Implementation_Guide.txt** — Comprehensive guide
- **RESKO_Faithful_Deployment_Summary.txt** — Local execution guide

---

## ❓ Common Issues

**"No common side-effects"?**
→ Normal! McGarry's code will prune seeds randomly and retry.

**DrugBank files not found?**
→ Download manually from drugbank.ca (requires registration)

**Very few candidates?**
→ Check that SIDER4 downloaded successfully (~500MB file)

---

## ✅ Status

✓ Faithful implementation of McGarry et al. (2018)  
✓ All 5 pipeline scripts ready  
✓ Git repo initialized  
✓ Ready to push to GitHub Faithful branch  

---

**Next**: Clone repo → download data → run pipeline → analyze results!
