"""
RESKO Faithful: Step 1 - Download DrugBank and SIDER4 Data

Downloads and parses:
1. DrugBank 4.0+ (drug indications, targets, adverse reactions)
2. SIDER4 (drug side-effect associations from FDA pharmacovigilance)

McGarry's original RESKO method (2018) uses these exact sources.
"""

import os
import urllib.request
import pandas as pd
import gzip
import shutil
from pathlib import Path

DATA_DIR = Path("data/raw")
DATA_DIR.mkdir(parents=True, exist_ok=True)

print("="*70)
print("RESKO FAITHFUL: DOWNLOADING DRUGBANK AND SIDER4")
print("="*70)

# ===== DRUGBANK =====
print("\n[1/4] DrugBank open-access data")
print("-" * 70)

# DrugBank 4.0+ release files (open-access)
drugbank_urls = {
    "drug_names_and_synonyms": "https://www.drugbank.ca/releases/latest/downloads/all-drug-links.csv.zip",
    "indications": "https://www.drugbank.ca/releases/latest/downloads/indications.csv.zip",
    "targets": "https://www.drugbank.ca/releases/latest/downloads/all-targets.csv.zip",
}

try:
    for name, url in drugbank_urls.items():
        try:
            print(f"  Downloading {name}...", end=" ", flush=True)
            filepath = DATA_DIR / f"drugbank_{name}.csv.zip"
            # Note: This will require DrugBank registration for full data
            # For now, we'll prepare the structure and note the manual step
            print(f"MANUAL STEP REQUIRED: Download from {url}")
        except Exception as e:
            print(f"ERROR: {e}")
except Exception as e:
    print(f"DrugBank download requires manual setup. See instructions.")

# ===== SIDER4 =====
print("\n[2/4] SIDER4 side-effect data (FDA pharmacovigilance)")
print("-" * 70)

sider_urls = {
    "meddra_all_se": "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz",
    "meddra_freq": "http://sideeffects.embl.de/media/download/meddra_freq.tsv.gz",
}

for name, url in sider_urls.items():
    print(f"  Downloading {name}...", end=" ", flush=True)
    try:
        filepath = DATA_DIR / f"sider_{name}.tsv.gz"
        urllib.request.urlretrieve(url, filepath)
        
        # Decompress
        output_path = DATA_DIR / f"sider_{name}.tsv"
        with gzip.open(filepath, 'rb') as f_in:
            with open(output_path, 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        print(f"✓ saved to {output_path}")
    except Exception as e:
        print(f"ERROR: {e}")

print("\n[3/4] Loading SIDER data into memory")
print("-" * 70)

try:
    # SIDER format: drugbank_id, drug_name, umls_id, meddra_id, side_effect_name, meddra_type
    sider_file = DATA_DIR / "sider_meddra_all_se.tsv"
    if sider_file.exists():
        sider_df = pd.read_csv(sider_file, sep='\t', header=None, 
                               names=['drugbank_id', 'drug_name', 'umls_id', 'meddra_id', 'side_effect_name', 'meddra_type'])
        print(f"  Loaded {len(sider_df)} drug-side effect pairs")
        print(f"  Unique drugs: {sider_df['drugbank_id'].nunique()}")
        print(f"  Unique side-effects: {sider_df['side_effect_name'].nunique()}")
        
        # Save summary
        sider_df.to_parquet(DATA_DIR / "sider_all_se.parquet")
    else:
        print(f"  File not found: {sider_file}")
except Exception as e:
    print(f"  ERROR loading SIDER: {e}")

print("\n[4/4] DrugBank Setup Instructions")
print("-" * 70)
print("""
DrugBank is available via registration at https://www.drugbank.ca/

For RESKO Faithful, you need:
1. Download the CSV releases:
   - drug_links.csv (or equivalent drug identifiers)
   - indications.csv (drug → therapeutic indication mapping)
   - all-targets.csv (drug → protein target mapping)
   
2. Place in: data/raw/drugbank_*.csv

The code will parse these in step 02_extract_seed_sideeffects.py
""")

# ===== DRUGBANK ADR FALLBACK =====
# Several seed compounds (didemnin B, metarrestin, ternatin-4, narciclasine,
# nannocystin Ax, BE-43547A2, plitidepsin) are investigational or were
# approved after SIDER4's ~2015 curation date, so SIDER4 has no entry for
# them. 02_extract_seed_sideeffects.py needs a second side-effect source for
# these — DrugBank's own adverse-reaction / clinical-trial-AE fields.
#
# This is manual-download, same as the indications/targets files above:
# DrugBank's structured ADR field is only exposed via the full XML database
# release (requires a DrugBank academic license), not the open CSV files.
print("\n[EXTRA] DrugBank ADR fallback data (for seeds absent from SIDER4)")
print("-" * 70)
print("""
MANUAL STEP REQUIRED — DrugBank ADR / adverse-event data is only in the full
XML release (drugbank_all_full_database.xml.zip), gated behind a DrugBank
academic license (https://go.drugbank.com/releases/latest).

1. Download + unzip the full database release.
2. Parse <drug><adverse-reactions> (or, if unavailable for an investigational
   compound, extract reported adverse events from its clinical-trial results
   / FDA orange-book entry, and record the source in seed_compounds.py).
3. Build a long-format table with columns [drugbank_id, side_effect_name]
   and save it to data/raw/drugbank_adr.parquet — 02_extract_seed_sideeffects.py
   reads this file directly and will proceed with SIDER4-only data (skipping
   the fallback) if it is absent, printing a warning per seed affected.

Until this file is populated, run 02_extract_seed_sideeffects.py anyway —
it will report per-seed coverage (data/processed/seed_coverage.csv) so you
can see exactly which seeds have no side-effect data at all rather than
silently proceeding on a false assumption.
""")

print("\n✓ Data download phase complete")
print("  Next: Run 02_extract_seed_sideeffects.py")
