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

# ===== openFDA/FAERS FALLBACK =====
# Several seed compounds (didemnin B, metarrestin, plitidepsin, and the
# 6 ChEMBL research compounds) are investigational or were approved after
# SIDER4's ~2015 curation date, so SIDER4 has no entry for them.
# 02_extract_seed_sideeffects.py needs a second side-effect source for these.
#
# DrugBank's own ADR field would be the McGarry-faithful fallback, but
# DrugBank's academic full-XML downloads (the only export carrying ADR data)
# are platform-wide paused as of this writing. openFDA/FAERS is used instead:
# free, public, no license required, drawn from real-world adverse-event
# reports. This is a deliberate deviation from McGarry's SIDER4-only method
# -- and a further deviation from the DrugBank-ADR fallback originally
# planned -- disclose both in any methods write-up. FAERS is voluntary-report
# data (reporting bias skews toward serious/unexpected events, not true
# incidence), a different bias profile than SIDER4's package-insert
# extraction.
#
# NOTE: 4 seeds (Ternatin-4, Narciclasine, Nannocystin Ax, BE-43547A2) have
# never been dosed in a human and so have zero FAERS records as a structural
# fact, not a query failure -- see seed_compounds.binding_evidence_only_seeds().
print("\n[EXTRA] openFDA/FAERS fallback data (for seeds absent from SIDER4)")
print("-" * 70)

import json
import sys
import time
import urllib.error

sys.path.insert(0, str(Path(__file__).parent))
from seed_compounds import SEED_COMPOUNDS, se_eligible_seeds

FAERS_ENDPOINT = "https://api.fda.gov/drug/event.json"

def query_faers(drug_name, limit=1000):
    """Query openFDA FAERS for adverse-event reaction terms mentioning
    `drug_name` as a suspect medicinal product. Returns a list of MedDRA
    PT (preferred term) reaction strings, or [] on no-match / error."""
    q = f'patient.drug.medicinalproduct:"{drug_name}"'
    url = (f"{FAERS_ENDPOINT}?search={urllib.parse.quote(q)}"
           f"&count=patient.reaction.reactionmeddrapt.exact")
    try:
        with urllib.request.urlopen(url, timeout=30) as resp:
            data = json.load(resp)
        return [row["term"] for row in data.get("results", [])]
    except urllib.error.HTTPError as e:
        if e.code == 404:
            return []  # openFDA returns 404 for zero-hit queries
        print(f"    HTTP error for {drug_name}: {e.code}")
        return []
    except Exception as e:
        print(f"    ERROR querying FAERS for {drug_name}: {e}")
        return []

import urllib.parse

faers_rows = []
eligible = se_eligible_seeds()
# Only query FAERS for seeds not expected to already have SIDER4 coverage,
# and only those with a queryable drug name.
query_targets = {cid: meta for cid, meta in eligible.items()
                  if not meta.get("sider4_expected", True) and meta.get("name")}

print(f"  Querying FAERS for {len(query_targets)} seed(s) expected to lack "
      f"SIDER4 coverage...")
for cid, meta in query_targets.items():
    name = meta["name"].split(" (")[0]  # strip parenthetical synonyms
    print(f"  {cid} ({name})...", end=" ", flush=True)
    terms = query_faers(name)
    print(f"{len(terms)} reaction terms")
    for t in terms:
        faers_rows.append({"seed_id": cid, "drugbank_id": meta["drugbank_id"],
                            "side_effect_name": t})
    time.sleep(0.5)  # be polite to the shared public API

faers_df = pd.DataFrame(faers_rows, columns=["seed_id", "drugbank_id", "side_effect_name"])
faers_df.to_parquet(DATA_DIR / "faers_adr.parquet")
print(f"  ✓ Saved {len(faers_df)} FAERS adverse-event terms to "
      f"data/raw/faers_adr.parquet ({faers_df['seed_id'].nunique() if len(faers_df) else 0} seeds covered)")

print("""
NOTE on the 5 unnamed ChEMBL Tier-1 seeds (name=None in seed_compounds.py):
these were checked directly against the live ChEMBL API and confirmed to
have no pref_name, no synonyms, and no max_phase -- i.e. they are pure
research/binding-assay compounds that never entered clinical development,
so they were correctly excluded from this FAERS run via
has_human_exposure=False (not merely skipped for lacking a name). FAERS is
keyed on medicinal-product free text, and no such text will ever exist for
these compounds.
""")

print("\n✓ Data download phase complete")
print("  Next: Run 02_extract_seed_sideeffects.py")
