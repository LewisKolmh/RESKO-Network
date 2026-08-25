"""
RESKO Faithful: Seed compound registry (single source of truth)

Every seed listed here is a CONFIRMED DIRECT eEF1A BINDER — evidence is either
(a) a quantitative ChEMBL binding/inhibition record against an eEF1A-family
target, or (b) a literature-reported direct-binding mechanism (crystal
structure, biochemical trapping of an eEF1A conformational state, etc.).

Side-effect data source is two-tier per McGarry's original SE-similarity step:
  1. SIDER4 (drug package-insert side effects) — used when the compound was a
     marketed drug at SIDER4's curation date (~2015).
  2. DrugBank ADR / clinical-trial adverse-event fields — fallback for
     compounds that are investigational or were approved after 2015, so they
     have no SIDER4 entry. This is a deliberate deviation from McGarry's
     SIDER4-only approach — flag it in any methods write-up.

`sider4_expected=False` marks compounds expected to need the DrugBank ADR
fallback; `01_download_drugbank_sider.py` / `02_extract_seed_sideeffects.py`
verify this empirically rather than assuming it.
"""

SEED_COMPOUNDS = {
    # ---- Tier 1: ChEMBL-verified (quantitative IC50/Kd/ED50 binding records
    #      from the eEF1A-centred hetnet interactome pull) ----
    "CHEMBL1232461": {
        "name": "Molibresib",
        "drugbank_id": None,  # not yet verified against a real DrugBank record
        "evidence": "chembl_binding",
        "notes": "IC50/Kd records against EEF1A1/EEF1G",
        "sider4_expected": False,  # investigational (BET inhibitor, Phase 2)
    },
    "CHEMBL1802814": {
        "name": None,
        "drugbank_id": None,
        "evidence": "chembl_binding",
        "notes": "Quantitative activity record vs eEF1A-family target",
        "sider4_expected": False,
    },
    "CHEMBL1802815": {
        "name": None,
        "drugbank_id": None,
        "evidence": "chembl_binding",
        "notes": "Quantitative activity record vs eEF1A-family target",
        "sider4_expected": False,
    },
    "CHEMBL1802973": {
        "name": None,
        "drugbank_id": None,
        "evidence": "chembl_binding",
        "notes": "Quantitative activity record vs eEF1A-family target",
        "sider4_expected": False,
    },
    "CHEMBL3752910": {
        "name": None,
        "drugbank_id": None,
        "evidence": "chembl_binding",
        "notes": "Quantitative activity record vs eEF1A-family target",
        "sider4_expected": False,
    },
    "CHEMBL5653589": {
        "name": None,
        "drugbank_id": None,
        "evidence": "chembl_binding",
        "notes": "EEF1G-complex targeting, quantitative record",
        "sider4_expected": False,
    },

    # ---- Tier 2: Literature-confirmed direct eEF1A binders ----
    "Plitidepsin": {
        "name": "Plitidepsin (Aplidin)",
        "drugbank_id": "DB04977",
        "evidence": "direct_binding_structural",
        "notes": "Approved (Australia, 2018); binds eEF1A directly; antiviral "
                 "activity (incl. SARS-CoV-2) attributed to eEF1A inhibition. "
                 "Approved AFTER SIDER4 curation -> likely no SIDER4 entry.",
        "sider4_expected": False,
    },
    "Didemnin_B": {
        "name": "Didemnin B",
        "drugbank_id": None,
        "evidence": "direct_binding_structural",
        "notes": "Binds eEF1A between domains I/III, trapping GTP-bound "
                 "conformation; site shared with ternatin-4 and nannocystin A. "
                 "No evidence of direct HIV-1 protein binding -> eEF1A-specific. "
                 "Never marketed -> no SIDER4 entry expected.",
        "sider4_expected": False,
    },
    "Metarrestin": {
        "name": "Metarrestin",
        "drugbank_id": None,
        "evidence": "direct_binding_functional",
        "notes": "eEF1A2-selective; Phase I clinical trial, never marketed.",
        "sider4_expected": False,
    },
    "Ternatin_4": {
        "name": "Ternatin-4",
        "drugbank_id": None,
        "evidence": "direct_binding_structural",
        "notes": "Synthetic analog; occupies same eEF1A site as didemnin B; "
                 "traps aminoacyl-tRNA-accommodation intermediate. Preclinical only.",
        "sider4_expected": False,
    },
    "Narciclasine": {
        "name": "Narciclasine",
        "drugbank_id": None,
        "evidence": "direct_binding_invivo",
        "notes": "Natural product, in vivo validated eEF1A inhibitor. Preclinical only.",
        "sider4_expected": False,
    },
    "Nannocystin_Ax": {
        "name": "Nannocystin Ax",
        "drugbank_id": None,
        "evidence": "direct_binding_invivo",
        "notes": "Shares didemnin-B/ternatin-4 eEF1A binding site. Preclinical only.",
        "sider4_expected": False,
    },
    "BE_43547A2": {
        "name": "BE-43547A2",
        "drugbank_id": None,
        "evidence": "direct_binding_invivo",
        "notes": "In vivo validated direct eEF1A inhibitor. Preclinical only.",
        "sider4_expected": False,
    },
}

# EXCLUDED — kept here with reasons so the exclusion survives code review
EXCLUDED_COMPOUNDS = {
    "Lactimidomycin": "Binds the ribosomal E-site to block the translocation "
                       "step that eEF1A drives; does not directly bind eEF1A "
                       "itself. Fails the direct-binder inclusion criterion.",
    "Efavirenz": "HIV-1 reverse-transcriptase inhibitor; no direct eEF1A "
                 "binding evidence. Side-effect profile reflects RT toxicity, "
                 "not eEF1A activity — would dilute the SE signature.",
    "Cycloheximide": "Binds the ribosomal E-site (same site class as "
                      "lactimidomycin), not eEF1A directly; also too toxic "
                      "for therapeutic use.",
    "Anisomycin": "Binds the ribosomal A-site / peptidyl transferase center, "
                  "not eEF1A.",
    "Aminoglycosides": "Broad-spectrum ribosomal/antibiotic mechanism, not "
                        "eEF1A-specific.",
}

if __name__ == "__main__":
    print(f"Active seeds: {len(SEED_COMPOUNDS)}")
    for k, v in SEED_COMPOUNDS.items():
        print(f"  {k:16s} evidence={v['evidence']:26s} sider4_expected={v['sider4_expected']}")
    print(f"\nExcluded: {len(EXCLUDED_COMPOUNDS)}")
    for k, v in EXCLUDED_COMPOUNDS.items():
        print(f"  {k}: {v}")
