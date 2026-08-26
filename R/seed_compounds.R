#' RESKO Faithful: Seed compound registry (single source of truth)
#'
#' Every seed listed here is a CONFIRMED DIRECT eEF1A BINDER -- evidence is
#' either (a) a quantitative ChEMBL binding/inhibition record against an
#' eEF1A-family target, or (b) a literature-reported direct-binding mechanism
#' (crystal structure, biochemical trapping of an eEF1A conformational
#' state, etc.).
#'
#' Side-effect data source is two-tier per McGarry's original SE-similarity
#' step:
#'   1. SIDER4 (drug package-insert side effects) -- used when the compound
#'      was a marketed drug at SIDER4's curation date (~2015).
#'   2. openFDA/FAERS (FDA Adverse Event Reporting System, via the free,
#'      public openFDA API) -- fallback for compounds investigational or
#'      approved after 2015, so absent from SIDER4. NOTE: this replaces an
#'      earlier plan to use DrugBank's own ADR field as the fallback --
#'      DrugBank's academic full-XML downloads (the only export carrying ADR
#'      data) are platform-wide paused as of this writing, so openFDA/FAERS
#'      is used instead. FAERS is voluntary-report data (reporting bias
#'      skews toward serious/unexpected events, not incidence), a different
#'      bias profile than SIDER4's package-insert extraction -- disclose
#'      both the substitution and the bias difference in any methods
#'      write-up.
#'
#' `sider4_expected = FALSE` marks compounds expected to need the FAERS
#' fallback; `01_download_sider_faers.R` / `02_extract_seed_sideeffects.R`
#' verify this empirically rather than assuming it.
#'
#' `has_human_exposure = FALSE` marks compounds that have NEVER been dosed
#' in a human (preclinical-only) -- no side-effect source (SIDER4, FAERS, or
#' DrugBank ADR) can ever have data for these, as a structural fact about
#' the compound, not a data-access gap. These seeds are retained in the
#' registry as binding-evidence-only seeds: they contribute to the seed
#' set's mechanistic/structural justification and to network diagrams, but
#' `02_extract_seed_sideeffects.R` excludes them from the SIDER4/FAERS
#' side-effect intersection step entirely (rather than letting a real
#' zero-SE row silently corrupt the intersection).
#'
#' IMPORTANT: only 4 of the 13 registered seeds are SE-eligible --
#' Molibresib, Plitidepsin, Didemnin_B, and Metarrestin. The 5
#' non-Molibresib ChEMBL Tier-1 seeds were checked directly against the live
#' ChEMBL API (molecule endpoint): each returns pref_name = NULL,
#' max_phase = NULL, and no synonyms, confirming they are pure
#' research/binding-assay compounds that never entered clinical
#' development. They join the 4 literature Tier-2 compounds already known
#' to be preclinical-only, so 9 of 13 seeds are binding-evidence-only. This
#' is a real reduction in the statistical power of McGarry's
#' side-effect-intersection step relative to his original (SIDER4-only,
#' marketed-drug) seed sets, and should be disclosed as such in any methods
#' write-up.

suppressPackageStartupMessages(library(tibble))

# ---- Tier 1: ChEMBL-verified (quantitative IC50/Kd/ED50 binding records
#      from the eEF1A-centred hetnet interactome pull) ----
SEED_COMPOUNDS <- list(
  CHEMBL1232461 = list(
    name = "Molibresib",
    drugbank_id = NA_character_,  # not yet verified against a real DrugBank record
    evidence = "chembl_binding",
    notes = "IC50/Kd records against EEF1A1/EEF1G",
    sider4_expected = FALSE,      # investigational (BET inhibitor, Phase 2)
    has_human_exposure = TRUE     # reached Phase 2 clinical trial
  ),
  CHEMBL1802814 = list(
    name = NA_character_,
    drugbank_id = NA_character_,
    evidence = "chembl_binding",
    notes = paste("Quantitative activity record vs eEF1A-family target. ChEMBL API",
                  "confirms no pref_name, no synonyms, no max_phase -- a pure",
                  "research/binding-assay compound, never in clinical",
                  "development. Structurally cannot appear in SIDER4 or FAERS",
                  "(both require a marketed/trialed drug name/identity) and",
                  "almost certainly has no DrugBank ID. Binding-evidence-only."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE    # confirmed via ChEMBL: no clinical phase, no name
  ),
  CHEMBL1802815 = list(
    name = NA_character_,
    drugbank_id = NA_character_,
    evidence = "chembl_binding",
    notes = paste("Quantitative activity record vs eEF1A-family target. Same",
                  "ChEMBL-confirmed profile as CHEMBL1802814: no pref_name, no",
                  "synonyms, no max_phase -- research-only, binding-evidence-only."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),
  CHEMBL1802973 = list(
    name = NA_character_,
    drugbank_id = NA_character_,
    evidence = "chembl_binding",
    notes = paste("Quantitative activity record vs eEF1A-family target. Same",
                  "ChEMBL-confirmed profile: no pref_name, no synonyms, no",
                  "max_phase -- research-only, binding-evidence-only."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),
  CHEMBL3752910 = list(
    name = NA_character_,
    drugbank_id = NA_character_,
    evidence = "chembl_binding",
    notes = paste("Quantitative activity record vs eEF1A-family target. Same",
                  "ChEMBL-confirmed profile: no pref_name, no synonyms, no",
                  "max_phase -- research-only, binding-evidence-only."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),
  CHEMBL5653589 = list(
    name = NA_character_,
    drugbank_id = NA_character_,
    evidence = "chembl_binding",
    notes = paste("EEF1G-complex targeting, quantitative record. Same",
                  "ChEMBL-confirmed profile: no pref_name, no synonyms, no",
                  "max_phase -- research-only, binding-evidence-only."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),

  # ---- Tier 2: Literature-confirmed direct eEF1A binders ----
  Plitidepsin = list(
    name = "Plitidepsin (Aplidin)",
    drugbank_id = "DB04977",
    evidence = "direct_binding_structural",
    notes = paste("Approved (Australia, 2018); binds eEF1A directly; antiviral",
                  "activity (incl. SARS-CoV-2) attributed to eEF1A inhibition.",
                  "Approved AFTER SIDER4 curation -> likely no SIDER4 entry;",
                  "try openFDA/FAERS fallback."),
    sider4_expected = FALSE,
    has_human_exposure = TRUE     # FDA-approved (Australia)
  ),
  Didemnin_B = list(
    name = "Didemnin B",
    drugbank_id = NA_character_,
    evidence = "direct_binding_structural",
    notes = paste("Binds eEF1A between domains I/III, trapping GTP-bound",
                  "conformation; site shared with ternatin-4 and nannocystin A.",
                  "No evidence of direct HIV-1 protein binding -> eEF1A-specific.",
                  "Never marketed but reached clinical trials -> try openFDA/FAERS",
                  "fallback; SIDER4 unlikely."),
    sider4_expected = FALSE,
    has_human_exposure = TRUE     # clinical trials (discontinued for toxicity)
  ),
  Metarrestin = list(
    name = "Metarrestin",
    drugbank_id = NA_character_,
    evidence = "direct_binding_functional",
    notes = paste("eEF1A2-selective; Phase I clinical trial, never marketed.",
                  "Try openFDA/FAERS fallback; SIDER4 unlikely."),
    sider4_expected = FALSE,
    has_human_exposure = TRUE     # Phase I clinical trial
  ),
  Ternatin_4 = list(
    name = "Ternatin-4",
    drugbank_id = NA_character_,
    evidence = "direct_binding_structural",
    notes = paste("Synthetic analog; occupies same eEF1A site as didemnin B;",
                  "traps aminoacyl-tRNA-accommodation intermediate. Preclinical",
                  "only -- NEVER dosed in humans, so no side-effect source",
                  "(SIDER4, FAERS, or DrugBank ADR) can ever have data for it.",
                  "Retained as a binding-evidence-only seed: contributes to the",
                  "seed set's structural justification and network diagrams,",
                  "excluded from the SE-intersection step."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE    # preclinical only -- structural, not a data gap
  ),
  Narciclasine = list(
    name = "Narciclasine",
    drugbank_id = NA_character_,
    evidence = "direct_binding_invivo",
    notes = paste("Natural product, in vivo validated eEF1A inhibitor. Preclinical",
                  "only -- never dosed in humans; binding-evidence-only seed,",
                  "excluded from the SE-intersection step."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),
  Nannocystin_Ax = list(
    name = "Nannocystin Ax",
    drugbank_id = NA_character_,
    evidence = "direct_binding_invivo",
    notes = paste("Shares didemnin-B/ternatin-4 eEF1A binding site. Preclinical",
                  "only -- never dosed in humans; binding-evidence-only seed,",
                  "excluded from the SE-intersection step."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  ),
  BE_43547A2 = list(
    name = "BE-43547A2",
    drugbank_id = NA_character_,
    evidence = "direct_binding_invivo",
    notes = paste("In vivo validated direct eEF1A inhibitor. Preclinical only --",
                  "never dosed in humans; binding-evidence-only seed, excluded",
                  "from the SE-intersection step."),
    sider4_expected = FALSE,
    has_human_exposure = FALSE
  )
)

# EXCLUDED -- kept here with reasons so the exclusion survives code review
EXCLUDED_COMPOUNDS <- list(
  Lactimidomycin = paste("Binds the ribosomal E-site to block the translocation",
                          "step that eEF1A drives; does not directly bind eEF1A",
                          "itself. Fails the direct-binder inclusion criterion."),
  Efavirenz = paste("HIV-1 reverse-transcriptase inhibitor; no direct eEF1A",
                     "binding evidence. Side-effect profile reflects RT toxicity,",
                     "not eEF1A activity -- would dilute the SE signature."),
  Cycloheximide = paste("Binds the ribosomal E-site (same site class as",
                         "lactimidomycin), not eEF1A directly; also too toxic",
                         "for therapeutic use."),
  Anisomycin = paste("Binds the ribosomal A-site / peptidyl transferase center,",
                      "not eEF1A."),
  Aminoglycosides = paste("Broad-spectrum ribosomal/antibiotic mechanism, not",
                           "eEF1A-specific.")
)

#' Seeds that can participate in the SIDER4/FAERS side-effect intersection
#' step -- i.e. NOT flagged has_human_exposure = FALSE. Compounds with
#' has_human_exposure = NA (clinical status not yet established) are still
#' attempted; only a confirmed preclinical-only compound is excluded
#' outright.
se_eligible_seeds <- function() {
  Filter(function(meta) !identical(meta$has_human_exposure, FALSE), SEED_COMPOUNDS)
}

#' Seeds retained for structural/mechanistic justification and network
#' diagrams, but excluded from the SE-intersection step because they have
#' never been dosed in a human (no side-effect source can ever have data
#' for them).
binding_evidence_only_seeds <- function() {
  Filter(function(meta) identical(meta$has_human_exposure, FALSE), SEED_COMPOUNDS)
}

#' Flatten the registry to a tibble for convenient joining/filtering
#' downstream (mirrors the shape used throughout the R pipeline).
seed_compounds_tbl <- function() {
  rows <- lapply(names(SEED_COMPOUNDS), function(id) {
    meta <- SEED_COMPOUNDS[[id]]
    tibble(
      seed_id = id,
      name = meta$name %||% NA_character_,
      drugbank_id = meta$drugbank_id %||% NA_character_,
      evidence = meta$evidence,
      notes = meta$notes,
      sider4_expected = meta$sider4_expected,
      has_human_exposure = meta$has_human_exposure
    )
  })
  do.call(rbind, rows)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

if (sys.nframe() == 0) {
  cat(sprintf("Active seeds: %d\n", length(SEED_COMPOUNDS)))
  for (id in names(SEED_COMPOUNDS)) {
    v <- SEED_COMPOUNDS[[id]]
    cat(sprintf("  %-16s evidence=%-26s sider4_expected=%s has_human_exposure=%s\n",
                id, v$evidence, v$sider4_expected, v$has_human_exposure))
  }
  se <- se_eligible_seeds()
  beo <- binding_evidence_only_seeds()
  cat(sprintf("\nSE-eligible (participate in SIDER4/FAERS intersection): %d\n", length(se)))
  cat(sprintf("Binding-evidence-only (never dosed in humans, no SE data possible): %d -> %s\n",
              length(beo), paste(names(beo), collapse = ", ")))
  cat(sprintf("\nExcluded: %d\n", length(EXCLUDED_COMPOUNDS)))
  for (id in names(EXCLUDED_COMPOUNDS)) {
    cat(sprintf("  %s: %s\n", id, EXCLUDED_COMPOUNDS[[id]]))
  }
}
