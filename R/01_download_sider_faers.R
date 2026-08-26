#' RESKO Faithful: Step 1 - Download DrugBank and SIDER4 Data
#'
#' Downloads and parses:
#' 1. DrugBank 4.0+ (drug indications, targets, adverse reactions) -- MANUAL
#'    STEP, see the printed instructions below.
#' 2. SIDER4 (drug side-effect associations from FDA pharmacovigilance)
#' 3. openFDA/FAERS fallback for seeds absent from SIDER4 (post-2015
#'    approvals / investigational compounds)
#'
#' McGarry's original RESKO method (2018) uses DrugBank + SIDER4. This
#' script mirrors src/01_download_drugbank_sider.py exactly, including its
#' disclosed deviation (openFDA/FAERS substituted for the DrugBank ADR
#' field, which is unavailable while DrugBank's academic XML export remains
#' paused).

suppressPackageStartupMessages({
  library(httr2)
  library(readr)
  library(dplyr)
  library(arrow)
})

DATA_DIR <- file.path("data", "raw")
dir.create(DATA_DIR, showWarnings = FALSE, recursive = TRUE)

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: DOWNLOADING DRUGBANK AND SIDER4\n")
cat(strrep("=", 70), "\n")

# ===== DRUGBANK =====
cat("\n[1/4] DrugBank open-access data\n")
cat(strrep("-", 70), "\n")
drugbank_urls <- c(
  drug_names_and_synonyms = "https://www.drugbank.ca/releases/latest/downloads/all-drug-links.csv.zip",
  indications = "https://www.drugbank.ca/releases/latest/downloads/indications.csv.zip",
  targets = "https://www.drugbank.ca/releases/latest/downloads/all-targets.csv.zip"
)
for (name in names(drugbank_urls)) {
  cat(sprintf("  MANUAL STEP REQUIRED: Download %s from %s\n", name, drugbank_urls[[name]]))
}

# ===== SIDER4 =====
cat("\n[2/4] SIDER4 side-effect data (FDA pharmacovigilance)\n")
cat(strrep("-", 70), "\n")

sider_urls <- c(
  meddra_all_se = "http://sideeffects.embl.de/media/download/meddra_all_se.tsv.gz",
  meddra_freq = "http://sideeffects.embl.de/media/download/meddra_freq.tsv.gz"
)

for (name in names(sider_urls)) {
  gz_path <- file.path(DATA_DIR, sprintf("sider_%s.tsv.gz", name))
  tsv_path <- file.path(DATA_DIR, sprintf("sider_%s.tsv", name))
  if (file.exists(tsv_path)) {
    cat(sprintf("  %s already present at %s, skipping download\n", name, tsv_path))
    next
  }
  cat(sprintf("  Downloading %s... ", name))
  result <- tryCatch({
    req <- request(sider_urls[[name]]) |> req_timeout(60)
    resp <- req_perform(req, path = gz_path)
    con_in <- gzfile(gz_path, "rb")
    writeBin(readBin(con_in, "raw", n = file.info(gz_path)$size * 10), tsv_path)
    close(con_in)
    cat(sprintf("saved to %s\n", tsv_path))
    TRUE
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", conditionMessage(e)))
    FALSE
  })
}

cat("\n[3/4] Loading SIDER data into memory\n")
cat(strrep("-", 70), "\n")

sider_file <- file.path(DATA_DIR, "sider_meddra_all_se.tsv")
if (file.exists(sider_file)) {
  sider_df <- read_tsv(sider_file, col_names = c(
    "drugbank_id", "drug_name", "umls_id", "meddra_id",
    "side_effect_name", "meddra_type"
  ), show_col_types = FALSE)
  cat(sprintf("  Loaded %d drug-side effect pairs\n", nrow(sider_df)))
  cat(sprintf("  Unique drugs: %d\n", n_distinct(sider_df$drugbank_id)))
  cat(sprintf("  Unique side-effects: %d\n", n_distinct(sider_df$side_effect_name)))
  write_parquet(sider_df, file.path(DATA_DIR, "sider_all_se.parquet"))
} else {
  cat(sprintf("  File not found: %s\n", sider_file))
}

cat("\n[4/4] DrugBank Setup Instructions\n")
cat(strrep("-", 70), "\n")
cat(paste(
  "",
  "DrugBank is available via registration at https://www.drugbank.ca/",
  "",
  "For RESKO Faithful, you need:",
  "1. Download the CSV releases:",
  "   - drug_links.csv (or equivalent drug identifiers)",
  "   - indications.csv (drug -> therapeutic indication mapping)",
  "   - all-targets.csv (drug -> protein target mapping)",
  "",
  "2. Place in: data/raw/drugbank_*.csv",
  "",
  "The code will parse these in step 02_extract_seed_sideeffects.R",
  "", sep = "\n"
))

# ===== openFDA/FAERS FALLBACK =====
# Several seed compounds (didemnin B, metarrestin, plitidepsin, and the
# 6 ChEMBL research compounds) are investigational or were approved after
# SIDER4's ~2015 curation date, so SIDER4 has no entry for them.
# 02_extract_seed_sideeffects.R needs a second side-effect source for these.
#
# DrugBank's own ADR field would be the McGarry-faithful fallback, but
# DrugBank's academic full-XML downloads (the only export carrying ADR data)
# are platform-wide paused as of this writing. openFDA/FAERS is used
# instead: free, public, no license required, drawn from real-world
# adverse-event reports. This is a deliberate deviation from McGarry's
# SIDER4-only method -- and a further deviation from the DrugBank-ADR
# fallback originally planned -- disclose both in any methods write-up.
# FAERS is voluntary-report data (reporting bias skews toward
# serious/unexpected events, not true incidence), a different bias profile
# than SIDER4's package-insert extraction.
#
# NOTE: 4 seeds (Ternatin-4, Narciclasine, Nannocystin Ax, BE-43547A2) have
# never been dosed in a human and so have zero FAERS records as a
# structural fact, not a query failure -- see
# seed_compounds.R::binding_evidence_only_seeds().
cat("\n[EXTRA] openFDA/FAERS fallback data (for seeds absent from SIDER4)\n")
cat(strrep("-", 70), "\n")

source(file.path("R", "seed_compounds.R"))

FAERS_ENDPOINT <- "https://api.fda.gov/drug/event.json"

#' Query openFDA FAERS for adverse-event reaction terms mentioning
#' `drug_name` as a suspect medicinal product. Returns a character vector of
#' MedDRA PT (preferred term) reaction strings, or character(0) on
#' no-match/error.
query_faers <- function(drug_name) {
  q <- sprintf('patient.drug.medicinalproduct:"%s"', drug_name)
  tryCatch({
    resp <- request(FAERS_ENDPOINT) |>
      req_url_query(search = q, count = "patient.reaction.reactionmeddrapt.exact") |>
      req_timeout(30) |>
      req_error(is_error = function(resp) FALSE) |>
      req_perform()
    if (resp_status(resp) == 404) return(character(0))  # zero-hit queries
    if (resp_status(resp) >= 400) {
      cat(sprintf("    HTTP error for %s: %d\n", drug_name, resp_status(resp)))
      return(character(0))
    }
    data <- resp_body_json(resp)
    vapply(data$results, function(row) row$term, character(1))
  }, error = function(e) {
    cat(sprintf("    ERROR querying FAERS for %s: %s\n", drug_name, conditionMessage(e)))
    character(0)
  })
}

eligible <- se_eligible_seeds()
# Only query FAERS for seeds not expected to already have SIDER4 coverage,
# and only those with a queryable drug name.
query_targets <- Filter(function(meta) {
  !isTRUE(meta$sider4_expected) && !is.na(meta$name)
}, eligible)

cat(sprintf("  Querying FAERS for %d seed(s) expected to lack SIDER4 coverage...\n",
            length(query_targets)))

faers_rows <- list()
for (cid in names(query_targets)) {
  meta <- query_targets[[cid]]
  name <- strsplit(meta$name, " \\(")[[1]][1]  # strip parenthetical synonyms
  cat(sprintf("  %s (%s)... ", cid, name))
  terms <- query_faers(name)
  cat(sprintf("%d reaction terms\n", length(terms)))
  if (length(terms) > 0) {
    faers_rows[[cid]] <- tibble::tibble(
      seed_id = cid, drugbank_id = meta$drugbank_id, side_effect_name = terms
    )
  }
  Sys.sleep(0.5)  # be polite to the shared public API
}

faers_df <- if (length(faers_rows) > 0) {
  dplyr::bind_rows(faers_rows)
} else {
  tibble::tibble(seed_id = character(0), drugbank_id = character(0),
                 side_effect_name = character(0))
}
write_parquet(faers_df, file.path(DATA_DIR, "faers_adr.parquet"))
cat(sprintf("  Saved %d FAERS adverse-event terms to data/raw/faers_adr.parquet (%d seeds covered)\n",
            nrow(faers_df), n_distinct(faers_df$seed_id)))

cat(paste(
  "",
  "NOTE on the 5 unnamed ChEMBL Tier-1 seeds (name=NA in seed_compounds.R):",
  "these were checked directly against the live ChEMBL API and confirmed to",
  "have no pref_name, no synonyms, and no max_phase -- i.e. they are pure",
  "research/binding-assay compounds that never entered clinical development,",
  "so they were correctly excluded from this FAERS run via",
  "has_human_exposure = FALSE (not merely skipped for lacking a name). FAERS",
  "is keyed on medicinal-product free text, and no such text will ever exist",
  "for these compounds.", "", sep = "\n"
))

cat("\nData download phase complete\n")
cat("  Next: Run 02_extract_seed_sideeffects.R\n")
