library(readr)
library(dplyr)
library(stringr)

# ==========================================
# A2 - Check whether eEF1A seed drugs exist
# in the SIDER database
# ==========================================

# Load seed drugs
seed_drugs <- read_csv(
  "eef1a_seed_drugs.csv",
  show_col_types = FALSE
)

# Load SIDER drug names
drug_names <- read_tsv(
  "SIDER/drug_names.tsv",
  col_names = c("CID", "DrugName"),
  show_col_types = FALSE
)

# Quick sanity check
cat("\nSeed drugs loaded:\n")
print(seed_drugs)

cat("\nSIDER drug count:\n")
print(nrow(drug_names))

# Search SIDER for each seed drug
results <- lapply(seed_drugs$Drug, function(drug) {

  matches <- drug_names %>%
    filter(
      str_detect(
        str_to_lower(DrugName),
        str_to_lower(drug)
      )
    )

  data.frame(
    SeedDrug = drug,
    PresentInSIDER = nrow(matches) > 0,
    NumberMatches = nrow(matches),
    SIDERMatch = ifelse(
      nrow(matches) > 0,
      paste(matches$DrugName, collapse = "; "),
      ""
    ),
    CID = ifelse(
      nrow(matches) > 0,
      paste(matches$CID, collapse = "; "),
      ""
    )
  )

})

# Combine results
results_df <- bind_rows(results)

# Create results folder if missing
if(!dir.exists("results")){
  dir.create("results")
}

# Save results
write_csv(
  results_df,
  "results/results_A2_seed_lookup.csv"
)

# Print results
cat("\nLookup Results:\n")
print(results_df)

cat("\nFinished.\n")