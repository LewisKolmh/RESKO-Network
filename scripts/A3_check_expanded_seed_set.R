library(readr)
library(dplyr)
library(stringr)

seed_drugs <- read_csv(
  "eef1a_master_seed_list.csv",
  show_col_types = FALSE
)

drug_names <- read_tsv(
  "SIDER/drug_names.tsv",
  col_names = c("CID", "DrugName"),
  show_col_types = FALSE
)

results <- lapply(seed_drugs$Drug, function(drug) {

  matches <- drug_names %>%
    filter(
      str_detect(
        str_to_lower(DrugName),
        str_to_lower(drug)
      )
    )

  data.frame(
    Drug = drug,
    PresentInSIDER = nrow(matches) > 0,
    NumberMatches = nrow(matches),
    Match = ifelse(
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

results_df <- bind_rows(results)

if (!dir.exists("results")) {
  dir.create("results")
}

write_csv(
  results_df,
  "results/A3_expanded_seed_lookup.csv"
)

print(results_df)