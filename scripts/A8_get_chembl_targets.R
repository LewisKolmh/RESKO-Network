library(jsonlite)
library(readr)
library(dplyr)

proteins <- read_csv(
  "results/A8_target_proteins.csv",
  show_col_types = FALSE
)

results <- list()

for (p in proteins$protein) {

  cat("Searching:", p, "\n")

  url <- paste0(
    "https://www.ebi.ac.uk/chembl/api/data/target/search.json?q=",
    p
  )

  res <- tryCatch(
    fromJSON(url),
    error = function(e) NULL
  )

  if (
    !is.null(res) &&
    "targets" %in% names(res) &&
    nrow(res$targets) > 0
  ) {

    results[[p]] <- data.frame(
      protein = p,
      chembl_target_id = res$targets$target_chembl_id[1],
      pref_name = res$targets$pref_name[1],
      organism = res$targets$organism[1]
    )

  }

}

target_map <- bind_rows(results)

write_csv(
  target_map,
  "results/A8_chembl_target_map.csv"
)

cat(
  "\nTargets Matched:",
  nrow(target_map),
  "\n"
)

print(target_map)