library(httr)
library(jsonlite)
library(dplyr)
library(readr)

# STRING Species ID
species <- 9606

# Seed proteins
proteins <- c("EEF1A1", "EEF1A2")

all_results <- list()

for(protein in proteins){

  cat("Querying:", protein, "\n")

  request_url <- paste0(
    "https://string-db.org/api/json/network?",
    "identifiers=", protein,
    "&species=", species
  )

  response <- GET(request_url)

  data <- fromJSON(
    content(response, "text", encoding = "UTF-8")
  )

  if(nrow(data) > 0){

    results <- data.frame(
      query_protein = protein,
      preferredName_A = data$preferredName_A,
      preferredName_B = data$preferredName_B,
      score = data$score
    )

    all_results[[protein]] <- results

  }

}

final_results <- bind_rows(all_results)

if(!dir.exists("results")){
  dir.create("results")
}

write_csv(
  final_results,
  "results/eef1a_string_interactions.csv"
)

print(head(final_results))

cat("\nInteractions saved.\n")