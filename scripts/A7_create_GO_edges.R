library(readr)
library(dplyr)
library(tidyr)

go_df <- read_csv(
  "results/A7_GO_Biological_Process.csv",
  show_col_types = FALSE
)

go_edges <- go_df %>%
  dplyr::select(
    Description,
    geneID
  ) %>%
  tidyr::separate_rows(
    geneID,
    sep = "/"
  ) %>%
  dplyr::transmute(
    protein = geneID,
    biological_process = Description,
    relationship = "INVOLVED_IN"
  )

write_csv(
  go_edges,
  "results/edges_protein_go.csv"
)

cat(
  "\nProtein-GO edges:",
  nrow(go_edges),
  "\n"
)

cat(
  "Unique GO terms:",
  length(unique(go_edges$biological_process)),
  "\n"
)