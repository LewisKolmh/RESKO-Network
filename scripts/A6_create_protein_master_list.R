library(readr)
library(dplyr)

nodes <- read_csv(
  "results/nodes_proteins.csv",
  show_col_types = FALSE
)

proteins <- nodes %>%
  distinct(protein)

write_csv(
  proteins,
  "results/A6_protein_master_list.csv"
)

cat(
  "Unique proteins:",
  nrow(proteins),
  "\n"
)