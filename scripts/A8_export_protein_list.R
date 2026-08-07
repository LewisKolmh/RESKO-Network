library(readr)
library(dplyr)

proteins <- read_csv(
  "results/nodes_proteins.csv",
  show_col_types = FALSE
)

proteins <- proteins %>%
  dplyr::select(protein) %>%
  distinct()

write_csv(
  proteins,
  "results/A8_target_proteins.csv"
)

cat(
  "\nProteins exported:",
  nrow(proteins),
  "\n"
)