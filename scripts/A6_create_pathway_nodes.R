library(readr)
library(dplyr)

pathways <- read_csv(
  "results/A6_reactome_pathways.csv",
  show_col_types = FALSE
)

pathway_nodes <- pathways %>%
  transmute(
    pathway = Description,
    type = "Pathway",
    p_adjust = p.adjust
  )

write_csv(
  pathway_nodes,
  "results/nodes_pathways.csv"
)

cat(
  "Pathway nodes:",
  nrow(pathway_nodes),
  "\n"
)