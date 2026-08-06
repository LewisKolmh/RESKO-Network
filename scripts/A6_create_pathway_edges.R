library(readr)
library(dplyr)
library(tidyr)

pathways <- read_csv(
  "results/A6_reactome_pathways.csv",
  show_col_types = FALSE
)

protein_pathway_edges <- pathways %>%
  dplyr::select(
    Description,
    geneID
  )

protein_pathway_edges <- protein_pathway_edges %>%
  tidyr::separate_rows(
    geneID,
    sep = "/"
  )

protein_pathway_edges <- protein_pathway_edges %>%
  dplyr::transmute(
    protein = geneID,
    pathway = Description,
    relationship = "PARTICIPATES_IN"
  )

write_csv(
  protein_pathway_edges,
  "results/edges_protein_pathway.csv"
)

cat(
  "\nProtein-pathway edges:",
  nrow(protein_pathway_edges),
  "\n"
)

cat(
  "Unique pathways:",
  length(unique(protein_pathway_edges$pathway)),
  "\n"
)

cat(
  "Unique proteins:",
  length(unique(protein_pathway_edges$protein)),
  "\n"
)