library(readr)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ReactomePA)

proteins <- read_csv(
  "results/A6_protein_master_list.csv",
  show_col_types = FALSE
)

# Convert symbols -> Entrez IDs
gene_map <- bitr(
  proteins$protein,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

print(head(gene_map))

# Reactome enrichment
result <- enrichPathway(
  gene = gene_map$ENTREZID,
  organism = "human",
  pvalueCutoff = 0.05,
  readable = TRUE
)

result_df <- as.data.frame(result)

write_csv(
  result_df,
  "results/A6_reactome_pathways.csv"
)

cat(
  "\nPathways identified:",
  nrow(result_df),
  "\n"
)