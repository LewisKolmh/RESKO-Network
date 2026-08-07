library(clusterProfiler)
library(org.Hs.eg.db)
library(readr)
library(dplyr)

# Load proteins

proteins <- read_csv(
  "results/A6_protein_master_list.csv",
  show_col_types = FALSE
)

# Convert Gene Symbols to Entrez IDs

gene_map <- bitr(
  proteins$protein,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

cat(
  "\nMapped genes:",
  nrow(gene_map),
  "\n"
)

# GO Biological Process enrichment

go_results <- enrichGO(
  gene = gene_map$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  readable = TRUE
)

go_df <- as.data.frame(go_results)

write_csv(
  go_df,
  "results/A7_GO_Biological_Process.csv"
)

cat(
  "\nGO Terms Identified:",
  nrow(go_df),
  "\n"
)