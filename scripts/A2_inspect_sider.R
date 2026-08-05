library(readr)

drug_names <- read_tsv(
  "SIDER/drug_names.tsv",
  col_names = FALSE
)

print(head(drug_names))
print(dim(drug_names))