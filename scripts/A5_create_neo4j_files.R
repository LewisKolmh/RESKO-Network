library(readr)
library(dplyr)

# Load STRING interactions
interactions <- read_csv(
  "results/A4_eef1a_string_interactions.csv",
  show_col_types = FALSE
)

# -------------------------
# Create protein nodes
# -------------------------

protein_nodes <- unique(
  c(
    interactions$preferredName_A,
    interactions$preferredName_B
  )
)

protein_nodes <- data.frame(
  protein = protein_nodes,
  type = "Protein"
)

# -------------------------
# Create protein-protein edges
# -------------------------

protein_edges <- interactions %>%
  transmute(
    source = preferredName_A,
    target = preferredName_B,
    score = score,
    relationship = "INTERACTS_WITH"
  )

# -------------------------
# Save files
# -------------------------

write_csv(
  protein_nodes,
  "results/nodes_proteins.csv"
)

write_csv(
  protein_edges,
  "results/edges_interacts_with.csv"
)

cat("\nUnique protein nodes:", nrow(protein_nodes), "\n")
cat("Interaction edges:", nrow(protein_edges), "\n")
cat("\nNeo4j files created successfully.\n")