library(readr)
library(dplyr)

go_df <- read_csv(
  "results/A7_GO_Biological_Process.csv",
  show_col_types = FALSE
)

go_nodes <- go_df %>%
  dplyr::transmute(
    go_term = Description,
    type = "BiologicalProcess",
    p_adjust = p.adjust
  )

write_csv(
  go_nodes,
  "results/nodes_biological_process.csv"
)

cat(
  "\nGO Nodes:",
  nrow(go_nodes),
  "\n"
)