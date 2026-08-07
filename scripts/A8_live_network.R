library(readr)
library(dplyr)
library(visNetwork)
library(htmlwidgets)

# =====================================
# LOAD DATA
# =====================================

proteins <- read_csv(
  "results/nodes_proteins.csv",
  show_col_types = FALSE
)

pathways <- read_csv(
  "results/nodes_pathways.csv",
  show_col_types = FALSE
)

go <- read_csv(
  "results/nodes_biological_process.csv",
  show_col_types = FALSE
)

ppi <- read_csv(
  "results/edges_interacts_with.csv",
  show_col_types = FALSE
)

ppath <- read_csv(
  "results/edges_protein_pathway.csv",
  show_col_types = FALSE
)

pgo <- read_csv(
  "results/edges_protein_go.csv",
  show_col_types = FALSE
)

# =====================================
# CALCULATE PROTEIN DEGREES
# =====================================

protein_degree <- table(
  c(
    ppi$source,
    ppi$target
  )
)

protein_degree <- data.frame(
  protein = names(protein_degree),
  degree = as.numeric(protein_degree)
)

# =====================================
# PROTEIN NODES
# =====================================

protein_nodes <- proteins %>%
  left_join(
    protein_degree,
    by = c("protein" = "protein")
  )

protein_nodes$degree[
  is.na(protein_nodes$degree)
] <- 1

protein_nodes <- data.frame(
  id = protein_nodes$protein,
  label = protein_nodes$protein,
  title = paste0(
    "<b>Protein:</b> ",
    protein_nodes$protein,
    "<br><b>Degree:</b> ",
    protein_nodes$degree
  ),
  group = "Protein",

  # node size scales with connectivity
  size = 15 + (protein_nodes$degree * 2),

  value = protein_nodes$degree
)

# =====================================
# PATHWAY NODES
# =====================================

pathway_nodes <- data.frame(
  id = pathways$pathway,
  label = pathways$pathway,

  title = paste0(
    "<b>Pathway:</b> ",
    pathways$pathway,
    "<br><b>Adjusted p:</b> ",
    signif(pathways$p_adjust, 3)
  ),

  group = "Pathway",
  size = 18
)

# =====================================
# GO NODES
# =====================================

go_nodes <- data.frame(
  id = go$go_term,
  label = go$go_term,

  title = paste0(
    "<b>GO Process:</b> ",
    go$go_term
  ),

  group = "GO",
  size = 12
)

# =====================================
# COMBINE NODES
# =====================================

nodes <- bind_rows(
  protein_nodes,
  pathway_nodes,
  go_nodes
)

# =====================================
# EDGES
# =====================================

ppi_edges <- data.frame(
  from = ppi$source,
  to = ppi$target,
  color = "#999999",
  width = 1 + (ppi$score * 2)
)

pathway_edges <- data.frame(
  from = ppath$protein,
  to = ppath$pathway,
  color = "#1A9850",
  width = 1
)

go_edges <- data.frame(
  from = pgo$protein,
  to = pgo$biological_process,
  color = "#4575B4",
  width = 1
)

edges <- bind_rows(
  ppi_edges,
  pathway_edges,
  go_edges
)

# =====================================
# BUILD NETWORK
# =====================================

network <- visNetwork(
  nodes,
  edges,
  width = "100%",
  height = "1000px",
  main = "EEF1A-Centred Biological Knowledge Graph"
) %>%

  # ---------------------------------
  # NODE GROUPS
  # ---------------------------------

  visGroups(
    groupname = "Protein",
    color = list(
      background = "#D73027",
      border = "#7F0000"
    )
  ) %>%

  visGroups(
    groupname = "Pathway",
    color = list(
      background = "#1A9850",
      border = "#006400"
    )
  ) %>%

  visGroups(
    groupname = "GO",
    color = list(
      background = "#4575B4",
      border = "#003399"
    )
  ) %>%

  # ---------------------------------
  # SEARCHING + FILTERING
  # ---------------------------------

  visOptions(
    highlightNearest = list(
      enabled = TRUE,
      degree = 1
    ),
    nodesIdSelection = TRUE,
    selectedBy = "group"
  ) %>%

  # ---------------------------------
  # INTERACTION
  # ---------------------------------

  visInteraction(
    hover = TRUE,
    navigationButtons = TRUE,
    keyboard = TRUE,
    zoomView = TRUE,
    dragView = TRUE
  ) %>%

  # ---------------------------------
  # EDGES
  # ---------------------------------

  visEdges(
    smooth = FALSE
  ) %>%

  # ---------------------------------
  # LAYOUT STABILISATION
  # ---------------------------------

  visPhysics(
    solver = "forceAtlas2Based",

    forceAtlas2Based = list(
      gravitationalConstant = -75,
      centralGravity = 0.002,
      springLength = 250,
      springConstant = 0.05
    ),

    stabilization = list(
      enabled = TRUE,
      iterations = 3000
    )
  ) %>%

  visLayout(
    randomSeed = 42
  )

# =====================================
# SAVE HTML
# =====================================

saveWidget(
  network,
  "results/eef1a_network_A7.html",
  selfcontained = TRUE
)

network