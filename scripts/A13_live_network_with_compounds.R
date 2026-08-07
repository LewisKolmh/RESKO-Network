# ============================================================
# A13 INTERACTIVE EEF1A NETWORK WITH COMPOUND LAYER
#
# Node types:
#   - Compound
#   - Protein
#   - Reactome pathway
#   - GO biological process
#
# Compound-protein relationship types:
#   - INHIBITS
#   - BINDS_TO
#   - HAS_FUNCTIONAL_ACTIVITY_AGAINST
#   - HAS_ACTIVITY_AGAINST
#
# Output:
#   results/eef1a_network_A13_compounds.html
# ============================================================

# ============================================================
# 1. CHECK REQUIRED PACKAGES
# ============================================================

required_packages <- c(
  "readr",
  "dplyr",
  "visNetwork",
  "htmlwidgets",
  "htmltools"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0) {
  stop(
    paste0(
      "Missing packages: ",
      paste(missing_packages, collapse = ", "),
      "\n\nInstall them with:\n",
      "install.packages(c(",
      paste0('"', missing_packages, '"', collapse = ", "),
      "))"
    )
  )
}

# ============================================================
# 2. FILE PATHS
# ============================================================

protein_nodes_file <- "results/nodes_proteins.csv"
pathway_nodes_file <- "results/nodes_pathways.csv"
go_nodes_file <- "results/nodes_biological_process.csv"
compound_nodes_file <- "results/nodes_drugs.csv"

ppi_edges_file <- "results/edges_interacts_with.csv"
pathway_edges_file <- "results/edges_protein_pathway.csv"
go_edges_file <- "results/edges_protein_go.csv"
compound_edges_file <- "results/edges_compound_protein.csv"

output_file <- "results/eef1a_network_A13_compounds.html"

input_files <- c(
  protein_nodes_file,
  pathway_nodes_file,
  go_nodes_file,
  compound_nodes_file,
  ppi_edges_file,
  pathway_edges_file,
  go_edges_file,
  compound_edges_file
)

missing_files <- input_files[!file.exists(input_files)]

if (length(missing_files) > 0) {
  stop(
    paste0(
      "Missing input files:\n",
      paste(missing_files, collapse = "\n"),
      "\n\nComplete A11 and A12 before running A13.",
      "\nCurrent working directory:\n",
      getwd()
    )
  )
}

if (!dir.exists("results")) {
  dir.create("results", recursive = TRUE)
}

# ============================================================
# 3. LOAD DATA
# ============================================================

proteins <- readr::read_csv(
  protein_nodes_file,
  show_col_types = FALSE,
  progress = FALSE
)

pathways <- readr::read_csv(
  pathway_nodes_file,
  show_col_types = FALSE,
  progress = FALSE
)

go_processes <- readr::read_csv(
  go_nodes_file,
  show_col_types = FALSE,
  progress = FALSE
)

compounds <- readr::read_csv(
  compound_nodes_file,
  show_col_types = FALSE,
  progress = FALSE
)

ppi <- readr::read_csv(
  ppi_edges_file,
  show_col_types = FALSE,
  progress = FALSE
)

protein_pathway <- readr::read_csv(
  pathway_edges_file,
  show_col_types = FALSE,
  progress = FALSE
)

protein_go <- readr::read_csv(
  go_edges_file,
  show_col_types = FALSE,
  progress = FALSE
)

compound_protein <- readr::read_csv(
  compound_edges_file,
  show_col_types = FALSE,
  progress = FALSE
)

# ============================================================
# 4. VALIDATE REQUIRED COLUMNS
# ============================================================

check_columns <- function(data, required_columns, file_name) {
  missing_columns <- setdiff(required_columns, names(data))

  if (length(missing_columns) > 0) {
    stop(
      paste0(
        "Missing columns in ",
        file_name,
        ": ",
        paste(missing_columns, collapse = ", ")
      )
    )
  }
}

check_columns(proteins, "protein", protein_nodes_file)
check_columns(pathways, c("pathway", "p_adjust"), pathway_nodes_file)
check_columns(go_processes, "go_term", go_nodes_file)
check_columns(
  compounds,
  c(
    "molecule_chembl_id",
    "name",
    "molecule_type",
    "development_status",
    "proteins",
    "activity_types",
    "strongest_record_nM",
    "total_measurements"
  ),
  compound_nodes_file
)
check_columns(ppi, c("source", "target", "score"), ppi_edges_file)
check_columns(protein_pathway, c("protein", "pathway"), pathway_edges_file)
check_columns(protein_go, c("protein", "biological_process"), go_edges_file)
check_columns(
  compound_protein,
  c(
    "molecule_chembl_id",
    "compound_name",
    "protein",
    "relationship",
    "activity_types",
    "minimum_activity_nM",
    "total_measurements",
    "development_status"
  ),
  compound_edges_file
)

# ============================================================
# 5. CLEAN DATA
# ============================================================

proteins <- proteins |>
  dplyr::filter(!is.na(protein), protein != "") |>
  dplyr::distinct(protein, .keep_all = TRUE)

pathways <- pathways |>
  dplyr::filter(!is.na(pathway), pathway != "") |>
  dplyr::distinct(pathway, .keep_all = TRUE)

go_processes <- go_processes |>
  dplyr::filter(!is.na(go_term), go_term != "") |>
  dplyr::distinct(go_term, .keep_all = TRUE)

compounds <- compounds |>
  dplyr::filter(
    !is.na(molecule_chembl_id),
    molecule_chembl_id != ""
  ) |>
  dplyr::mutate(
    name = dplyr::if_else(
      is.na(name) | name == "",
      molecule_chembl_id,
      name
    ),
    strongest_record_nM = suppressWarnings(
      as.numeric(strongest_record_nM)
    ),
    total_measurements = suppressWarnings(
      as.integer(total_measurements)
    )
  ) |>
  dplyr::distinct(molecule_chembl_id, .keep_all = TRUE)

ppi <- ppi |>
  dplyr::filter(
    !is.na(source),
    !is.na(target),
    source != "",
    target != ""
  ) |>
  dplyr::mutate(
    score = suppressWarnings(as.numeric(score)),
    score = dplyr::coalesce(score, 0)
  ) |>
  dplyr::distinct(source, target, .keep_all = TRUE)

protein_pathway <- protein_pathway |>
  dplyr::filter(
    !is.na(protein),
    !is.na(pathway),
    protein != "",
    pathway != ""
  ) |>
  dplyr::distinct(protein, pathway, .keep_all = TRUE)

protein_go <- protein_go |>
  dplyr::filter(
    !is.na(protein),
    !is.na(biological_process),
    protein != "",
    biological_process != ""
  ) |>
  dplyr::distinct(protein, biological_process, .keep_all = TRUE)

compound_protein <- compound_protein |>
  dplyr::filter(
    !is.na(molecule_chembl_id),
    !is.na(protein),
    molecule_chembl_id != "",
    protein != ""
  ) |>
  dplyr::mutate(
    compound_name = dplyr::if_else(
      is.na(compound_name) | compound_name == "",
      molecule_chembl_id,
      compound_name
    ),
    minimum_activity_nM = suppressWarnings(
      as.numeric(minimum_activity_nM)
    ),
    total_measurements = suppressWarnings(
      as.integer(total_measurements)
    )
  ) |>
  dplyr::distinct(
    molecule_chembl_id,
    protein,
    .keep_all = TRUE
  )

# ============================================================
# 6. CONNECTIVITY METRICS
# ============================================================

ppi_degree <- table(c(ppi$source, ppi$target))
ppi_degree <- data.frame(
  protein = names(ppi_degree),
  ppi_degree = as.numeric(ppi_degree),
  stringsAsFactors = FALSE
)

pathway_degree <- protein_pathway |>
  dplyr::count(protein, name = "pathway_degree")

go_degree <- protein_go |>
  dplyr::count(protein, name = "go_degree")

compound_degree <- compound_protein |>
  dplyr::count(protein, name = "compound_degree")

protein_metrics <- proteins |>
  dplyr::left_join(ppi_degree, by = "protein") |>
  dplyr::left_join(pathway_degree, by = "protein") |>
  dplyr::left_join(go_degree, by = "protein") |>
  dplyr::left_join(compound_degree, by = "protein") |>
  dplyr::mutate(
    ppi_degree = dplyr::coalesce(ppi_degree, 0),
    pathway_degree = dplyr::coalesce(pathway_degree, 0L),
    go_degree = dplyr::coalesce(go_degree, 0L),
    compound_degree = dplyr::coalesce(compound_degree, 0L),
    total_degree =
      ppi_degree +
      pathway_degree +
      go_degree +
      compound_degree,
    node_size = 18 + 3.0 * sqrt(total_degree)
  )

pathway_metrics <- protein_pathway |>
  dplyr::count(pathway, name = "protein_count")

pathways <- pathways |>
  dplyr::left_join(pathway_metrics, by = "pathway") |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 10 + 2.3 * sqrt(protein_count)
  )

go_metrics <- protein_go |>
  dplyr::count(biological_process, name = "protein_count")

go_processes <- go_processes |>
  dplyr::left_join(
    go_metrics,
    by = c("go_term" = "biological_process")
  ) |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 9 + 2.1 * sqrt(protein_count)
  )

compound_target_counts <- compound_protein |>
  dplyr::count(molecule_chembl_id, name = "target_count")

compounds <- compounds |>
  dplyr::left_join(
    compound_target_counts,
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    target_count = dplyr::coalesce(target_count, 0L),
    node_size = 16 + 3.5 * sqrt(target_count + total_measurements)
  )

# ============================================================
# 7. CREATE NODES
# ============================================================

protein_nodes <- protein_metrics |>
  dplyr::transmute(
    id = paste0("protein::", protein),
    label = protein,
    group = "Protein",
    shape = "dot",
    size = node_size,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>Protein</strong><br>", protein, "<br><br>",
      "<strong>PPI connections:</strong> ", ppi_degree, "<br>",
      "<strong>Compound connections:</strong> ", compound_degree, "<br>",
      "<strong>Pathway connections:</strong> ", pathway_degree, "<br>",
      "<strong>GO connections:</strong> ", go_degree, "<br>",
      "<strong>Total connectivity:</strong> ", total_degree,
      "</div>"
    )
  )

pathway_nodes <- pathways |>
  dplyr::transmute(
    id = paste0("pathway::", pathway),
    label = pathway,
    group = "Pathway",
    shape = "dot",
    size = node_size,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>Reactome pathway</strong><br>", pathway, "<br><br>",
      "<strong>Adjusted p-value:</strong> ",
      format(p_adjust, scientific = TRUE, digits = 3), "<br>",
      "<strong>Associated proteins:</strong> ", protein_count,
      "</div>"
    )
  )

go_nodes <- go_processes |>
  dplyr::transmute(
    id = paste0("go::", go_term),
    label = go_term,
    group = "BiologicalProcess",
    shape = "dot",
    size = node_size,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>GO biological process</strong><br>", go_term, "<br><br>",
      "<strong>Associated proteins:</strong> ", protein_count,
      "</div>"
    )
  )

compound_nodes <- compounds |>
  dplyr::transmute(
    id = paste0("compound::", molecule_chembl_id),
    label = name,
    group = "Compound",
    shape = "diamond",
    size = node_size,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>Compound</strong><br>", name, "<br>",
      "<strong>ChEMBL ID:</strong> ", molecule_chembl_id, "<br><br>",
      "<strong>Molecule type:</strong> ", dplyr::coalesce(molecule_type, "Unknown"), "<br>",
      "<strong>Development status:</strong> ", dplyr::coalesce(development_status, "Unknown"), "<br>",
      "<strong>Associated proteins:</strong> ", dplyr::coalesce(proteins, "Unknown"), "<br>",
      "<strong>Activity types:</strong> ", dplyr::coalesce(activity_types, "Unknown"), "<br>",
      "<strong>Strongest quantitative record:</strong> ",
      ifelse(
        is.na(strongest_record_nM),
        "Unknown",
        paste0(signif(strongest_record_nM, 4), " nM")
      ), "<br>",
      "<strong>Measurements:</strong> ", dplyr::coalesce(total_measurements, 0L),
      "</div>"
    )
  )

nodes <- dplyr::bind_rows(
  protein_nodes,
  pathway_nodes,
  go_nodes,
  compound_nodes
)

if (anyDuplicated(nodes$id) > 0) {
  stop(
    paste0(
      "Duplicate node IDs detected:\n",
      paste(unique(nodes$id[duplicated(nodes$id)]), collapse = "\n")
    )
  )
}

# ============================================================
# 8. CREATE EDGES
# ============================================================

ppi_edges <- ppi |>
  dplyr::transmute(
    from = paste0("protein::", source),
    to = paste0("protein::", target),
    relationship = "INTERACTS_WITH",
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;'>",
      "<strong>Protein interaction</strong><br>", source, " - ", target, "<br>",
      "<strong>STRING score:</strong> ", round(score, 3),
      "</div>"
    ),
    color = "#9A9A9A",
    width = 0.9 + 2.7 * score,
    dashes = FALSE,
    arrows = ""
  )

pathway_edges <- protein_pathway |>
  dplyr::transmute(
    from = paste0("protein::", protein),
    to = paste0("pathway::", pathway),
    relationship = "PARTICIPATES_IN",
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;'>",
      "<strong>Pathway association</strong><br>", protein, " - ", pathway,
      "</div>"
    ),
    color = "#40916C",
    width = 1.15,
    dashes = FALSE,
    arrows = ""
  )

go_edges <- protein_go |>
  dplyr::transmute(
    from = paste0("protein::", protein),
    to = paste0("go::", biological_process),
    relationship = "INVOLVED_IN",
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;'>",
      "<strong>GO association</strong><br>", protein, " - ", biological_process,
      "</div>"
    ),
    color = "#5B8FD1",
    width = 1.05,
    dashes = FALSE,
    arrows = ""
  )

compound_edges <- compound_protein |>
  dplyr::mutate(
    edge_color = dplyr::case_when(
      relationship == "INHIBITS" ~ "#C43C39",
      relationship == "BINDS_TO" ~ "#8E63B0",
      relationship == "HAS_FUNCTIONAL_ACTIVITY_AGAINST" ~ "#C99720",
      TRUE ~ "#C77729"
    ),
    edge_width = dplyr::case_when(
      relationship == "INHIBITS" ~ 3.0,
      relationship == "BINDS_TO" ~ 2.4,
      relationship == "HAS_FUNCTIONAL_ACTIVITY_AGAINST" ~ 2.0,
      TRUE ~ 1.7
    )
  ) |>
  dplyr::transmute(
    from = paste0("compound::", molecule_chembl_id),
    to = paste0("protein::", protein),
    relationship,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>", relationship, "</strong><br>",
      compound_name, " - ", protein, "<br>",
      "<strong>Activity types:</strong> ", activity_types, "<br>",
      "<strong>Minimum quantitative activity:</strong> ",
      ifelse(
        is.na(minimum_activity_nM),
        "Unknown",
        paste0(signif(minimum_activity_nM, 4), " nM")
      ), "<br>",
      "<strong>Measurements:</strong> ", dplyr::coalesce(total_measurements, 0L), "<br>",
      "<strong>Development status:</strong> ", dplyr::coalesce(development_status, "Unknown"),
      "</div>"
    ),
    color = edge_color,
    width = edge_width,
    dashes = FALSE,
    arrows = ""
  )

edges <- dplyr::bind_rows(
  ppi_edges,
  pathway_edges,
  go_edges,
  compound_edges
)

missing_source_nodes <- setdiff(unique(edges$from), nodes$id)
missing_target_nodes <- setdiff(unique(edges$to), nodes$id)

if (length(missing_source_nodes) > 0 || length(missing_target_nodes) > 0) {
  stop(
    paste0(
      "Some edge endpoints do not have matching nodes.",
      "\n\nMissing source nodes:\n",
      paste(missing_source_nodes, collapse = "\n"),
      "\n\nMissing target nodes:\n",
      paste(missing_target_nodes, collapse = "\n")
    )
  )
}

# ============================================================
# 9. BUILD NETWORK
# ============================================================

network <- visNetwork::visNetwork(
  nodes = nodes,
  edges = edges,
  width = "100%",
  height = "900px",
  main = "EEF1A Compound-Protein Network Explorer"
) |>
  visNetwork::visGroups(
    groupname = "Compound",
    shape = "diamond",
    color = list(
      background = "#E58B3A",
      border = "#9A4F12",
      highlight = list(background = "#F4AF6A", border = "#7A3D0D"),
      hover = list(background = "#F4AF6A", border = "#7A3D0D")
    ),
    font = list(face = "Calibri", size = 15, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "Protein",
    shape = "dot",
    color = list(
      background = "#D95F59",
      border = "#8E2F2A",
      highlight = list(background = "#FF8A80", border = "#641D19"),
      hover = list(background = "#F28E85", border = "#641D19")
    ),
    font = list(face = "Calibri", size = 16, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "Pathway",
    shape = "dot",
    color = list(
      background = "#4E9F6E",
      border = "#25633A",
      highlight = list(background = "#78BE92", border = "#174728"),
      hover = list(background = "#78BE92", border = "#174728")
    ),
    font = list(face = "Calibri", size = 14, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "BiologicalProcess",
    shape = "dot",
    color = list(
      background = "#5B8FD1",
      border = "#2B5797",
      highlight = list(background = "#83AEE0", border = "#1E3E70"),
      hover = list(background = "#83AEE0", border = "#1E3E70")
    ),
    font = list(face = "Calibri", size = 13, color = "#202020")
  ) |>
  visNetwork::visNodes(
    borderWidth = 1.5,
    borderWidthSelected = 3,
    shadow = list(
      enabled = TRUE,
      color = "rgba(0,0,0,0.14)",
      size = 5,
      x = 2,
      y = 2
    ),
    font = list(
      face = "Calibri",
      strokeWidth = 3,
      strokeColor = "#FFFFFF"
    )
  ) |>
  visNetwork::visEdges(
    arrows = "",
    smooth = list(enabled = FALSE),
    selectionWidth = 3,
    hoverWidth = 2,
    font = list(face = "Calibri", size = 12)
  ) |>
  visNetwork::visOptions(
    highlightNearest = list(
      enabled = TRUE,
      degree = 1,
      hover = FALSE,
      algorithm = "all"
    ),
    nodesIdSelection = list(
      enabled = TRUE,
      useLabels = TRUE,
      main = "Find a node"
    ),
    selectedBy = list(
      variable = "group",
      main = "Filter by node type",
      multiple = TRUE
    )
  ) |>
  visNetwork::visInteraction(
    hover = TRUE,
    navigationButtons = TRUE,
    keyboard = TRUE,
    zoomView = TRUE,
    dragView = TRUE,
    dragNodes = TRUE,
    hideEdgesOnDrag = TRUE,
    tooltipDelay = 150
  ) |>
  visNetwork::visPhysics(
    enabled = TRUE,
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(
      gravitationalConstant = -105,
      centralGravity = 0.006,
      springLength = 260,
      springConstant = 0.033,
      damping = 0.64,
      avoidOverlap = 0.8
    ),
    stabilization = list(
      enabled = TRUE,
      iterations = 3500,
      updateInterval = 50,
      onlyDynamicEdges = FALSE,
      fit = TRUE
    ),
    minVelocity = 0.7,
    timestep = 0.35,
    adaptiveTimestep = TRUE
  ) |>
  visNetwork::visLayout(
    randomSeed = 42,
    improvedLayout = TRUE
  ) |>
  visNetwork::visEvents(
    stabilizationIterationsDone = "
      function () {
        this.setOptions({physics: {enabled: false}});
        this.fit({animation: {duration: 700, easingFunction: 'easeInOutQuad'}});
      }
    "
  )

# ============================================================
# 10. CREATE UPDATED LEGEND
#
# The legend explicitly includes all node and edge types that
# are present in the graph. It is placed above the graph so it
# cannot cover any network elements.
# ============================================================

node_key <- function(color, border, label, shape = "circle") {
  radius <- if (shape == "diamond") "2px" else "50%"
  transform <- if (shape == "diamond") "transform:rotate(45deg);" else ""

  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(
      style = paste0(
        "display:inline-block;",
        "width:14px;",
        "height:14px;",
        "border-radius:", radius, ";",
        "background:", color, ";",
        "border:1px solid ", border, ";",
        "margin-right:8px;",
        transform
      )
    ),
    label
  )
}

edge_key <- function(color, label, width = 3) {
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(
      style = paste0(
        "display:inline-block;",
        "width:25px;",
        "height:0;",
        "border-top:", width, "px solid ", color, ";",
        "margin-right:8px;"
      )
    ),
    label
  )
}

legend_panel <- htmltools::tags$div(
  style = paste0(
    "box-sizing:border-box;",
    "width:100%;",
    "margin:0 0 10px 0;",
    "padding:11px 14px;",
    "background:#F8F9FB;",
    "border:1px solid #D7DADE;",
    "border-radius:5px;",
    "font-family:Calibri,Arial,sans-serif;",
    "font-size:13px;",
    "color:#222222;"
  ),
  htmltools::tags$div(
    style = "display:flex;flex-wrap:wrap;align-items:center;gap:10px 22px;",
    htmltools::tags$span(
      style = "font-weight:600;font-size:14px;margin-right:3px;",
      "Node types"
    ),
    node_key("#E58B3A", "#9A4F12", "Compound", "diamond"),
    node_key("#D95F59", "#8E2F2A", "Protein"),
    node_key("#4E9F6E", "#25633A", "Reactome pathway"),
    node_key("#5B8FD1", "#2B5797", "GO biological process")
  ),
  htmltools::tags$div(
    style = "border-top:1px solid #D7DADE;margin:10px 0;"
  ),
  htmltools::tags$div(
    style = "display:flex;flex-wrap:wrap;align-items:center;gap:10px 22px;",
    htmltools::tags$span(
      style = "font-weight:600;font-size:14px;margin-right:3px;",
      "Relationship types"
    ),
    edge_key("#C43C39", "Inhibitory evidence"),
    edge_key("#8E63B0", "Binding evidence"),
    edge_key("#C99720", "Functional-response evidence"),
    edge_key("#C77729", "Other compound activity"),
    edge_key("#9A9A9A", "Protein interaction"),
    edge_key("#40916C", "Pathway association"),
    edge_key("#5B8FD1", "GO association"),
    htmltools::tags$span(
      style = "color:#666666;font-size:12px;margin-left:auto;",
      "Select a node to fade unrelated elements"
    )
  )
)

# ============================================================
# 11. ADD PAGE STYLING AND LEGEND
# ============================================================

page_style <- htmltools::tags$style(
  htmltools::HTML(
    "
    body {
      margin: 0;
      padding: 14px;
      background: #FFFFFF;
      font-family: Calibri, Arial, sans-serif;
    }

    h1, h2, h3 {
      margin-top: 4px;
      margin-bottom: 10px;
      color: #222222;
      font-family: Calibri, Arial, sans-serif;
      font-weight: 600;
    }

    select, input, button, label {
      font-family: Calibri, Arial, sans-serif !important;
      font-size: 13px !important;
    }

    .vis-network {
      border: 1px solid #E0E0E0;
      border-radius: 5px;
      background: #FFFFFF;
    }
    "
  )
)

network <- htmlwidgets::prependContent(
  network,
  page_style,
  legend_panel
)

# ============================================================
# 12. SAVE, VERIFY AND OPEN
# ============================================================

htmlwidgets::saveWidget(
  widget = network,
  file = output_file,
  selfcontained = TRUE
)

if (!file.exists(output_file)) {
  stop("The A13 network HTML file was not created.")
}

if (is.na(file.info(output_file)$size) || file.info(output_file)$size <= 0) {
  stop("The A13 network HTML file is empty or invalid.")
}

cat("\n")
cat("A13 compound network created successfully.\n")
cat("------------------------------------------\n")
cat("Compound nodes:             ", nrow(compound_nodes), "\n")
cat("Protein nodes:              ", nrow(protein_nodes), "\n")
cat("Reactome pathway nodes:     ", nrow(pathway_nodes), "\n")
cat("GO biological processes:    ", nrow(go_nodes), "\n")
cat("Compound-protein edges:     ", nrow(compound_edges), "\n")
cat("Protein interaction edges:  ", nrow(ppi_edges), "\n")
cat("Protein-pathway edges:      ", nrow(pathway_edges), "\n")
cat("Protein-GO edges:           ", nrow(go_edges), "\n")
cat("Total nodes:                ", nrow(nodes), "\n")
cat("Total edges:                ", nrow(edges), "\n")
cat("Output file:                ", output_file, "\n\n")

output_path <- normalizePath(
  output_file,
  winslash = "/",
  mustWork = TRUE
)

utils::browseURL(output_path)
cat("Opened saved network:\n", output_path, "\n")
