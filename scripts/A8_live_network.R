# ============================================================
# A8 INTERACTIVE EEF1A NETWORK EXPLORER
# ============================================================

required_packages <- c("readr", "dplyr", "visNetwork", "htmlwidgets", "htmltools")
missing_packages <- required_packages[!vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_packages) > 0) {
  stop(paste0(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall them with:\ninstall.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))"
  ))
}

# ---- Paths ----
protein_nodes_file <- "results/nodes_proteins.csv"
pathway_nodes_file <- "results/nodes_pathways.csv"
go_nodes_file <- "results/nodes_biological_process.csv"
ppi_edges_file <- "results/edges_interacts_with.csv"
pathway_edges_file <- "results/edges_protein_pathway.csv"
go_edges_file <- "results/edges_protein_go.csv"
output_file <- "results/eef1a_network_A8.html"

input_files <- c(
  protein_nodes_file, pathway_nodes_file, go_nodes_file,
  ppi_edges_file, pathway_edges_file, go_edges_file
)
missing_files <- input_files[!file.exists(input_files)]
if (length(missing_files) > 0) {
  stop(paste0(
    "Missing input files:\n", paste(missing_files, collapse = "\n"),
    "\n\nCurrent working directory:\n", getwd(),
    "\n\nRun this script from the RESKO project root."
  ))
}
if (!dir.exists("results")) dir.create("results", recursive = TRUE)

# ---- Load data ----
proteins <- readr::read_csv(protein_nodes_file, show_col_types = FALSE)
pathways <- readr::read_csv(pathway_nodes_file, show_col_types = FALSE)
go_processes <- readr::read_csv(go_nodes_file, show_col_types = FALSE)
ppi <- readr::read_csv(ppi_edges_file, show_col_types = FALSE)
protein_pathway <- readr::read_csv(pathway_edges_file, show_col_types = FALSE)
protein_go <- readr::read_csv(go_edges_file, show_col_types = FALSE)

check_columns <- function(data, required, file_name) {
  missing <- setdiff(required, names(data))
  if (length(missing) > 0) {
    stop(paste0("Missing columns in ", file_name, ": ", paste(missing, collapse = ", ")))
  }
}
check_columns(proteins, "protein", protein_nodes_file)
check_columns(pathways, c("pathway", "p_adjust"), pathway_nodes_file)
check_columns(go_processes, "go_term", go_nodes_file)
check_columns(ppi, c("source", "target", "score"), ppi_edges_file)
check_columns(protein_pathway, c("protein", "pathway"), pathway_edges_file)
check_columns(protein_go, c("protein", "biological_process"), go_edges_file)

# ---- Clean data ----
proteins <- proteins |>
  dplyr::filter(!is.na(protein), protein != "") |>
  dplyr::distinct(protein, .keep_all = TRUE)

pathways <- pathways |>
  dplyr::filter(!is.na(pathway), pathway != "") |>
  dplyr::distinct(pathway, .keep_all = TRUE)

go_processes <- go_processes |>
  dplyr::filter(!is.na(go_term), go_term != "") |>
  dplyr::distinct(go_term, .keep_all = TRUE)

ppi <- ppi |>
  dplyr::filter(!is.na(source), !is.na(target), source != "", target != "") |>
  dplyr::mutate(
    score = suppressWarnings(as.numeric(score)),
    score = dplyr::coalesce(score, 0)
  ) |>
  dplyr::distinct(source, target, .keep_all = TRUE)

protein_pathway <- protein_pathway |>
  dplyr::filter(!is.na(protein), !is.na(pathway), protein != "", pathway != "") |>
  dplyr::distinct(protein, pathway, .keep_all = TRUE)

protein_go <- protein_go |>
  dplyr::filter(!is.na(protein), !is.na(biological_process), protein != "", biological_process != "") |>
  dplyr::distinct(protein, biological_process, .keep_all = TRUE)

# ---- Connectivity metrics ----
ppi_degree <- table(c(ppi$source, ppi$target))
ppi_degree <- data.frame(
  protein = names(ppi_degree),
  ppi_degree = as.numeric(ppi_degree),
  stringsAsFactors = FALSE
)
pathway_degree <- protein_pathway |> dplyr::count(protein, name = "pathway_degree")
go_degree <- protein_go |> dplyr::count(protein, name = "go_degree")

protein_metrics <- proteins |>
  dplyr::left_join(ppi_degree, by = "protein") |>
  dplyr::left_join(pathway_degree, by = "protein") |>
  dplyr::left_join(go_degree, by = "protein") |>
  dplyr::mutate(
    ppi_degree = dplyr::coalesce(ppi_degree, 0),
    pathway_degree = dplyr::coalesce(pathway_degree, 0L),
    go_degree = dplyr::coalesce(go_degree, 0L),
    total_degree = ppi_degree + pathway_degree + go_degree,
    node_size = 18 + 3.2 * sqrt(total_degree)
  )

pathway_metrics <- protein_pathway |> dplyr::count(pathway, name = "protein_count")
pathways <- pathways |>
  dplyr::left_join(pathway_metrics, by = "pathway") |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 10 + 2.3 * sqrt(protein_count)
  )

go_metrics <- protein_go |> dplyr::count(biological_process, name = "protein_count")
go_processes <- go_processes |>
  dplyr::left_join(go_metrics, by = c("go_term" = "biological_process")) |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 9 + 2.1 * sqrt(protein_count)
  )

# ---- Nodes ----
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
      "<strong>Adjusted p-value:</strong> ", format(p_adjust, scientific = TRUE, digits = 3), "<br>",
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

nodes <- dplyr::bind_rows(protein_nodes, pathway_nodes, go_nodes)
if (anyDuplicated(nodes$id) > 0) {
  stop(paste0("Duplicate node IDs detected:\n", paste(unique(nodes$id[duplicated(nodes$id)]), collapse = "\n")))
}

# ---- Edges ----
ppi_edges <- ppi |>
  dplyr::transmute(
    from = paste0("protein::", source),
    to = paste0("protein::", target),
    relationship = "INTERACTS_WITH",
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;'>",
      "<strong>Protein interaction</strong><br>", source, " - ", target, "<br>",
      "<strong>STRING score:</strong> ", round(score, 3), "</div>"
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
      "<strong>Pathway association</strong><br>", protein, " - ", pathway, "</div>"
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
      "<strong>GO association</strong><br>", protein, " - ", biological_process, "</div>"
    ),
    color = "#5B8FD1",
    width = 1.05,
    dashes = FALSE,
    arrows = ""
  )

edges <- dplyr::bind_rows(ppi_edges, pathway_edges, go_edges)

missing_source_nodes <- setdiff(unique(edges$from), nodes$id)
missing_target_nodes <- setdiff(unique(edges$to), nodes$id)
if (length(missing_source_nodes) > 0 || length(missing_target_nodes) > 0) {
  stop(paste0(
    "Some edge endpoints do not have matching nodes.",
    "\n\nMissing source nodes:\n", paste(missing_source_nodes, collapse = "\n"),
    "\n\nMissing target nodes:\n", paste(missing_target_nodes, collapse = "\n")
  ))
}

# ---- Network ----
network <- visNetwork::visNetwork(
  nodes = nodes,
  edges = edges,
  width = "100%",
  height = "900px",
  main = "EEF1A Network Explorer"
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
    shadow = list(enabled = TRUE, color = "rgba(0,0,0,0.14)", size = 5, x = 2, y = 2),
    font = list(face = "Calibri", strokeWidth = 3, strokeColor = "#FFFFFF")
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
      gravitationalConstant = -95,
      centralGravity = 0.006,
      springLength = 250,
      springConstant = 0.035,
      damping = 0.62,
      avoidOverlap = 0.75
    ),
    stabilization = list(
      enabled = TRUE,
      iterations = 3000,
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

# ---- Reliable legend above graph ----
legend_panel <- htmltools::tags$div(
  style = paste0(
    "box-sizing:border-box;width:100%;margin:0 0 10px 0;padding:10px 14px;",
    "background:#F8F9FB;border:1px solid #D7DADE;border-radius:5px;",
    "font-family:Calibri,Arial,sans-serif;font-size:13px;color:#222222;",
    "display:flex;flex-wrap:wrap;align-items:center;gap:10px 22px;"
  ),
  htmltools::tags$span(
    style = "font-weight:600;font-size:14px;margin-right:4px;",
    "Network key"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:14px;height:14px;border-radius:50%;background:#D95F59;border:1px solid #8E2F2A;margin-right:7px;"),
    "Protein"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:14px;height:14px;border-radius:50%;background:#4E9F6E;border:1px solid #25633A;margin-right:7px;"),
    "Reactome pathway"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:14px;height:14px;border-radius:50%;background:#5B8FD1;border:1px solid #2B5797;margin-right:7px;"),
    "GO biological process"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:24px;height:0;border-top:3px solid #9A9A9A;margin-right:7px;"),
    "Protein interaction"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:24px;height:0;border-top:3px solid #40916C;margin-right:7px;"),
    "Pathway association"
  ),
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(style = "display:inline-block;width:24px;height:0;border-top:3px solid #5B8FD1;margin-right:7px;"),
    "GO association"
  ),
  htmltools::tags$span(
    style = "color:#666666;font-size:12px;margin-left:auto;",
    "Select a node to fade unrelated elements"
  )
)

# ---- Page styling ----
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

# ---- Save ----
htmlwidgets::saveWidget(
  widget = network,
  file = output_file,
  selfcontained = TRUE
)

# ---- Summary ----
cat("\n")
cat("EEF1A network created successfully.\n")
cat("-----------------------------------\n")
cat("Protein nodes:             ", nrow(protein_nodes), "\n")
cat("Reactome pathway nodes:    ", nrow(pathway_nodes), "\n")
cat("GO biological processes:   ", nrow(go_nodes), "\n")
cat("Protein interaction edges: ", nrow(ppi_edges), "\n")
cat("Protein-pathway edges:     ", nrow(pathway_edges), "\n")
cat("Protein-GO edges:          ", nrow(go_edges), "\n")
cat("Total nodes:               ", nrow(nodes), "\n")
cat("Total edges:               ", nrow(edges), "\n")
cat("Output file:               ", output_file, "\n\n")

# ---- Open only the saved file ----
output_path <- normalizePath(
  output_file,
  winslash = "/",
  mustWork = TRUE
)

utils::browseURL(output_path)
cat("Opened saved network:\n", output_path, "\n")
