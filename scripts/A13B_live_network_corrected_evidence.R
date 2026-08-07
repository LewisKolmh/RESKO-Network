# ============================================================
# A13B LIVE NETWORK WITH PROVENANCE-CORRECTED COMPOUND EVIDENCE
#
# Inputs:
#   results/nodes_proteins.csv
#   results/nodes_pathways.csv
#   results/nodes_biological_process.csv
#   results/nodes_drugs.csv
#   results/edges_interacts_with.csv
#   results/edges_protein_pathway.csv
#   results/edges_protein_go.csv
#   results/edges_compound_protein_corrected.csv
#
# Output:
#   results/eef1a_network_A13B_corrected.html
# ============================================================

required_packages <- c("readr", "dplyr", "visNetwork", "htmlwidgets", "htmltools")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop(paste0(
    "Missing packages: ", paste(missing_packages, collapse = ", "),
    "\nInstall with: install.packages(c(",
    paste0('"', missing_packages, '"', collapse = ", "), "))"
  ))
}

files <- list(
  proteins = "results/nodes_proteins.csv",
  pathways = "results/nodes_pathways.csv",
  go = "results/nodes_biological_process.csv",
  compounds = "results/nodes_drugs.csv",
  ppi = "results/edges_interacts_with.csv",
  protein_pathway = "results/edges_protein_pathway.csv",
  protein_go = "results/edges_protein_go.csv",
  compound_protein = "results/edges_compound_protein_corrected.csv"
)
output_file <- "results/eef1a_network_A13B_corrected.html"

missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) {
  stop(paste0(
    "Missing inputs:\n", paste(missing_files, collapse = "\n"),
    "\nRun A12B before A13B."
  ))
}

read_data <- function(path) {
  readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
}

proteins <- read_data(files$proteins)
pathways <- read_data(files$pathways)
go_processes <- read_data(files$go)
compounds <- read_data(files$compounds)
ppi <- read_data(files$ppi)
protein_pathway <- read_data(files$protein_pathway)
protein_go <- read_data(files$protein_go)
compound_protein <- read_data(files$compound_protein)

check_columns <- function(data, columns, name) {
  missing <- setdiff(columns, names(data))
  if (length(missing) > 0) {
    stop(paste0("Missing columns in ", name, ": ", paste(missing, collapse = ", ")))
  }
}

check_columns(proteins, "protein", files$proteins)
check_columns(pathways, c("pathway", "p_adjust"), files$pathways)
check_columns(go_processes, "go_term", files$go)
check_columns(compounds, c(
  "molecule_chembl_id", "name", "molecule_type", "development_status",
  "proteins", "activity_types", "strongest_record_nM", "total_measurements"
), files$compounds)
check_columns(ppi, c("source", "target", "score"), files$ppi)
check_columns(protein_pathway, c("protein", "pathway"), files$protein_pathway)
check_columns(protein_go, c("protein", "biological_process"), files$protein_go)
check_columns(compound_protein, c(
  "molecule_chembl_id", "display_name", "queried_protein", "relationship",
  "corrected_activity_types", "original_activity_types", "minimum_activity_nM",
  "raw_record_count", "unique_assay_count", "unique_document_count",
  "confidence_score", "duplicate_kd_ed50_collapsed", "development_status"
), files$compound_protein)

# Clean and deduplicate.
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
  dplyr::filter(!is.na(molecule_chembl_id), molecule_chembl_id != "") |>
  dplyr::mutate(
    name = dplyr::if_else(is.na(name) | name == "", molecule_chembl_id, name),
    strongest_record_nM = suppressWarnings(as.numeric(strongest_record_nM)),
    total_measurements = suppressWarnings(as.integer(total_measurements))
  ) |>
  dplyr::distinct(molecule_chembl_id, .keep_all = TRUE)

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
  dplyr::filter(
    !is.na(protein), !is.na(biological_process),
    protein != "", biological_process != ""
  ) |>
  dplyr::distinct(protein, biological_process, .keep_all = TRUE)

compound_protein <- compound_protein |>
  dplyr::filter(
    !is.na(molecule_chembl_id), !is.na(queried_protein),
    molecule_chembl_id != "", queried_protein != ""
  ) |>
  dplyr::mutate(
    display_name = dplyr::if_else(
      is.na(display_name) | display_name == "",
      molecule_chembl_id,
      display_name
    ),
    minimum_activity_nM = suppressWarnings(as.numeric(minimum_activity_nM)),
    raw_record_count = suppressWarnings(as.integer(raw_record_count)),
    unique_assay_count = suppressWarnings(as.integer(unique_assay_count)),
    unique_document_count = suppressWarnings(as.integer(unique_document_count)),
    confidence_score = suppressWarnings(as.numeric(confidence_score))
  ) |>
  dplyr::distinct(molecule_chembl_id, queried_protein, .keep_all = TRUE)

# Connectivity metrics.
ppi_degree <- table(c(ppi$source, ppi$target))
ppi_degree <- data.frame(
  protein = names(ppi_degree),
  ppi_degree = as.numeric(ppi_degree),
  stringsAsFactors = FALSE
)
pathway_degree <- protein_pathway |> dplyr::count(protein, name = "pathway_degree")
go_degree <- protein_go |> dplyr::count(protein, name = "go_degree")
compound_degree <- compound_protein |>
  dplyr::count(queried_protein, name = "compound_degree") |>
  dplyr::rename(protein = queried_protein)

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
    total_degree = ppi_degree + pathway_degree + go_degree + compound_degree,
    node_size = 18 + 3.0 * sqrt(total_degree)
  )

pathways <- pathways |>
  dplyr::left_join(
    protein_pathway |> dplyr::count(pathway, name = "protein_count"),
    by = "pathway"
  ) |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 10 + 2.3 * sqrt(protein_count)
  )

go_processes <- go_processes |>
  dplyr::left_join(
    protein_go |> dplyr::count(biological_process, name = "protein_count"),
    by = c("go_term" = "biological_process")
  ) |>
  dplyr::mutate(
    protein_count = dplyr::coalesce(protein_count, 0L),
    node_size = 9 + 2.1 * sqrt(protein_count)
  )

compounds <- compounds |>
  dplyr::left_join(
    compound_protein |> dplyr::count(molecule_chembl_id, name = "target_count"),
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    target_count = dplyr::coalesce(target_count, 0L),
    total_measurements = dplyr::coalesce(total_measurements, 0L),
    node_size = 16 + 3.2 * sqrt(target_count + total_measurements)
  )

# Nodes.
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
      "<strong>GO connections:</strong> ", go_degree,
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
      "<strong>Candidate proteins:</strong> ", dplyr::coalesce(proteins, "Unknown"), "<br>",
      "<strong>Original activity types:</strong> ", dplyr::coalesce(activity_types, "Unknown"), "<br>",
      "<strong>Strongest quantitative record:</strong> ",
      ifelse(is.na(strongest_record_nM), "Unknown", paste0(signif(strongest_record_nM, 4), " nM")),
      "</div>"
    )
  )

nodes <- dplyr::bind_rows(protein_nodes, pathway_nodes, go_nodes, compound_nodes)
if (anyDuplicated(nodes$id) > 0) stop("Duplicate node IDs detected.")

# Edges.
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
    arrows = ""
  )

compound_edges <- compound_protein |>
  dplyr::mutate(
    edge_color = dplyr::case_when(
      relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ "#C43C39",
      relationship == "BINDS_TO" ~ "#8E63B0",
      TRUE ~ "#C77729"
    ),
    edge_width = dplyr::case_when(
      relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ 3.0,
      relationship == "BINDS_TO" ~ 2.4,
      TRUE ~ 1.7
    )
  ) |>
  dplyr::transmute(
    from = paste0("compound::", molecule_chembl_id),
    to = paste0("protein::", queried_protein),
    relationship,
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;font-size:13px;line-height:1.35;'>",
      "<strong>", relationship, "</strong><br>",
      display_name, " - ", queried_protein, "<br>",
      "<strong>Corrected activity types:</strong> ", corrected_activity_types, "<br>",
      "<strong>Original activity types:</strong> ", original_activity_types, "<br>",
      "<strong>Minimum activity:</strong> ",
      ifelse(is.na(minimum_activity_nM), "Unknown", paste0(signif(minimum_activity_nM, 4), " nM")), "<br>",
      "<strong>Independent assays:</strong> ", dplyr::coalesce(unique_assay_count, 0L), "<br>",
      "<strong>Independent documents:</strong> ", dplyr::coalesce(unique_document_count, 0L), "<br>",
      "<strong>Target confidence:</strong> ", dplyr::coalesce(confidence_score, 0), "<br>",
      "<strong>Kd/ED50 duplicate collapsed:</strong> ",
      dplyr::coalesce(duplicate_kd_ed50_collapsed, FALSE),
      "</div>"
    ),
    color = edge_color,
    width = edge_width,
    arrows = ""
  )

edges <- dplyr::bind_rows(ppi_edges, pathway_edges, go_edges, compound_edges)
missing_sources <- setdiff(unique(edges$from), nodes$id)
missing_targets <- setdiff(unique(edges$to), nodes$id)
if (length(missing_sources) > 0 || length(missing_targets) > 0) {
  stop("Some edge endpoints do not have matching nodes.")
}

# Network.
network <- visNetwork::visNetwork(
  nodes = nodes,
  edges = edges,
  width = "100%",
  height = "900px",
  main = "EEF1A Network Explorer: Provenance-Corrected Compound Evidence"
) |>
  visNetwork::visGroups(
    groupname = "Compound", shape = "diamond",
    color = list(
      background = "#E58B3A", border = "#9A4F12",
      highlight = list(background = "#F4AF6A", border = "#7A3D0D"),
      hover = list(background = "#F4AF6A", border = "#7A3D0D")
    ),
    font = list(face = "Calibri", size = 15, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "Protein", shape = "dot",
    color = list(
      background = "#D95F59", border = "#8E2F2A",
      highlight = list(background = "#FF8A80", border = "#641D19"),
      hover = list(background = "#F28E85", border = "#641D19")
    ),
    font = list(face = "Calibri", size = 16, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "Pathway", shape = "dot",
    color = list(
      background = "#4E9F6E", border = "#25633A",
      highlight = list(background = "#78BE92", border = "#174728"),
      hover = list(background = "#78BE92", border = "#174728")
    ),
    font = list(face = "Calibri", size = 14, color = "#202020")
  ) |>
  visNetwork::visGroups(
    groupname = "BiologicalProcess", shape = "dot",
    color = list(
      background = "#5B8FD1", border = "#2B5797",
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
    highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE, algorithm = "all"),
    nodesIdSelection = list(enabled = TRUE, useLabels = TRUE, main = "Find a node"),
    selectedBy = list(variable = "group", main = "Filter by node type", multiple = TRUE)
  ) |>
  visNetwork::visInteraction(
    hover = TRUE, navigationButtons = TRUE, keyboard = TRUE,
    zoomView = TRUE, dragView = TRUE, dragNodes = TRUE,
    hideEdgesOnDrag = TRUE, tooltipDelay = 150
  ) |>
  visNetwork::visPhysics(
    enabled = TRUE,
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(
      gravitationalConstant = -105, centralGravity = 0.006,
      springLength = 260, springConstant = 0.033,
      damping = 0.64, avoidOverlap = 0.8
    ),
    stabilization = list(enabled = TRUE, iterations = 3500, updateInterval = 50, fit = TRUE),
    minVelocity = 0.7, timestep = 0.35, adaptiveTimestep = TRUE
  ) |>
  visNetwork::visLayout(randomSeed = 42, improvedLayout = TRUE) |>
  visNetwork::visEvents(
    stabilizationIterationsDone = "
      function () {
        this.setOptions({physics: {enabled: false}});
        this.fit({animation: {duration: 700, easingFunction: 'easeInOutQuad'}});
      }
    "
  )

# Reliable updated legend above graph.
node_key <- function(color, border, label, diamond = FALSE) {
  htmltools::tags$span(
    style = "display:flex;align-items:center;",
    htmltools::tags$span(
      style = paste0(
        "display:inline-block;width:14px;height:14px;",
        "border-radius:", if (diamond) "2px" else "50%", ";",
        "background:", color, ";border:1px solid ", border, ";",
        "margin-right:8px;", if (diamond) "transform:rotate(45deg);" else ""
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
        "display:inline-block;width:25px;height:0;",
        "border-top:", width, "px solid ", color, ";margin-right:8px;"
      )
    ),
    label
  )
}

legend_panel <- htmltools::tags$div(
  style = paste0(
    "box-sizing:border-box;width:100%;margin:0 0 10px 0;padding:11px 14px;",
    "background:#F8F9FB;border:1px solid #D7DADE;border-radius:5px;",
    "font-family:Calibri,Arial,sans-serif;font-size:13px;color:#222222;"
  ),
  htmltools::tags$div(
    style = "display:flex;flex-wrap:wrap;align-items:center;gap:10px 22px;",
    htmltools::tags$span(style = "font-weight:600;font-size:14px;", "Node types"),
    node_key("#E58B3A", "#9A4F12", "Compound", TRUE),
    node_key("#D95F59", "#8E2F2A", "Protein"),
    node_key("#4E9F6E", "#25633A", "Reactome pathway"),
    node_key("#5B8FD1", "#2B5797", "GO biological process")
  ),
  htmltools::tags$div(style = "border-top:1px solid #D7DADE;margin:10px 0;"),
  htmltools::tags$div(
    style = "display:flex;flex-wrap:wrap;align-items:center;gap:10px 22px;",
    htmltools::tags$span(style = "font-weight:600;font-size:14px;", "Relationship types"),
    edge_key("#C43C39", "Inhibitory-activity evidence"),
    edge_key("#8E63B0", "Binding evidence"),
    edge_key("#C77729", "Other compound activity"),
    edge_key("#9A9A9A", "Protein interaction"),
    edge_key("#40916C", "Pathway association"),
    edge_key("#5B8FD1", "GO association"),
    htmltools::tags$span(
      style = "color:#666666;font-size:12px;margin-left:auto;",
      "Kd/ED50 duplicates collapsed using assay and document provenance"
    )
  )
)

page_style <- htmltools::tags$style(
  htmltools::HTML(
    "
    body {margin:0;padding:14px;background:#FFFFFF;font-family:Calibri,Arial,sans-serif;}
    h1,h2,h3 {margin-top:4px;margin-bottom:10px;color:#222222;font-family:Calibri,Arial,sans-serif;font-weight:600;}
    select,input,button,label {font-family:Calibri,Arial,sans-serif !important;font-size:13px !important;}
    .vis-network {border:1px solid #E0E0E0;border-radius:5px;background:#FFFFFF;}
    "
  )
)

network <- htmlwidgets::prependContent(network, page_style, legend_panel)
htmlwidgets::saveWidget(network, output_file, selfcontained = TRUE)

if (!file.exists(output_file) || is.na(file.info(output_file)$size) || file.info(output_file)$size <= 0) {
  stop("The corrected A13B network HTML file was not created correctly.")
}

cat("\n")
cat("A13B corrected compound network created successfully.\n")
cat("----------------------------------------------------\n")
cat("Compound nodes:             ", nrow(compound_nodes), "\n")
cat("Protein nodes:              ", nrow(protein_nodes), "\n")
cat("Corrected compound edges:   ", nrow(compound_edges), "\n")
cat("Total nodes:                ", nrow(nodes), "\n")
cat("Total edges:                ", nrow(edges), "\n")
cat("Output file:                ", output_file, "\n\n")

output_path <- normalizePath(output_file, winslash = "/", mustWork = TRUE)
utils::browseURL(output_path)
cat("Opened saved network:\n", output_path, "\n")
