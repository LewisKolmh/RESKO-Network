# ============================================================
# A15B RANKED INTERACTIVE eEF1A NETWORK
#
# Uses the provenance-corrected A13B network and A15 rankings.
# Adds rank-sensitive compound styling, tooltips, two ranking
# tables below the graph, and a dedicated ranking-visual key.
#
# Output:
#   results/eef1a_network_A15_ranked.html
# ============================================================

required_packages <- c("readr", "dplyr", "visNetwork", "htmlwidgets", "htmltools")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) stop(paste("Missing packages:", paste(missing_packages, collapse = ", ")))

results_dir <- file.path(normalizePath(getwd(), winslash = "/", mustWork = TRUE), "results")
files <- list(
  proteins = file.path(results_dir, "nodes_proteins.csv"),
  pathways = file.path(results_dir, "nodes_pathways.csv"),
  go = file.path(results_dir, "nodes_biological_process.csv"),
  compounds = file.path(results_dir, "nodes_drugs.csv"),
  ppi = file.path(results_dir, "edges_interacts_with.csv"),
  protein_pathway = file.path(results_dir, "edges_protein_pathway.csv"),
  protein_go = file.path(results_dir, "edges_protein_go.csv"),
  compound_edges = file.path(results_dir, "edges_compound_protein_corrected.csv"),
  scores = file.path(results_dir, "A15_candidate_score_components.csv"),
  direct = file.path(results_dir, "A15_direct_eef1a1_ranking.csv"),
  network = file.path(results_dir, "A15_translation_network_ranking.csv")
)
output_file <- file.path(results_dir, "eef1a_network_A15_ranked.html")
missing_files <- unlist(files)[!file.exists(unlist(files))]
if (length(missing_files) > 0) stop(paste0("Missing inputs:\n", paste(missing_files, collapse = "\n")))

read_data <- function(path) readr::read_csv(path, show_col_types = FALSE, progress = FALSE)
proteins <- read_data(files$proteins)
pathways <- read_data(files$pathways)
go_processes <- read_data(files$go)
compounds <- read_data(files$compounds)
ppi <- read_data(files$ppi)
protein_pathway <- read_data(files$protein_pathway)
protein_go <- read_data(files$protein_go)
compound_edges_data <- read_data(files$compound_edges)
scores <- read_data(files$scores)
direct_ranking <- read_data(files$direct)
network_ranking <- read_data(files$network)

# Merge ranks and scores into compound nodes.
compound_rank_data <- compounds |>
  dplyr::left_join(scores, by = "molecule_chembl_id") |>
  dplyr::left_join(
    direct_ranking |> dplyr::select(molecule_chembl_id, direct_rank),
    by = "molecule_chembl_id"
  ) |>
  dplyr::left_join(
    network_ranking |> dplyr::select(molecule_chembl_id, network_rank),
    by = "molecule_chembl_id"
  ) |>
  dplyr::mutate(
    name = dplyr::coalesce(name, molecule_chembl_id),
    overall_priority_score = dplyr::coalesce(overall_priority_score, 0),
    compound_size = 18 + 1.6 * overall_priority_score,
    is_top_direct = direct_rank == min(direct_rank, na.rm = TRUE),
    is_top_network = network_rank == min(network_rank, na.rm = TRUE),
    compound_border_width = dplyr::case_when(
      is_top_direct & is_top_network ~ 7,
      is_top_direct ~ 6,
      is_top_network ~ 4,
      TRUE ~ 2
    ),
    compound_shadow_size = dplyr::case_when(
      is_top_network ~ 18,
      TRUE ~ 5
    ),
    compound_shadow_color = dplyr::case_when(
      is_top_network ~ "rgba(70,35,120,0.55)",
      TRUE ~ "rgba(0,0,0,0.16)"
    ),
    rank_label = paste0("D", direct_rank, " / N", network_rank, "  ", name)
  )

# Connectivity calculations for biological nodes.
ppi <- ppi |> dplyr::mutate(score = dplyr::coalesce(suppressWarnings(as.numeric(score)), 0))
ppi_degree <- table(c(ppi$source, ppi$target))
ppi_degree <- data.frame(protein = names(ppi_degree), ppi_degree = as.numeric(ppi_degree))
pathway_degree <- protein_pathway |> dplyr::count(protein, name = "pathway_degree")
go_degree <- protein_go |> dplyr::count(protein, name = "go_degree")
compound_degree <- compound_edges_data |> dplyr::count(queried_protein, name = "compound_degree") |> dplyr::rename(protein = queried_protein)

protein_metrics <- proteins |>
  dplyr::left_join(ppi_degree, by = "protein") |>
  dplyr::left_join(pathway_degree, by = "protein") |>
  dplyr::left_join(go_degree, by = "protein") |>
  dplyr::left_join(compound_degree, by = "protein") |>
  dplyr::mutate(
    dplyr::across(c(ppi_degree, pathway_degree, go_degree, compound_degree), ~ dplyr::coalesce(.x, 0)),
    node_size = 18 + 3 * sqrt(ppi_degree + pathway_degree + go_degree + compound_degree)
  )

pathways <- pathways |>
  dplyr::left_join(protein_pathway |> dplyr::count(pathway, name = "protein_count"), by = "pathway") |>
  dplyr::mutate(protein_count = dplyr::coalesce(protein_count, 0L), node_size = 10 + 2.3 * sqrt(protein_count))

go_processes <- go_processes |>
  dplyr::left_join(protein_go |> dplyr::count(biological_process, name = "protein_count"), by = c("go_term" = "biological_process")) |>
  dplyr::mutate(protein_count = dplyr::coalesce(protein_count, 0L), node_size = 9 + 2.1 * sqrt(protein_count))

protein_nodes <- protein_metrics |>
  dplyr::transmute(
    id = paste0("protein::", protein), label = protein, group = "Protein", shape = "dot", size = node_size,
    title = paste0("<b>Protein:</b> ", protein, "<br><b>Compound links:</b> ", compound_degree)
  )
pathway_nodes <- pathways |>
  dplyr::transmute(
    id = paste0("pathway::", pathway), label = pathway, group = "Pathway", shape = "dot", size = node_size,
    title = paste0("<b>Reactome pathway:</b> ", pathway, "<br><b>Adjusted p:</b> ", signif(p_adjust, 3))
  )
go_nodes <- go_processes |>
  dplyr::transmute(
    id = paste0("go::", go_term), label = go_term, group = "BiologicalProcess", shape = "dot", size = node_size,
    title = paste0("<b>GO biological process:</b> ", go_term)
  )
compound_nodes <- compound_rank_data |>
  dplyr::transmute(
    id = paste0("compound::", molecule_chembl_id),
    label = rank_label,
    group = "Compound",
    shape = "diamond",
    size = compound_size,
    borderWidth = compound_border_width,
    shadow = Map(function(enabled, color, size) list(enabled = enabled, color = color, size = size, x = 0, y = 0), TRUE, compound_shadow_color, compound_shadow_size),
    title = paste0(
      "<div style='font-family:Calibri,Arial,sans-serif;line-height:1.4;'>",
      "<strong>Compound:</strong> ", name, "<br>",
      "<strong>ChEMBL ID:</strong> ", molecule_chembl_id, "<br><br>",
      "<strong>Direct EEF1A1 rank:</strong> ", direct_rank, "<br>",
      "<strong>Translation-network rank:</strong> ", network_rank, "<br>",
      "<strong>Direct score:</strong> ", round(direct_total_score, 2), " / 16<br>",
      "<strong>Network score:</strong> ", round(network_total_score, 2), "<br>",
      "<strong>Priority category:</strong> ", priority_category, "<br>",
      "<strong>Direct primary target:</strong> ", direct_primary_target, "<br>",
      "<strong>Direct relationship:</strong> ", direct_primary_relationship, "<br>",
      "<strong>Direct minimum activity:</strong> ", signif(direct_minimum_activity_nM, 4), " nM<br>",
      "<strong>Strongest network target:</strong> ", strongest_network_target, "<br>",
      "<strong>Network target count:</strong> ", network_target_count,
      "</div>"
    )
  )

nodes <- dplyr::bind_rows(protein_nodes, pathway_nodes, go_nodes, compound_nodes)

ppi_edges <- ppi |>
  dplyr::transmute(from = paste0("protein::", source), to = paste0("protein::", target), color = "#9A9A9A", width = 0.9 + 2.7 * score, arrows = "")
pathway_edges <- protein_pathway |>
  dplyr::transmute(from = paste0("protein::", protein), to = paste0("pathway::", pathway), color = "#40916C", width = 1.15, arrows = "")
go_edges <- protein_go |>
  dplyr::transmute(from = paste0("protein::", protein), to = paste0("go::", biological_process), color = "#5B8FD1", width = 1.05, arrows = "")
compound_edges <- compound_edges_data |>
  dplyr::mutate(
    edge_color = dplyr::case_when(
      relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ "#C43C39",
      relationship == "BINDS_TO" ~ "#8E63B0",
      TRUE ~ "#C77729"
    ),
    edge_width = dplyr::case_when(
      relationship == "HAS_INHIBITORY_ACTIVITY_AGAINST" ~ 3,
      relationship == "BINDS_TO" ~ 2.4,
      TRUE ~ 1.7
    )
  ) |>
  dplyr::transmute(
    from = paste0("compound::", molecule_chembl_id),
    to = paste0("protein::", queried_protein),
    color = edge_color, width = edge_width, arrows = "",
    title = paste0("<b>", relationship, "</b><br>", display_name, " - ", queried_protein, "<br>Minimum activity: ", signif(minimum_activity_nM, 4), " nM")
  )
edges <- dplyr::bind_rows(ppi_edges, pathway_edges, go_edges, compound_edges)

network <- visNetwork::visNetwork(nodes, edges, width = "100%", height = "900px", main = "Ranked eEF1A Compound Network") |>
  visNetwork::visGroups(groupname = "Compound", shape = "diamond", color = list(background = "#E58B3A", border = "#713900"), font = list(face = "Calibri", size = 14)) |>
  visNetwork::visGroups(groupname = "Protein", shape = "dot", color = list(background = "#D95F59", border = "#8E2F2A"), font = list(face = "Calibri", size = 16)) |>
  visNetwork::visGroups(groupname = "Pathway", shape = "dot", color = list(background = "#4E9F6E", border = "#25633A"), font = list(face = "Calibri", size = 14)) |>
  visNetwork::visGroups(groupname = "BiologicalProcess", shape = "dot", color = list(background = "#5B8FD1", border = "#2B5797"), font = list(face = "Calibri", size = 13)) |>
  visNetwork::visNodes(font = list(face = "Calibri", strokeWidth = 3, strokeColor = "#FFFFFF")) |>
  visNetwork::visEdges(arrows = "", smooth = FALSE, selectionWidth = 3, hoverWidth = 2) |>
  visNetwork::visOptions(
    highlightNearest = list(enabled = TRUE, degree = 1, hover = FALSE, algorithm = "all"),
    nodesIdSelection = list(enabled = TRUE, useLabels = TRUE, main = "Find a node"),
    selectedBy = list(variable = "group", main = "Filter by node type", multiple = TRUE)
  ) |>
  visNetwork::visInteraction(hover = TRUE, navigationButtons = TRUE, keyboard = TRUE, hideEdgesOnDrag = TRUE) |>
  visNetwork::visPhysics(
    solver = "forceAtlas2Based",
    forceAtlas2Based = list(gravitationalConstant = -110, centralGravity = 0.006, springLength = 265, springConstant = 0.033, damping = 0.65, avoidOverlap = 0.8),
    stabilization = list(enabled = TRUE, iterations = 3500, fit = TRUE)
  ) |>
  visNetwork::visLayout(randomSeed = 42, improvedLayout = TRUE) |>
  visNetwork::visEvents(stabilizationIterationsDone = "function(){this.setOptions({physics:{enabled:false}});this.fit();}")

# Rendering helpers for table panels.
html_table <- function(data, columns, labels) {
  shown <- data |> dplyr::select(dplyr::all_of(columns))
  names(shown) <- labels
  htmltools::tags$table(
    class = "ranking-table",
    htmltools::tags$thead(htmltools::tags$tr(lapply(names(shown), htmltools::tags$th))),
    htmltools::tags$tbody(lapply(seq_len(nrow(shown)), function(i) {
      htmltools::tags$tr(lapply(shown[i, , drop = FALSE], function(x) htmltools::tags$td(as.character(x))))
    }))
  )
}

rank_visual_key <- htmltools::tags$div(
  class = "rank-key",
  htmltools::tags$h3("How ranking is represented in the graph"),
  htmltools::tags$div(class = "rank-key-grid",
    htmltools::tags$div(class = "symbol diamond orange"), htmltools::tags$div("Orange diamond"), htmltools::tags$div("Compound"),
    htmltools::tags$div(class = "symbol diamond orange large"), htmltools::tags$div("Larger diamond"), htmltools::tags$div("Higher overall priority"),
    htmltools::tags$div(class = "symbol diamond orange thick"), htmltools::tags$div("Thick dark border"), htmltools::tags$div("Top-ranked direct candidate"),
    htmltools::tags$div(class = "symbol diamond orange halo"), htmltools::tags$div("Dark halo"), htmltools::tags$div("Top-ranked translation-network candidate"),
    htmltools::tags$div(class = "tooltip-symbol", "i"), htmltools::tags$div("Tooltip"), htmltools::tags$div("Full ranking evidence"),
    htmltools::tags$div(class = "label-symbol", "D1 / N2"), htmltools::tags$div("Label prefix"), htmltools::tags$div("Direct and network rank"),
    htmltools::tags$div(class = "opacity-symbol"), htmltools::tags$div("Reduced opacity"), htmltools::tags$div("Unrelated nodes after selection")
  )
)

leaderboards <- htmltools::tags$div(
  class = "ranking-section",
  htmltools::tags$h2("Candidate rankings"),
  htmltools::tags$div(class = "leaderboard-grid",
    htmltools::tags$div(class = "leaderboard-card",
      htmltools::tags$h3("Direct EEF1A1 ranking"),
      html_table(
        direct_ranking,
        c("direct_rank", "display_name", "direct_total_score", "direct_primary_target", "direct_primary_relationship", "direct_minimum_activity_nM"),
        c("Rank", "Compound", "Score", "Primary target", "Evidence", "Minimum nM")
      )
    ),
    htmltools::tags$div(class = "leaderboard-card",
      htmltools::tags$h3("Translation-network ranking"),
      html_table(
        network_ranking,
        c("network_rank", "display_name", "network_total_score", "strongest_network_target", "network_target_count", "strongest_network_activity_nM"),
        c("Rank", "Compound", "Score", "Strongest target", "Targets", "Minimum nM")
      )
    )
  ),
  rank_visual_key,
  htmltools::tags$p(class = "report-link", "Detailed score decomposition: results/A15_candidate_rankings.html")
)

page_style <- htmltools::tags$style(htmltools::HTML("\nbody{margin:0;padding:14px;background:#fff;font-family:Calibri,Arial,sans-serif;color:#222;}\nh1,h2,h3{font-weight:600;} select,input,button,label{font-family:Calibri,Arial,sans-serif!important;}\n.vis-network{border:1px solid #e0e0e0;border-radius:5px;background:#fff;}\n.ranking-section{margin-top:18px;}\n.leaderboard-grid{display:grid;grid-template-columns:repeat(auto-fit,minmax(520px,1fr));gap:18px;}\n.leaderboard-card,.rank-key{border:1px solid #d7dade;border-radius:6px;background:#f8f9fb;padding:14px;}\n.ranking-table{border-collapse:collapse;width:100%;font-size:12px;background:#fff;}\n.ranking-table th{background:#e9eef3;text-align:left;padding:7px;border:1px solid #ccd3da;}\n.ranking-table td{padding:7px;border:1px solid #d9dee3;}\n.rank-key{margin-top:18px;}\n.rank-key-grid{display:grid;grid-template-columns:70px 150px 1fr;gap:10px 12px;align-items:center;}\n.symbol{width:18px;height:18px;margin-left:14px;}\n.diamond{transform:rotate(45deg);border-radius:2px;} .orange{background:#E58B3A;border:2px solid #713900;}\n.large{width:28px;height:28px;margin-left:9px;} .thick{border-width:5px;}\n.halo{box-shadow:0 0 0 5px rgba(70,35,120,.45),0 0 14px 7px rgba(70,35,120,.4);}\n.tooltip-symbol{width:24px;height:24px;border-radius:50%;background:#eef2f6;border:1px solid #8996a3;text-align:center;line-height:24px;font-weight:bold;margin-left:10px;}\n.label-symbol{font-weight:600;color:#713900;} .opacity-symbol{width:30px;height:5px;background:rgba(154,154,154,.2);margin-left:6px;}\n.report-link{font-size:12px;color:#555;}\n"))

network <- htmlwidgets::appendContent(network, leaderboards)
network <- htmlwidgets::prependContent(network, page_style)
htmlwidgets::saveWidget(network, output_file, selfcontained = TRUE)

if (!file.exists(output_file) || file.info(output_file)$size <= 0) stop("Ranked network HTML was not created.")
cat("\nA15B ranked network created successfully.\n")
cat("Output: ", output_file, "\n", sep = "")
utils::browseURL(normalizePath(output_file, winslash = "/", mustWork = TRUE))
