#' RESKO Faithful: Step 6 - Interactive HTML network graph (Variant A & B)
#'
#' Builds an explorable HTML network graph via visNetwork, covering both
#' variants in one page (toggle between them), reusing the same node set
#' and per-variant score coloring already validated in the static figure
#' (05_visualize_network_3d.R): eEF1A1/eEF1A2 at center, top-20-per-variant
#' union of candidates as surrounding nodes, node fill = that variant's own
#' composite score (viridis-like scale), Variant B additionally drawing
#' solid edges for interactome-evidenced links (pathway_score > 0) vs.
#' dashed for unsupported.
#'
#' Hover tooltips report drug name, composite score, rank, target count,
#' and (Variant B only) pathway_score -- se_score_undefined is deliberately
#' OMITTED from tooltips rather than labeled as a real side-effect score
#' (see METHODS_DEVIATIONS.md: it is fixed at 0.0 for every candidate and
#' carries no information).

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(visNetwork)
  library(htmlwidgets)
})

PROCESSED_DIR <- file.path("data", "processed")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: INTERACTIVE HTML NETWORK GRAPH (VARIANT A & B)\n")
cat(strrep("=", 70), "\n")

variantB <- read_csv(file.path(PROCESSED_DIR, "resko_variantB_interactome_candidates.csv"), show_col_types = FALSE)
N_TOP <- 20
topA <- variantB %>% slice_max(composite_score_variantA, n = N_TOP)
topB <- variantB %>% slice_max(composite_score_variantB, n = N_TOP)
union_names <- sort(union(topA$drug_name, topB$drug_name))
cat(sprintf("  Top %d per variant, union = %d unique candidates\n", N_TOP, length(union_names)))

union_df <- variantB %>% filter(drug_name %in% union_names)

# Shared ring layout, matching 05_visualize_network_3d.R exactly, in
# visNetwork's coordinate convention (y flipped: vis draws +y downward).
n <- length(union_names)
angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
RADIUS <- 400
layout_df <- tibble(
  drug_name = union_names,
  x = cos(angles) * RADIUS,
  y = -sin(angles) * RADIUS
)

viridis_hex <- function(t) {
  # Approximate viridis palette without importing scales/viridisLite twice.
  grDevices::colorRampPalette(c("#440154", "#3b528b", "#21918c", "#5ec962", "#fde725"))(101)[
    pmin(100, pmax(0, round(t * 100))) + 1
  ]
}

build_variant_network <- function(score_col, show_interactome_support, group_label) {
  d <- union_df %>%
    select(drug_id, drug_name, score = all_of(score_col), pathway_score,
           n_interactome_targets_hit, max_stage, moa) %>%
    left_join(layout_df, by = "drug_name")
  vmax <- max(max(d$score, na.rm = TRUE), 1e-6)
  d$score_norm <- d$score / vmax
  d$color <- viridis_hex(d$score_norm)

  tooltip_extra <- if (show_interactome_support) {
    sprintf("<br>pathway_score: %.3f", d$pathway_score)
  } else {
    "<br>pathway_score: n/a (Variant A has no interactome term)"
  }

  nodes_candidates <- tibble(
    id = d$drug_id,
    label = d$drug_name,
    x = d$x, y = d$y,
    color = d$color,
    shape = "dot",
    size = 16,
    title = paste0(
      "<b>", d$drug_name, "</b>",
      "<br>", group_label, " composite score: ", sprintf("%.4f", d$score),
      "<br>interactome targets hit: ", d$n_interactome_targets_hit,
      "<br>development stage: ", d$max_stage,
      "<br>MoA: ", d$moa,
      tooltip_extra
    ),
    group = "candidate"
  )
  nodes_target <- tibble(
    id = "EEF1A_TARGET", label = "eEF1A1/\neEF1A2", x = 0, y = 0,
    color = "#d62728", shape = "dot", size = 34,
    title = "eEF1A1/eEF1A2 -- the shared, non-canonical-function target this screen is built around",
    group = "target"
  )
  nodes <- bind_rows(nodes_target, nodes_candidates)

  if (show_interactome_support) {
    edges <- tibble(
      from = "EEF1A_TARGET", to = d$drug_id,
      dashes = !(d$pathway_score > 0),
      width = ifelse(d$pathway_score > 0, 2.5, 1),
      color = ifelse(d$pathway_score > 0, "#444444", "#aaaaaa"),
      title = ifelse(d$pathway_score > 0,
                      "interactome-evidenced direct link (pathway_score > 0)",
                      "no direct interactome evidence for this candidate's targets")
    )
  } else {
    edges <- tibble(
      from = "EEF1A_TARGET", to = d$drug_id,
      dashes = FALSE, width = 1.5, color = "#999999",
      title = "edge weight = composite_score_variantA (se_score fixed at 0/undefined)"
    )
  }
  list(nodes = nodes, edges = edges)
}

cat("\n[1/3] Building Variant A network (indication-breadth + on-target-promiscuity)\n")
netA <- build_variant_network("composite_score_variantA", FALSE, "Variant A")

cat("[2/3] Building Variant B network (+ interactome direct-coupling term)\n")
netB <- build_variant_network("composite_score_variantB", TRUE, "Variant B")

cat("[3/3] Rendering visNetwork widgets and writing standalone HTML\n")

make_widget <- function(net, title_text) {
  visNetwork(net$nodes, net$edges, main = title_text, width = "100%", height = "650px") %>%
    visNodes(font = list(size = 13)) %>%
    visEdges(smooth = FALSE) %>%
    visOptions(highlightNearest = TRUE, nodesIdSelection = TRUE) %>%
    visInteraction(hover = TRUE, dragNodes = FALSE, zoomView = TRUE) %>%
    visPhysics(enabled = FALSE) %>%
    visLayout(randomSeed = 42)
}

widgetA <- make_widget(netA, "Variant A: RESKO minus side-effects (indication-breadth + on-target-promiscuity)")
widgetB <- make_widget(netB, "Variant B: RESKO minus side-effects, interactome-supplemented")

saveWidget(widgetA, file.path(normalizePath("."), "resko_variantA_network_interactive.html"), selfcontained = TRUE)
saveWidget(widgetB, file.path(normalizePath("."), "resko_variantB_network_interactive.html"), selfcontained = TRUE)

cat("\n  Saved resko_variantA_network_interactive.html\n")
cat("  Saved resko_variantB_network_interactive.html\n")
cat("\nInteractive HTML network graphs complete\n")
cat("  NOTE: se_score_undefined is deliberately omitted from tooltips (fixed at 0.0,\n")
cat("  carries no information) -- see METHODS_DEVIATIONS.md.\n")
