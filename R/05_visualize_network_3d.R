#' RESKO Faithful: Step 5 - Render comparison network graphs for Variant A
#' and Variant B.
#'
#' Both figures share the SAME node layout (computed once from the union of
#' top-20 candidates across both variants, seeds, and eEF1A) so a reader can
#' visually track any candidate's position between the two panels. What
#' differs between the figures is what is real and different between the
#' variants:
#'   - Variant A: edge weight = composite_score_variantA (indication-breadth +
#'     on-target-promiscuity only; se_score fixed at 0/undefined)
#'   - Variant B: edge weight = composite_score_variantB (adds pathway_score,
#'     the interactome direct-coupling-strength term); edges from candidates
#'     whose pathway_score > 0 (i.e., an interactome-supported link exists)
#'     are drawn as solid, all others as dashed -- this is the one visual
#'     encoding with no equivalent in Variant A, because pathway_score does
#'     not exist there.
#'
#' Layout is a ring (eEF1A at center, candidates evenly spaced on a circle)
#' -- identical construction to the Python original, chosen there because a
#' spring layout produced two overlapping label positions for this
#' star-topology graph; a ring guarantees zero label collision.

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(ggplot2)
  library(grid)
})

PROCESSED_DIR <- file.path("data", "processed")

cat(strrep("=", 70), "\n")
cat("RESKO FAITHFUL: RENDER VARIANT A vs VARIANT B NETWORK GRAPHS\n")
cat(strrep("=", 70), "\n")

variantB <- read_csv(file.path(PROCESSED_DIR, "resko_variantB_interactome_candidates.csv"), show_col_types = FALSE)
N_TOP <- 20
topA <- variantB %>% slice_max(composite_score_variantA, n = N_TOP)
topB <- variantB %>% slice_max(composite_score_variantB, n = N_TOP)
union_names <- sort(union(topA$drug_name, topB$drug_name))
cat(sprintf("  Top %d per variant, union = %d unique candidates\n", N_TOP, length(union_names)))

union_df <- variantB %>% filter(drug_name %in% union_names)

# Build ONE shared ring layout: eEF1A at center, candidates evenly spaced.
n <- length(union_names)
angles <- seq(0, 2 * pi, length.out = n + 1)[1:n]
layout_df <- tibble(
  drug_name = union_names,
  x = cos(angles) * 1.6,
  y = sin(angles) * 1.6,
  theta = angles
)

build_panel_data <- function(score_col, show_interactome_support) {
  d <- union_df %>%
    select(drug_name, score = all_of(score_col), pathway_score) %>%
    left_join(layout_df, by = "drug_name") %>%
    mutate(
      supported = if (show_interactome_support) (pathway_score > 0) else NA,
      theta_deg = theta * 180 / pi,
      ha = ifelse(theta_deg > -90 & theta_deg <= 90, "left", "right"),
      rot = ifelse(ha == "right", theta_deg + 180, theta_deg),
      lx = x * 1.14, ly = y * 1.14
    )
  d
}

render_panel <- function(score_col, title, show_interactome_support, panel_letter_txt) {
  d <- build_panel_data(score_col, show_interactome_support)
  vmax <- max(max(d$score, na.rm = TRUE), 1e-6)

  p <- ggplot() +
    { if (show_interactome_support) {
        list(
          geom_segment(data = d %>% filter(supported), aes(x = 0, y = 0, xend = x, yend = y),
                       color = "#444444", linewidth = 0.5, lineend = "round"),
          geom_segment(data = d %>% filter(!supported), aes(x = 0, y = 0, xend = x, yend = y),
                       color = "#aaaaaa", linewidth = 0.35, linetype = "dashed")
        )
      } else {
        geom_segment(data = d, aes(x = 0, y = 0, xend = x, yend = y), color = "#999999", linewidth = 0.4)
      } } +
    geom_point(aes(x = 0, y = 0), size = 9, shape = 21, fill = "#d62728", color = "black", stroke = 1.0) +
    annotate("text", x = 0, y = 0, label = "eEF1A1/\neEF1A2", size = 2.4, fontface = "bold") +
    geom_point(data = d, aes(x = x, y = y, fill = score), size = 3.3, shape = 21, color = "black", stroke = 0.4) +
    scale_fill_viridis_c(name = gsub("composite_score_variant", "composite score, variant ", score_col),
                          limits = c(0, vmax)) +
    ggtitle(title) +
    coord_fixed(xlim = c(-2.3, 2.3), ylim = c(-2.3, 2.3), clip = "off") +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(size = 9, hjust = 0),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 6),
      plot.margin = margin(20, 30, 10, 30)
    )

  for (i in seq_len(nrow(d))) {
    r <- d[i, ]
    p <- p + annotate("text", x = r$lx, y = r$ly, label = r$drug_name, size = 2.1,
                       angle = r$rot, hjust = ifelse(r$ha == "left", 0, 1), vjust = 0.5)
  }
  p
}

pA <- render_panel("composite_score_variantA",
                    "Variant A: indication-breadth + on-target-promiscuity only\n(se_score undefined -- no real side-effect data available for this seed set)",
                    FALSE, "a")
pB <- render_panel("composite_score_variantB",
                    "Variant B: Variant A + interactome direct-coupling term\n(solid edge = candidate target has a direct, evidenced eEF1A interactome link)",
                    TRUE, "b")

if (!requireNamespace("patchwork", quietly = TRUE)) install.packages("patchwork", repos = "https://cloud.r-project.org")
library(patchwork)

fig <- (pA | pB) +
  plot_annotation(
    title = "RESKO-Faithful candidate ranking: real-data variants without side-effect similarity",
    tag_levels = "a"
  ) &
  theme(plot.tag = element_text(size = 13, face = "bold"))

out_png <- "resko_variantA_vs_variantB_network_R.png"
ggsave(out_png, plot = fig, width = 15, height = 7.5, dpi = 300, bg = "white")
cat(sprintf("\n  Saved comparison figure to %s\n", out_png))
cat("\nVisualization complete\n")
