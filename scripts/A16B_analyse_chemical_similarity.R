#!/usr/bin/env Rscript

# ============================================================
# RESKO A16B: Chemical similarity analysis in R
# ============================================================
# Run from the RESKO project root while resko-a16 is active:
#   conda activate resko-a16
#   Rscript scripts/A16B_analyse_chemical_similarity.R
#
# Scientific implementation:
#   R controls validation, analysis, clustering, plotting, reporting,
#   and output verification. RDKit is accessed through reticulate so
#   the requested Morgan radius-2, 2048-bit fingerprints and chemical
#   depictions remain consistent with the original A16 specification.
# ============================================================

options(stringsAsFactors = FALSE, warn = 1)

required_r_packages <- c("readr", "dplyr", "tibble", "ggplot2", "reticulate", "base64enc")
missing_r_packages <- required_r_packages[
  !vapply(required_r_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))
]
if (length(missing_r_packages) > 0L) {
  stop(
    "Missing required R package(s): ", paste(missing_r_packages, collapse = ", "),
    ". Install from the terminal with Rscript -e 'install.packages(c(",
    paste(sprintf("\"%s\"", missing_r_packages), collapse = ", "),
    "), repos=\"https://cloud.r-project.org\")'.",
    call. = FALSE
  )
}

project_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
if (!dir.exists(file.path(project_root, "scripts")) || !dir.exists(file.path(project_root, "results"))) {
  stop("Run this script from the RESKO project root containing scripts/ and results/.", call. = FALSE)
}

results_dir <- file.path(project_root, "results")
input_file <- file.path(results_dir, "A16_structures_combined.csv")

outputs <- c(
  pairwise = file.path(results_dir, "A16_pairwise_tanimoto.csv"),
  candidate_reference = file.path(results_dir, "A16_candidate_reference_similarity.csv"),
  clusters = file.path(results_dir, "A16_structural_clusters.csv"),
  heatmap = file.path(results_dir, "A16_similarity_heatmap.png"),
  chemical_space = file.path(results_dir, "A16_chemical_space.png"),
  structure_grid = file.path(results_dir, "A16_structure_grid.png"),
  summary = file.path(results_dir, "A16_similarity_summary.csv"),
  report = file.path(results_dir, "A16_chemical_similarity_report.html")
)

if (!file.exists(input_file) || is.na(file.info(input_file)$size) || file.info(input_file)$size <= 0L) {
  stop("Required input is missing or empty: ", input_file, call. = FALSE)
}

clean_text <- function(x) {
  output <- as.character(x)
  output[is.na(output) | trimws(output) == ""] <- NA_character_
  output
}

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;", x, fixed = TRUE)
  x <- gsub(">", "&gt;", x, fixed = TRUE)
  x <- gsub('"', "&quot;", x, fixed = TRUE)
  x
}

backup_existing <- function(path) {
  if (!file.exists(path)) return(invisible(NULL))
  stamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  extension <- tools::file_ext(path)
  stem <- tools::file_path_sans_ext(path)
  backup_path <- paste0(stem, "_previous_", stamp, ".", extension)
  if (!file.copy(path, backup_path, overwrite = FALSE)) {
    stop("Could not preserve previous output: ", path, call. = FALSE)
  }
  if (!file.exists(backup_path) || file.info(backup_path)$size <= 0L) {
    stop("Backup verification failed: ", backup_path, call. = FALSE)
  }
  invisible(backup_path)
}

safe_write_csv <- function(data, path) {
  temp_path <- paste0(path, ".tmp")
  if (file.exists(temp_path)) unlink(temp_path)
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  list_columns <- names(data)[vapply(data, is.list, logical(1))]
  if (length(list_columns) > 0L) {
    stop("Refusing to write list-column(s): ", paste(list_columns, collapse = ", "), call. = FALSE)
  }
  readr::write_csv(data, temp_path, na = "")
  if (!file.exists(temp_path) || file.info(temp_path)$size <= 0L) {
    stop("Failed to create temporary output: ", temp_path, call. = FALSE)
  }
  backup_existing(path)
  if (!file.rename(temp_path, path)) stop("Could not move output into place: ", path, call. = FALSE)
  if (!file.exists(path) || file.info(path)$size <= 0L) stop("Output verification failed: ", path, call. = FALSE)
  invisible(path)
}

safe_ggsave <- function(plot, path, width, height, dpi = 220) {
  temp_path <- tempfile(pattern = "A16_", tmpdir = results_dir, fileext = ".png")
  ggplot2::ggsave(temp_path, plot = plot, width = width, height = height, dpi = dpi, units = "in", bg = "white")
  if (!file.exists(temp_path) || file.info(temp_path)$size <= 0L) {
    stop("Failed to create plot: ", path, call. = FALSE)
  }
  backup_existing(path)
  if (!file.rename(temp_path, path)) stop("Could not move plot into place: ", path, call. = FALSE)
  if (!file.exists(path) || file.info(path)$size <= 0L) stop("Plot verification failed: ", path, call. = FALSE)
  invisible(path)
}

table_to_html <- function(data, digits = 3L) {
  data <- as.data.frame(data, stringsAsFactors = FALSE)
  numeric_columns <- vapply(data, is.numeric, logical(1))
  data[numeric_columns] <- lapply(data[numeric_columns], function(x) {
    ifelse(is.na(x), "", format(round(x, digits), nsmall = digits, trim = TRUE))
  })
  data[!numeric_columns] <- lapply(data[!numeric_columns], function(x) {
    ifelse(is.na(x), "", html_escape(x))
  })
  header <- paste0("<tr>", paste0("<th>", html_escape(names(data)), "</th>", collapse = ""), "</tr>")
  body <- paste(vapply(seq_len(nrow(data)), function(i) {
    paste0("<tr>", paste0("<td>", unlist(data[i, , drop = FALSE]), "</td>", collapse = ""), "</tr>")
  }, character(1)), collapse = "\n")
  paste0("<table class='data'><thead>", header, "</thead><tbody>", body, "</tbody></table>")
}

# -----------------------------
# Read and validate input
# -----------------------------
compounds <- tryCatch(
  readr::read_csv(input_file, show_col_types = FALSE, progress = FALSE),
  error = function(e) stop("Failed to read A16 input: ", conditionMessage(e), call. = FALSE)
)

required_columns <- c(
  "compound_id", "compound_name", "compound_class", "structure_source",
  "analysis_smiles", "inchi_key", "structure_available"
)
missing_columns <- setdiff(required_columns, names(compounds))
if (length(missing_columns) > 0L) {
  stop("Input is missing required column(s): ", paste(missing_columns, collapse = ", "), call. = FALSE)
}
if (nrow(compounds) != 10L) stop("Expected 10 compounds, found ", nrow(compounds), ".", call. = FALSE)
if (anyDuplicated(compounds$compound_id) > 0L) stop("Duplicate compound_id values detected.", call. = FALSE)
if (sum(compounds$compound_class == "network_candidate") != 6L) stop("Expected 6 network candidates.", call. = FALSE)
if (sum(compounds$compound_class == "reference_ligand") != 4L) stop("Expected 4 reference ligands.", call. = FALSE)
if (any(!compounds$compound_class %in% c("network_candidate", "reference_ligand"))) {
  stop("Unexpected compound_class values detected.", call. = FALSE)
}
compounds$analysis_smiles <- clean_text(compounds$analysis_smiles)
if (any(is.na(compounds$analysis_smiles))) stop("One or more analysis_smiles values are missing.", call. = FALSE)

structure_available <- tolower(as.character(compounds$structure_available)) %in% c("true", "t", "1", "yes")
if (any(!structure_available)) stop("One or more compounds are marked as lacking a structure.", call. = FALSE)

# -----------------------------
# Bind reticulate to active Conda Python
# -----------------------------
active_python <- Sys.which("python")
if (active_python == "") stop("Python was not found. Activate resko-a16 before running A16B.", call. = FALSE)
Sys.setenv(RETICULATE_PYTHON = active_python)

if (!reticulate::py_module_available("rdkit")) {
  stop(
    "RDKit is not available to reticulate through: ", active_python,
    ". Activate the resko-a16 Conda environment and rerun.",
    call. = FALSE
  )
}

Chem <- reticulate::import("rdkit.Chem", convert = FALSE)
DataStructs <- reticulate::import("rdkit.DataStructs", convert = FALSE)
rdFingerprintGenerator <- reticulate::import("rdkit.Chem.rdFingerprintGenerator", convert = FALSE)
Descriptors <- reticulate::import("rdkit.Chem.Descriptors", convert = FALSE)
Crippen <- reticulate::import("rdkit.Chem.Crippen", convert = FALSE)
Draw <- reticulate::import("rdkit.Chem.Draw", convert = FALSE)

# -----------------------------
# Parse structures, calculate fingerprints and descriptors
# -----------------------------
molecules <- lapply(compounds$analysis_smiles, function(smiles) Chem$MolFromSmiles(smiles))
valid_structure <- !vapply(molecules, reticulate::py_is_null_xptr, logical(1))
if (any(!valid_structure)) {
  invalid_names <- compounds$compound_name[!valid_structure]
  stop("RDKit could not parse: ", paste(invalid_names, collapse = "; "), call. = FALSE)
}

# ECFP4-equivalent Morgan fingerprint: radius 2, 2048 bits, chirality enabled.
fp_generator <- rdFingerprintGenerator$GetMorganGenerator(
  radius = 2L,
  fpSize = 2048L,
  includeChirality = TRUE
)
fingerprints <- lapply(molecules, function(molecule) fp_generator$GetFingerprint(molecule))

compound_count <- nrow(compounds)
similarity_matrix <- diag(1, nrow = compound_count, ncol = compound_count)
for (i in seq_len(compound_count)) {
  if (i < compound_count) {
    for (j in seq.int(i + 1L, compound_count)) {
      similarity_value <- reticulate::py_to_r(
        DataStructs$TanimotoSimilarity(fingerprints[[i]], fingerprints[[j]])
      )
      similarity_matrix[i, j] <- as.numeric(similarity_value)
      similarity_matrix[j, i] <- as.numeric(similarity_value)
    }
  }
}

if (!isTRUE(all.equal(similarity_matrix, t(similarity_matrix), tolerance = 1e-12))) {
  stop("Similarity matrix is not symmetric.", call. = FALSE)
}
if (any(similarity_matrix < 0 | similarity_matrix > 1)) {
  stop("Tanimoto values outside the interval 0 to 1 were detected.", call. = FALSE)
}

molecular_weight <- vapply(molecules, function(molecule) {
  as.numeric(reticulate::py_to_r(Descriptors$MolWt(molecule)))
}, numeric(1))
clogp <- vapply(molecules, function(molecule) {
  as.numeric(reticulate::py_to_r(Crippen$MolLogP(molecule)))
}, numeric(1))
tpsa <- vapply(molecules, function(molecule) {
  as.numeric(reticulate::py_to_r(Descriptors$TPSA(molecule)))
}, numeric(1))

# -----------------------------
# Pairwise similarities
# -----------------------------
pair_indices <- utils::combn(seq_len(compound_count), 2L)
pairwise_table <- tibble::tibble(
  compound_id_1 = compounds$compound_id[pair_indices[1, ]],
  compound_name_1 = compounds$compound_name[pair_indices[1, ]],
  compound_class_1 = compounds$compound_class[pair_indices[1, ]],
  compound_id_2 = compounds$compound_id[pair_indices[2, ]],
  compound_name_2 = compounds$compound_name[pair_indices[2, ]],
  compound_class_2 = compounds$compound_class[pair_indices[2, ]],
  tanimoto_similarity = similarity_matrix[cbind(pair_indices[1, ], pair_indices[2, ])]
) |>
  dplyr::mutate(tanimoto_distance = 1 - .data$tanimoto_similarity) |>
  dplyr::arrange(dplyr::desc(.data$tanimoto_similarity), .data$compound_name_1, .data$compound_name_2)

candidate_indices <- which(compounds$compound_class == "network_candidate")
reference_indices <- which(compounds$compound_class == "reference_ligand")

candidate_reference_rows <- vector("list", length(candidate_indices))
for (candidate_position in seq_along(candidate_indices)) {
  candidate_index <- candidate_indices[candidate_position]
  values <- similarity_matrix[candidate_index, reference_indices]
  ordering <- order(-values, compounds$compound_name[reference_indices])
  ordered_references <- reference_indices[ordering]
  candidate_reference_rows[[candidate_position]] <- tibble::tibble(
    candidate_id = compounds$compound_id[candidate_index],
    candidate_name = compounds$compound_name[candidate_index],
    reference_id = compounds$compound_id[ordered_references],
    reference_name = compounds$compound_name[ordered_references],
    tanimoto_similarity = similarity_matrix[candidate_index, ordered_references],
    reference_rank_for_candidate = seq_along(ordered_references),
    is_nearest_reference = seq_along(ordered_references) == 1L
  )
}
candidate_reference_table <- dplyr::bind_rows(candidate_reference_rows) |>
  dplyr::arrange(.data$candidate_name, .data$reference_rank_for_candidate)

# -----------------------------
# Hierarchical clustering and chemical-space projection
# -----------------------------
distance_matrix <- 1 - similarity_matrix
diag(distance_matrix) <- 0
rownames(distance_matrix) <- compounds$compound_name
colnames(distance_matrix) <- compounds$compound_name

distance_object <- stats::as.dist(distance_matrix)
hierarchical_model <- stats::hclust(distance_object, method = "average")
cluster_ids <- stats::cutree(hierarchical_model, h = 0.60)

set.seed(42)
mds_coordinates <- stats::cmdscale(distance_object, k = 2L, eig = TRUE, add = TRUE)
if (is.null(mds_coordinates$points) || nrow(mds_coordinates$points) != compound_count) {
  stop("Two-dimensional chemical-space projection failed.", call. = FALSE)
}

cluster_table <- tibble::tibble(
  compound_id = compounds$compound_id,
  compound_name = compounds$compound_name,
  compound_class = compounds$compound_class,
  structure_source = compounds$structure_source,
  inchi_key = compounds$inchi_key,
  rdkit_parse_status = "valid",
  cluster_id = as.integer(cluster_ids),
  mds_dimension_1 = as.numeric(mds_coordinates$points[, 1]),
  mds_dimension_2 = as.numeric(mds_coordinates$points[, 2]),
  rdkit_molecular_weight = molecular_weight,
  rdkit_clogp = clogp,
  rdkit_tpsa = tpsa
) |>
  dplyr::arrange(.data$cluster_id, .data$compound_class, .data$compound_name)

nearest_reference <- candidate_reference_table |>
  dplyr::filter(.data$is_nearest_reference) |>
  dplyr::transmute(
    candidate_id,
    candidate_name,
    nearest_reference_id = .data$reference_id,
    nearest_reference_name = .data$reference_name,
    nearest_reference_tanimoto = .data$tanimoto_similarity
  )

summary_table <- nearest_reference |>
  dplyr::left_join(
    cluster_table |>
      dplyr::select(
        candidate_id = .data$compound_id,
        cluster_id,
        rdkit_molecular_weight,
        rdkit_clogp,
        rdkit_tpsa
      ),
    by = "candidate_id"
  ) |>
  dplyr::arrange(dplyr::desc(.data$nearest_reference_tanimoto), .data$candidate_name)

# -----------------------------
# Heatmap
# -----------------------------
heatmap_order <- hierarchical_model$order
heatmap_long <- expand.grid(
  row_index = seq_len(compound_count),
  column_index = seq_len(compound_count),
  KEEP.OUT.ATTRS = FALSE
)
heatmap_long$row_original <- heatmap_order[heatmap_long$row_index]
heatmap_long$column_original <- heatmap_order[heatmap_long$column_index]
heatmap_long$row_name <- factor(
  compounds$compound_name[heatmap_long$row_original],
  levels = rev(compounds$compound_name[heatmap_order])
)
heatmap_long$column_name <- factor(
  compounds$compound_name[heatmap_long$column_original],
  levels = compounds$compound_name[heatmap_order]
)
heatmap_long$similarity <- similarity_matrix[cbind(heatmap_long$row_original, heatmap_long$column_original)]

heatmap_plot <- ggplot2::ggplot(
  heatmap_long,
  ggplot2::aes(x = .data$column_name, y = .data$row_name, fill = .data$similarity)
) +
  ggplot2::geom_tile(color = "white", linewidth = 0.3) +
  ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$similarity)), size = 2.5) +
  ggplot2::scale_fill_gradient(low = "#440154", high = "#FDE725", limits = c(0, 1), name = "Tanimoto") +
  ggplot2::labs(
    title = "RESKO A16 Morgan fingerprint Tanimoto similarity",
    x = NULL,
    y = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    axis.text.x = ggplot2::element_text(angle = 55, hjust = 1),
    panel.grid = ggplot2::element_blank(),
    plot.title = ggplot2::element_text(face = "bold")
  )

safe_ggsave(heatmap_plot, outputs[["heatmap"]], width = 10, height = 8)

# -----------------------------
# Chemical-space projection
# -----------------------------
space_table <- cluster_table |>
  dplyr::mutate(
    display_class = dplyr::if_else(
      .data$compound_class == "network_candidate",
      "Network candidate",
      "Reference ligand"
    )
  )

space_plot <- ggplot2::ggplot(
  space_table,
  ggplot2::aes(
    x = .data$mds_dimension_1,
    y = .data$mds_dimension_2,
    color = .data$display_class,
    shape = .data$display_class
  )
) +
  ggplot2::geom_point(size = 4, stroke = 0.8) +
  ggplot2::geom_text(
    ggplot2::aes(label = .data$compound_name),
    hjust = -0.08,
    vjust = -0.3,
    size = 3,
    show.legend = FALSE
  ) +
  ggplot2::scale_color_manual(values = c("Network candidate" = "#E67E22", "Reference ligand" = "#6F42C1")) +
  ggplot2::scale_shape_manual(values = c("Network candidate" = 18, "Reference ligand" = 16)) +
  ggplot2::labs(
    title = "RESKO A16 chemical space based on Tanimoto distance",
    x = "Classical MDS dimension 1",
    y = "Classical MDS dimension 2",
    color = NULL,
    shape = NULL
  ) +
  ggplot2::theme_minimal(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(face = "bold"),
    legend.position = "bottom"
  )

safe_ggsave(space_plot, outputs[["chemical_space"]], width = 10, height = 7)

# -----------------------------
# RDKit structure grid
# -----------------------------
legends <- paste0(
  compounds$compound_name,
  "\n",
  gsub("_", " ", compounds$compound_class)
)
python_molecule_list <- reticulate::r_to_py(molecules, convert = FALSE)
python_legend_list <- reticulate::r_to_py(as.list(legends), convert = FALSE)
grid_image <- Draw$MolsToGridImage(
  python_molecule_list,
  molsPerRow = 2L,
  subImgSize = reticulate::tuple(500L, 360L),
  legends = python_legend_list,
  useSVG = FALSE
)

grid_temp <- tempfile(pattern = "A16_grid_", tmpdir = results_dir, fileext = ".png")
grid_image$save(grid_temp)
if (!file.exists(grid_temp) || file.info(grid_temp)$size <= 0L) {
  stop("RDKit structure grid was not created.", call. = FALSE)
}
backup_existing(outputs[["structure_grid"]])
if (!file.rename(grid_temp, outputs[["structure_grid"]])) {
  stop("Could not move the structure grid into place.", call. = FALSE)
}

# -----------------------------
# Write tables
# -----------------------------
safe_write_csv(pairwise_table, outputs[["pairwise"]])
safe_write_csv(candidate_reference_table, outputs[["candidate_reference"]])
safe_write_csv(cluster_table, outputs[["clusters"]])
safe_write_csv(summary_table, outputs[["summary"]])

# -----------------------------
# Standalone HTML report
# -----------------------------
top_pairs <- pairwise_table |>
  dplyr::slice_head(n = 10L) |>
  dplyr::select(.data$compound_name_1, .data$compound_name_2, .data$tanimoto_similarity)

cluster_report <- cluster_table |>
  dplyr::select(
    .data$cluster_id,
    .data$compound_name,
    .data$compound_class,
    .data$rdkit_molecular_weight,
    .data$rdkit_clogp,
    .data$rdkit_tpsa
  )

heatmap_uri <- base64enc::dataURI(file = outputs[["heatmap"]], mime = "image/png")
space_uri <- base64enc::dataURI(file = outputs[["chemical_space"]], mime = "image/png")
grid_uri <- base64enc::dataURI(file = outputs[["structure_grid"]], mime = "image/png")

report_html <- paste0(
  "<!doctype html><html lang='en'><head><meta charset='utf-8'>",
  "<meta name='viewport' content='width=device-width, initial-scale=1'>",
  "<title>RESKO A16 Chemical Similarity Report</title>",
  "<style>",
  "body{font-family:-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif;max-width:1200px;margin:0 auto;padding:28px;color:#202124;line-height:1.5}",
  "h1,h2{color:#263238}.note{background:#fff4e5;border-left:5px solid #e67e22;padding:14px;margin:16px 0}",
  "img{max-width:100%;height:auto;border:1px solid #ddd;border-radius:8px}",
  "table.data{border-collapse:collapse;width:100%;font-size:.9rem}table.data th,table.data td{border:1px solid #ddd;padding:7px;text-align:left}",
  "table.data th{background:#f3f4f6}small{color:#5f6368}</style></head><body>",
  "<h1>RESKO A16 Chemical Similarity Report</h1>",
  "<p><small>Generated ", html_escape(format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")), "</small></p>",
  "<h2>Analysis summary</h2>",
  "<p>Ten compounds were analysed: six RESKO network candidates and four established eEF1A reference ligands. All structures parsed successfully with RDKit.</p>",
  "<div class='note'><strong>Interpretive boundary:</strong> Chemical similarity is a prioritisation aid. It does not establish a shared binding site, mechanism, efficacy, selectivity, or toxicity.</div>",
  "<h2>Methodological note</h2>",
  "<p>Analysis SMILES preserved encoded stereochemistry. Morgan fingerprints used radius 2, 2,048 bits, and chirality. Pairwise similarity used the Tanimoto coefficient. Average-linkage hierarchical clustering used Tanimoto distance, with clusters cut at distance 0.60. The two-dimensional projection used classical multidimensional scaling on the same distance matrix. Molecular weight, cLogP, and TPSA were calculated with RDKit.</p>",
  "<h2>Nearest reference ligand for each candidate</h2>", table_to_html(summary_table),
  "<h2>Pairwise similarity interpretation</h2>",
  "<p>The ten highest observed pairwise similarities are shown below. Values are comparative within this ten-compound set and do not prove shared pharmacology.</p>",
  table_to_html(top_pairs),
  "<h2>Structural clusters</h2>", table_to_html(cluster_report),
  "<h2>Similarity heatmap</h2><img src='", heatmap_uri, "' alt='Tanimoto similarity heatmap'>",
  "<h2>Chemical-space projection</h2><img src='", space_uri, "' alt='Chemical-space projection'>",
  "<h2>Compound structure grid</h2><img src='", grid_uri, "' alt='Compound structure grid'>",
  "<h2>Limitations</h2><ul>",
  "<li>Only ten compounds are included, so clusters and chemical space describe this local comparison set.</li>",
  "<li>Two-dimensional projection loses information from the complete fingerprint-distance matrix.</li>",
  "<li>Fingerprint similarity depends on chemical representation, stereochemistry, salts, protonation, tautomerism, and source standardisation.</li>",
  "<li>Large macrocycles and smaller drug-like compounds differ strongly in size and topology, which may dominate comparisons.</li>",
  "<li>No conclusion about binding pose, target engagement, inhibition, selectivity, efficacy, or safety follows from similarity alone.</li>",
  "</ul></body></html>"
)

report_temp <- tempfile(pattern = "A16_report_", tmpdir = results_dir, fileext = ".html")
writeLines(report_html, report_temp, useBytes = TRUE)
if (!file.exists(report_temp) || file.info(report_temp)$size <= 0L) stop("HTML report was not created.", call. = FALSE)
backup_existing(outputs[["report"]])
if (!file.rename(report_temp, outputs[["report"]])) stop("Could not move HTML report into place.", call. = FALSE)

# -----------------------------
# Final output validation
# -----------------------------
for (path in outputs) {
  if (!file.exists(path) || is.na(file.info(path)$size) || file.info(path)$size <= 0L) {
    stop("Output verification failed: ", path, call. = FALSE)
  }
}

if (nrow(readr::read_csv(outputs[["pairwise"]], show_col_types = FALSE)) != 45L) {
  stop("Pairwise output must contain 45 rows.", call. = FALSE)
}
if (nrow(readr::read_csv(outputs[["candidate_reference"]], show_col_types = FALSE)) != 24L) {
  stop("Candidate-reference output must contain 24 rows.", call. = FALSE)
}
if (nrow(readr::read_csv(outputs[["clusters"]], show_col_types = FALSE)) != 10L) {
  stop("Cluster output must contain 10 rows.", call. = FALSE)
}
if (nrow(readr::read_csv(outputs[["summary"]], show_col_types = FALSE)) != 6L) {
  stop("Similarity summary must contain 6 rows.", call. = FALSE)
}

cat("A16B chemical similarity analysis completed in R.\n")
cat("RDKit Python: ", active_python, "\n", sep = "")
cat("Compounds analysed: ", compound_count, "\n", sep = "")
cat("Pairwise comparisons: ", nrow(pairwise_table), "\n", sep = "")
cat("Candidate-reference comparisons: ", nrow(candidate_reference_table), "\n", sep = "")
cat("Structural clusters at distance 0.60: ", length(unique(cluster_ids)), "\n", sep = "")
cat("Verified outputs:\n")
for (path in outputs) {
  cat("  ", sub(paste0("^", project_root, "/"), "", path), " (", file.info(path)$size, " bytes)\n", sep = "")
}
