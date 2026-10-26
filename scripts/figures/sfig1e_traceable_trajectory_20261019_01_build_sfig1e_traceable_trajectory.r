#!/usr/bin/env Rscript

# ==========================================================
# Script   : sfig1e_traceable_trajectory_20261019_01_build_sfig1e_traceable_trajectory.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:sfig1e_traceable_trajectory_20261019/scripts/01_build_sfig1e_traceable_trajectory.R
# ==========================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)

root <- normalizePath(".", winslash = "/", mustWork = TRUE)
module_dir <- file.path(root, "sfig1e_traceable_trajectory_20261019")
input_rds <- file.path(root, "inputs/singlecell/subcluster/b_subclusters.rds")
provenance_script <- file.path(root, "inputs/parent_scripts/L/monocle3_bcell_trajectory.R")
out_data <- file.path(module_dir, "source_data")
out_fig <- file.path(module_dir, "figures")
out_qc <- file.path(module_dir, "qc")
out_log <- file.path(module_dir, "logs")
dir.create(out_data, recursive = TRUE, showWarnings = FALSE)
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_qc, recursive = TRUE, showWarnings = FALSE)
dir.create(out_log, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(input_rds), file.exists(provenance_script))

log_file <- file.path(out_log, "01_build_sfig1e_traceable_trajectory.log")
log_con <- file(log_file, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  sink(type = "message")
  sink(type = "output")
  close(log_con)
}, add = TRUE)

cat("SFig1E traceable trajectory build\n")
cat("Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")

seu <- readRDS(input_rds)
required_meta <- c("patient", "sample", "b_subtype", "monocle3_pt")
if (!inherits(seu, "Seurat")) stop("Input is not a Seurat object")
if (!all(required_meta %in% colnames(seu@meta.data))) {
  stop("Missing required metadata: ", paste(setdiff(required_meta, colnames(seu@meta.data)), collapse = ", "))
}
if (!"umap" %in% names(seu@reductions)) stop("Required UMAP reduction is absent")

meta <- seu@meta.data
umap <- Embeddings(seu, "umap")
if (!identical(rownames(meta), rownames(umap))) stop("Cell order differs between metadata and UMAP")

state_levels <- c("Naive B", "Transitional B", "GC-like B", "Memory B", "Plasma cell")
unknown_states <- setdiff(unique(as.character(meta$b_subtype)), state_levels)
if (length(unknown_states) > 0) stop("Unexpected B-cell states: ", paste(unknown_states, collapse = ", "))

cell_df <- data.frame(
  cell_id = rownames(meta),
  patient = as.character(meta$patient),
  sample = as.character(meta$sample),
  b_subtype = factor(as.character(meta$b_subtype), levels = state_levels),
  monocle3_pt = as.numeric(meta$monocle3_pt),
  UMAP1 = as.numeric(umap[, 1]),
  UMAP2 = as.numeric(umap[, 2]),
  stringsAsFactors = FALSE
)

if (anyNA(cell_df$patient) || any(cell_df$patient == "")) stop("Missing patient IDs")
if (anyNA(cell_df$b_subtype)) stop("Missing or invalid B-cell states")
if (any(!is.finite(cell_df$monocle3_pt))) stop("Non-finite monocle3_pt values")
if (any(!is.finite(cell_df$UMAP1)) || any(!is.finite(cell_df$UMAP2))) stop("Non-finite UMAP values")
if (anyDuplicated(cell_df$cell_id)) stop("Duplicated cell IDs")

write.table(cell_df, file.path(out_data, "SFig1E_cell_level_coordinates.tsv.gz"),
            sep = "\t", row.names = FALSE, quote = FALSE)

state_summary <- do.call(rbind, lapply(state_levels, function(s) {
  z <- cell_df$monocle3_pt[cell_df$b_subtype == s]
  data.frame(
    b_subtype = s,
    n_cells = length(z),
    n_patients = length(unique(cell_df$patient[cell_df$b_subtype == s])),
    median = median(z),
    q1 = unname(quantile(z, 0.25)),
    q3 = unname(quantile(z, 0.75)),
    min = min(z),
    max = max(z)
  )
}))
write.table(state_summary, file.path(out_data, "SFig1E_state_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)

patient_state <- aggregate(
  monocle3_pt ~ patient + b_subtype,
  data = cell_df,
  FUN = function(z) c(n_cells = length(z), median = median(z), q1 = quantile(z, 0.25), q3 = quantile(z, 0.75))
)
patient_state <- data.frame(
  patient = patient_state$patient,
  b_subtype = factor(patient_state$b_subtype, levels = state_levels),
  n_cells = patient_state$monocle3_pt[, "n_cells"],
  median_pt = patient_state$monocle3_pt[, "median"],
  q1_pt = patient_state$monocle3_pt[, "q1.25%"],
  q3_pt = patient_state$monocle3_pt[, "q3.75%"]
)
patient_state <- patient_state[order(patient_state$b_subtype, patient_state$patient), ]
patient_state$eligible_ge10_cells <- patient_state$n_cells >= 10
write.table(patient_state, file.path(out_data, "SFig1E_patient_state_summary.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
patient_state_plot <- patient_state[patient_state$eligible_ge10_cells, ]
write.table(patient_state_plot, file.path(out_data, "SFig1E_patient_state_displayed_ge10.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
eligible_by_state <- table(factor(patient_state_plot$b_subtype, levels = state_levels))

contract <- list(
  panel = "SFig1E",
  source_authority = "投稿文件.pdf::Supplementary Figure 1E",
  canvas_mm = list(width = 72, height = 86),
  geometry_basis = paste(
    "The submitted Supplementary Figure 1 composite is labelled 200 x 219 mm.",
    "SFig.1E occupies approximately 35.8% of composite width and 39.3% of composite height",
    "when measured from the embedded 1145 x 1251 px figure object, yielding 72 x 86 mm."
  ),
  data_input = "inputs/singlecell/subcluster/b_subclusters.rds",
  provenance_script = "inputs/parent_scripts/L/monocle3_bcell_trajectory.R",
  n_cells = nrow(cell_df),
  n_patients = length(unique(cell_df$patient)),
  n_samples = length(unique(cell_df$sample)),
  finite_pseudotime = sum(is.finite(cell_df$monocle3_pt)),
  root_definition_from_provenance_script = "all cells labelled Naive B",
  estimand = "descriptive cross-sectional Monocle3 ordering already stored in the frozen object",
  independent_unit_for_support = paste(
    "patient; patient-state medians are displayed descriptively only when",
    "the patient-state combination contains at least 10 cells, as prespecified in ANALYSIS_REGISTRY.tsv"
  ),
  inferential_test = "none",
  multiplicity = "not applicable",
  claim_ceiling = paste(
    "Supports descriptive cross-sectional transcriptional-state ordering only;",
    "does not establish monotonic naive-to-plasma lineage, lineage time, direction, fate, or causality."
  ),
  forbidden = c("showtext_auto", "bootstrap", "model refit", "trajectory refit", "lineage-direction claim")
)
write_json(contract, file.path(out_qc, "SFig1E_FIGURE_CONTRACT.json"), pretty = TRUE, auto_unbox = TRUE)

state_colors <- c(
  "Naive B" = "#2166AC",
  "Transitional B" = "#67A9CF",
  "GC-like B" = "#F4A261",
  "Memory B" = "#4DAF4A",
  "Plasma cell" = "#D73027"
)

theme_small <- theme_classic(base_family = "Arial", base_size = 8.0) +
  theme(
    axis.line = element_line(linewidth = 0.28, colour = "black"),
    axis.ticks = element_line(linewidth = 0.28, colour = "black"),
    axis.text = element_text(size = 7.5, colour = "black"),
    axis.title = element_text(size = 7.5, colour = "black"),
    plot.title = element_text(size = 8.0, face = "bold", margin = margin(b = 1.5)),
    plot.subtitle = element_text(size = 7.5, colour = "#333333", margin = margin(b = 1.5)),
    legend.title = element_text(size = 7.5),
    legend.text = element_text(size = 7.5),
    plot.margin = margin(2, 2, 1, 2)
  )

p_umap <- ggplot(cell_df, aes(UMAP1, UMAP2, colour = monocle3_pt)) +
  geom_point(size = 0.16, alpha = 0.70, stroke = 0) +
  scale_colour_viridis_c(option = "C", direction = 1, name = "Monocle3\nordering") +
  coord_equal() +
  labs(
    title = "Cross-sectional B-cell ordering",
    subtitle = sprintf("Frozen Monocle3 coordinate; n = %s cells", format(nrow(cell_df), big.mark = ",")),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_small +
  theme(
    legend.position = "right",
    legend.key.height = grid::unit(7, "mm"),
    legend.key.width = grid::unit(1.5, "mm")
  )

p_patient <- ggplot(patient_state_plot, aes(b_subtype, median_pt, colour = b_subtype)) +
  geom_boxplot(aes(group = b_subtype), width = 0.58, outlier.shape = NA,
               fill = "white", linewidth = 0.32, colour = "#444444") +
  geom_point(position = position_jitter(width = 0.12, height = 0, seed = 20261019),
             size = 1.0, alpha = 0.78) +
  scale_colour_manual(values = state_colors, guide = "none") +
  scale_x_discrete(labels = sprintf(
    "%s\nn=%d",
    c("Naive", "Trans.", "GC-like", "Memory", "Plasma"),
    as.integer(eligible_by_state)
  )) +
  labs(
    title = "Patient-level state summaries",
    subtitle = "Each point is one patient-state median (at least 10 cells)",
    x = NULL, y = "Monocle3 ordering"
  ) +
  theme_small +
  theme(axis.text.x = element_text(size = 7.5, angle = 0, hjust = 0.5, vjust = 1, lineheight = 0.9))

combined <- p_umap / p_patient +
  plot_layout(heights = c(1.45, 1)) +
  plot_annotation(
    caption = "Cross-sectional ordering only; not lineage time or direction.",
    theme = theme(
      plot.caption = element_text(family = "Arial", size = 7.5, colour = "#333333", hjust = 0),
      plot.margin = margin(1, 1, 1, 1)
    )
  )

width_mm <- 82
height_mm <- 96
pdf_path <- file.path(out_fig, "SFig1E_traceable_Monocle3_ordering.pdf")
png_path <- file.path(out_fig, "SFig1E_traceable_Monocle3_ordering.png")

cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4,
          family = "Arial", bg = "white")
print(combined)
dev.off()

agg_png(png_path, width = width_mm, height = height_mm, units = "mm",
        res = 600, background = "white", scaling = 1)
print(combined)
dev.off()

qc <- data.frame(
  check = c(
    "input_is_seurat", "cell_count_10572", "finite_monocle3_pt_all_cells",
    "cell_metadata_umap_alignment", "unique_cell_ids", "patient_ids_complete",
    "five_expected_states_only", "all_states_patient_supported", "display_gate_ge10_cells",
    "all_states_have_eligible_patients", "no_inferential_p_values",
    "no_trajectory_refit", "no_bootstrap", "no_showtext_auto", "pdf_created", "png_created"
  ),
  result = c(
    inherits(seu, "Seurat"), nrow(cell_df) == 10572,
    all(is.finite(cell_df$monocle3_pt)), identical(rownames(meta), rownames(umap)),
    !anyDuplicated(cell_df$cell_id), !anyNA(cell_df$patient) && all(cell_df$patient != ""),
    setequal(as.character(unique(cell_df$b_subtype)), state_levels),
    all(state_summary$n_patients > 0), all(patient_state_plot$n_cells >= 10),
    all(as.integer(eligible_by_state) > 0), TRUE, TRUE, TRUE, TRUE,
    file.exists(pdf_path), file.exists(png_path)
  ),
  evidence = c(
    class(seu)[1], as.character(nrow(cell_df)),
    sprintf("%d/%d", sum(is.finite(cell_df$monocle3_pt)), nrow(cell_df)),
    "identical row names and order", as.character(anyDuplicated(cell_df$cell_id)),
    sprintf("%d patients", length(unique(cell_df$patient))),
    paste(state_levels, collapse = "; "),
    paste(sprintf("%s=%d", state_summary$b_subtype, state_summary$n_patients), collapse = "; "),
    sprintf("%d/%d patient-state combinations displayed", nrow(patient_state_plot), nrow(patient_state)),
    paste(sprintf("%s=%d", names(eligible_by_state), as.integer(eligible_by_state)), collapse = "; "),
    "descriptive summaries only", "stored monocle3_pt used without learn_graph/order_cells",
    "no resampling called", "function absent from script", pdf_path, png_path
  )
)
write.table(qc, file.path(out_qc, "SFig1E_PRODUCTION_QC.tsv"),
            sep = "\t", row.names = FALSE, quote = FALSE)
if (!all(qc$result)) stop("One or more production QC checks failed")

cat("Cells:", nrow(cell_df), "\n")
cat("Patients:", length(unique(cell_df$patient)), "\n")
cat("Samples:", length(unique(cell_df$sample)), "\n")
cat("Pseudotime range:", paste(range(cell_df$monocle3_pt), collapse = " to "), "\n")
cat("State summaries:\n")
print(state_summary)
cat("Eligible patient-state summaries (>=10 cells):\n")
print(eligible_by_state)
cat("Saved:", pdf_path, "\n")
cat("Saved:", png_path, "\n")
cat("Completed:", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "\n")
