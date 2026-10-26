#!/usr/bin/env Rscript

# ==========================================================
# Script   : 05_bcell_cross_sectional_ordering.R
# Purpose  : Reproduce the cross-sectional B-cell ordering analysis
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: cross-sectional ordering; Supplementary Fig. 1E
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/01_single_cell/monocle3_bcell_trajectory.R
# ==========================================================
# monocle3_bcell_trajectory.R
# B cell maturation trajectory with Monocle3 (replaces Slingshot)
# Input:  results/subcluster/b_subclusters.rds
# Output: results/revision_restructure/sfig2D_bcell_monocle3.{png,pdf,tiff}
#         results/subcluster/b_subclusters.rds  (adds monocle3_pt column)

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(SeuratWrappers)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(ragg)
  library(igraph)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "results/revision_restructure"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── B cell subtype color palette (consistent with Fig1E) ──────────────────────
b_colors <- c(
  "Naive B"        = "#2171B5",
  "Transitional B" = "#6BAED6",
  "GC-like B"      = "#E28E2C",
  "Memory B"       = "#74C476",
  "Plasma cell"    = "#D24B40"
)

base_theme <- theme_classic(base_size = 6.5, base_family = "sans") +
  theme(
    axis.line    = element_line(linewidth = 0.35),
    axis.ticks   = element_line(linewidth = 0.35),
    axis.title   = element_text(size = 6.0),
    axis.text    = element_text(size = 5.5),
    legend.title = element_text(size = 6.0),
    legend.text  = element_text(size = 5.5),
    legend.key.size = unit(3, "mm"),
    plot.title   = element_text(size = 7.0, face = "bold"),
    plot.subtitle = element_text(size = 5.5, colour = "grey45"),
    panel.grid   = element_blank()
  )

# ── 1. Load ────────────────────────────────────────────────────────────────────
log_msg("Loading B cell subclusters...")
seu <- readRDS("results/subcluster/b_subclusters.rds")
log_msg(sprintf("  n=%d cells | subtypes: %s",
                ncol(seu), paste(names(table(seu$b_subtype)), collapse=", ")))

# ── 2. Convert to Monocle3 CellDataSet ───────────────────────────────────────
log_msg("Converting to Monocle3 CellDataSet...")
cds <- as.cell_data_set(seu)

# Inject Seurat UMAP into the CDS reducedDims
reducedDims(cds)[["UMAP"]] <- Embeddings(seu, "umap")

# Force single partition (prevents disconnected mini-trajectories)
one_partition <- factor(rep(1L, ncol(cds)), levels = "1")
names(one_partition) <- colnames(cds)
cds@clusters[["UMAP"]]$partitions <- one_partition

# Pass Seurat cluster IDs so Monocle3 has a cluster reference
seurat_clust <- seu$seurat_clusters
names(seurat_clust) <- colnames(cds)
cds@clusters[["UMAP"]]$clusters <- seurat_clust

# ── 3. Learn trajectory graph ─────────────────────────────────────────────────
log_msg("Learning principal graph (use_partition=FALSE, prune short branches)...")
cds <- learn_graph(
  cds,
  use_partition = FALSE,
  close_loop    = FALSE,
  learn_graph_control = list(
    ncenter               = 12,    # 12 principal nodes for 5 subtypes → smooth trajectory
    minimal_branch_len    = 5,     # post-pruning short branch threshold
    geodesic_distance_ratio = 0.33,
    prune_graph           = TRUE
  )
)

# ── 4. Order cells: root = Naive B centroid ───────────────────────────────────
log_msg("Setting root to Naive B cells and ordering...")
naive_cells <- colnames(cds)[colData(cds)$b_subtype == "Naive B"]
cds <- order_cells(cds, root_cells = naive_cells)

# ── 5. Extract pseudotime → write back to Seurat ─────────────────────────────
pt <- pseudotime(cds)
seu$monocle3_pt <- pt[colnames(seu)]
finite_pt <- pt[is.finite(pt)]
log_msg(sprintf("  Pseudotime: min=%.2f, max=%.2f, median=%.2f  (%.1f%% finite)",
                min(finite_pt), max(finite_pt), median(finite_pt),
                100 * length(finite_pt) / length(pt)))

log_msg("Saving updated Seurat object (monocle3_pt added)...")
saveRDS(seu, "results/subcluster/b_subclusters.rds")

# Validate trajectory direction
pt_by_subtype <- tapply(seu$monocle3_pt[is.finite(seu$monocle3_pt)],
                        seu$b_subtype[is.finite(seu$monocle3_pt)], median)
log_msg("Median pseudotime by subtype (expected: Naive B lowest, Plasma cell highest):")
for (nm in names(sort(pt_by_subtype))) {
  log_msg(sprintf("  %-20s %.4f", nm, pt_by_subtype[nm]))
}

# ── 6. Extract principal graph edges for ggplot overlay ──────────────────────
dp_mst     <- t(cds@principal_graph_aux[["UMAP"]]$dp_mst)   # nodes × 2
edge_list  <- igraph::as_data_frame(cds@principal_graph[["UMAP"]], what = "edges")
seg_df <- data.frame(
  x    = dp_mst[edge_list$from, 1], y    = dp_mst[edge_list$from, 2],
  xend = dp_mst[edge_list$to,   1], yend = dp_mst[edge_list$to,   2]
)

# Cell data frame for ggplot
cell_df <- data.frame(
  UMAP1     = Embeddings(seu, "umap")[, 1],
  UMAP2     = Embeddings(seu, "umap")[, 2],
  b_subtype = factor(seu$b_subtype, levels = names(b_colors)),
  pt        = seu$monocle3_pt,
  stringsAsFactors = FALSE
)
cell_df_finite <- cell_df[is.finite(cell_df$pt), ]
# Clamp pseudotime to 5–95th percentile for visible colour gradient
pt_lo <- quantile(cell_df_finite$pt, 0.05)
pt_hi <- quantile(cell_df_finite$pt, 0.95)
cell_df_finite$pt_clamped <- pmax(pt_lo, pmin(pt_hi, cell_df_finite$pt))

# ── 7. Build figure ───────────────────────────────────────────────────────────
log_msg("Drawing panels...")

# Panel A: pseudotime gradient
pA <- ggplot(cell_df_finite, aes(UMAP1, UMAP2, colour = pt_clamped)) +
  geom_point(size = 0.20, alpha = 0.78, stroke = 0) +
  geom_segment(data = seg_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "black", linewidth = 0.55, inherit.aes = FALSE) +
  scale_colour_viridis_c(
    option = "D", name = "Pseudotime",
    limits = c(pt_lo, pt_hi),
    guide  = guide_colorbar(barwidth = unit(2.2, "cm"), barheight = unit(1.8, "mm"),
                            title.position = "top", title.hjust = 0.5)
  ) +
  labs(title    = "Monocle3 pseudotime",
       subtitle = "root: Naive B → GC-like / Memory → Plasma cell",
       x = "UMAP 1", y = "UMAP 2") +
  base_theme +
  theme(legend.position = "bottom",
        legend.margin   = margin(t = 0))

# Panel B: subtype identity
pB <- ggplot(cell_df, aes(UMAP1, UMAP2, colour = b_subtype)) +
  geom_point(size = 0.20, alpha = 0.78, stroke = 0) +
  geom_segment(data = seg_df,
               aes(x = x, y = y, xend = xend, yend = yend),
               colour = "black", linewidth = 0.55, inherit.aes = FALSE) +
  scale_colour_manual(
    values = b_colors, name = "B-cell subtype",
    guide  = guide_legend(override.aes = list(size = 2.2, alpha = 1),
                          ncol = 1, byrow = FALSE)
  ) +
  labs(title    = "B-cell subtypes",
       subtitle = sprintf("n = %d cells  |  5 maturation states", nrow(cell_df)),
       x = "UMAP 1", y = "UMAP 2") +
  base_theme +
  theme(legend.position = "right")

p_out <- pA + pB +
  plot_annotation(
    title = "B-cell maturation trajectory (Monocle3, principal graph)",
    theme = theme(plot.title = element_text(size = 7.5, face = "bold",
                                            family = "sans"))
  )

# ── 8. Save ───────────────────────────────────────────────────────────────────
base <- file.path(OUT, "sfig2D_bcell_monocle3")
w_in <- 183 / 25.4
h_in <-  90 / 25.4

log_msg("Saving PNG...")
ragg::agg_png(paste0(base, ".png"), width = w_in, height = h_in,
              units = "in", res = 300, bg = "white")
print(p_out); dev.off()

log_msg("Saving PDF...")
grDevices::cairo_pdf(paste0(base, ".pdf"), width = w_in, height = h_in)
print(p_out); dev.off()

log_msg("Saving TIFF...")
ragg::agg_tiff(paste0(base, ".tiff"), width = w_in, height = h_in,
               units = "in", res = 600, bg = "white")
print(p_out); dev.off()

log_msg(sprintf("Saved: %s.{png,pdf,tiff}", base))
log_msg("=== monocle3_bcell_trajectory.R complete ===")
