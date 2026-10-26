#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1a_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1A_source.R
# ==========================================================
# fig1A_final.R — 基于 fig1a_umap_v2.R 修改
# 修改项：输出路径→final/fig1/、画布55×65mm、字体统一（标题9pt居中无副标题、坐标轴7-8pt、legend7pt）
# 其他代码严格保持原脚本不变
#
# Input:  results/01_scRNAseq/annotation/seurat_annotated.rds
# Output: final/fig1/fig1A_umap_celltype.{png,pdf}

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(ggrepel); library(ragg); library(ggrastr)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "final/fig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── More saturated canonical palette ─────────────────────────────────────────
CELL_COLORS <- c(
  "B_cell"          = "#1560A8",
  "Plasma_cell"     = "#4D9DD4",
  "CD8_T_Cytotoxic" = "#C8351A",
  "CD4_T_Helper"    = "#D97E12",
  "Treg"            = "#D4A800",
  "Cycling"         = "#8B5EA8",
  "Macrophage"      = "#7B3520",
  "DC_pDC"          = "#86B530",
  "Mast_cell"       = "#F08030",
  "CAF"             = "#A8895A",
  "Endothelial"     = "#38A038",
  "Pericyte"        = "#909090",
  "Tumor_cell"      = "#8080C0"
)
ct_order <- names(CELL_COLORS)

# ── Load Seurat ───────────────────────────────────────────────────────────────
log_msg("Loading Seurat object...")
seu <- readRDS("results/01_scRNAseq/annotation/seurat_annotated.rds")
seu$cell_type <- factor(seu$cell_type, levels = ct_order)
n <- ncol(seu)
log_msg(sprintf("  %d cells | %d patients", n, length(unique(seu$patient))))

# ── Point / alpha parameters ──────────────────────────────────────────────────
if (n < 30000) {
  pt <- 0.5; al <- 0.90
} else if (n < 100000) {
  pt <- 0.35; al <- 0.82
} else {
  pt <- 0.22; al <- 0.75
}

# ── Cluster centroids for labels ─────────────────────────────────────────────
umap_df <- as.data.frame(Embeddings(seu, "umap"))
umap_df$cell_type <- seu$cell_type
centroids <- umap_df %>%
  group_by(cell_type) %>%
  summarise(U1 = median(umap_1), U2 = median(umap_2), .groups = "drop")

# ── Plot ──────────────────────────────────────────────────────────────────────
p1A <- ggplot(umap_df, aes(umap_1, umap_2, colour = cell_type)) +
  geom_point(size = pt, alpha = al, stroke = 0) +
  ggrepel::geom_text_repel(
    data = centroids,
    aes(U1, U2, label = gsub("_", " ", cell_type)),
    colour = "black",
    size = 2.0,
    fontface = "bold",
    family = "Arial",
    segment.size = 0.25,
    segment.colour = "grey50",
    min.segment.length = 0.2,
    box.padding = 0.25,
    force = 1.5,
    inherit.aes = FALSE
  ) +
  scale_colour_manual(values = CELL_COLORS, name = "Cell type",
                      guide = guide_legend(
                        override.aes = list(size = 3, alpha = 1),
                        ncol = 1
                      )) +
  labs(
    title = sprintf("UMAP: %s high-quality cells, 13 annotated types",
                    format(n, big.mark = ",")),
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line      = element_line(linewidth = 0.35),
    axis.ticks     = element_blank(),
    axis.text      = element_blank(),
    axis.title     = element_text(size = 7),
    legend.position = "right",
    legend.title   = element_text(size = 7),
    legend.text    = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.title     = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle  = element_blank(),
    panel.grid     = element_blank()
  )
p1A$layers[[1]]$aes_params$alpha <- al

# 光栅化点云层
p1A_rast <- ggrastr::rasterize(p1A, layers = "Point", dpi = 600)

# ── Save: 55×65 mm ───────────────────────────────────────────────────────────
w_mm <- 65; h_mm <- 55
w <- w_mm / 25.4; h <- h_mm / 25.4
ragg::agg_png(file.path(OUT, "fig1A_umap_celltype.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p1A_rast); dev.off()
grDevices::cairo_pdf(file.path(OUT, "fig1A_umap_celltype.pdf"), width = w, height = h, family = "Arial")
print(p1A_rast); dev.off()
log_msg("Saved: final/fig1/fig1A_umap_celltype.{png,pdf} [55×65mm]")
log_msg("=== fig1A_final.R complete ===")
