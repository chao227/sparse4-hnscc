#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1b_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1B_source.R
# ==========================================================
# fig1B_final.R — 基于 fig1bf_violin_hv.R 修改（仅 Fig1B vertical 版本）
# 修改项：输出路径→final/fig1/、画布135×95mm、字体统一（标题9pt居中无副标题、坐标轴7-8pt、legend7pt、基因名7pt）
# 其他代码严格保持原脚本不变
#
# Input:  results/01_scRNAseq/annotation/seurat_annotated.rds
# Output: final/fig1/fig1B_stacked_violin_celltypes.{png,pdf}

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(ragg)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "final/fig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── Canonical palette ────────────────────────────────────────────────────────
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

# ── Fig1B — 13-cell-type stacked violin (vertical) ──────────────────────────
log_msg("Loading main Seurat object...")
seu <- readRDS("results/01_scRNAseq/annotation/seurat_annotated.rds")
seu$cell_type <- factor(seu$cell_type, levels = ct_order)

markers <- c(
  "MS4A1","CD79A","JCHAIN","MZB1",
  "CD8A","GZMB","PRF1",
  "CD4","IL7R","FOXP3",
  "MKI67","TOP2A",
  "CD14","LYZ","FCGR3A",
  "LILRA4","CLEC4C",
  "TPSAB1","CPA3",
  "FAP","ACTA2",
  "PECAM1","VWF",
  "RGS5","PDGFRB",
  "EPCAM","KRT5"
)
avail_genes <- intersect(markers, rownames(seu))
log_msg(sprintf("Available marker genes: %d / %d", length(avail_genes), length(markers)))

log_msg("Extracting expression data for Fig1B...")
expr_mat <- as.matrix(GetAssayData(seu, layer = "data")[avail_genes, ])
expr_long <- as.data.frame(t(expr_mat))
expr_long$cell_type <- as.character(seu$cell_type)

df_vln <- expr_long %>%
  pivot_longer(-cell_type, names_to = "gene", values_to = "expr") %>%
  mutate(
    gene      = factor(gene, levels = avail_genes),
    cell_type = factor(cell_type, levels = ct_order)
  )

# ── Vertical: genes = facet rows, cell types = x-axis ─────────────────────────
log_msg("Fig1B: vertical version...")
p1B_v <- ggplot(df_vln, aes(x = cell_type, y = expr, fill = cell_type)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
  facet_wrap(~gene, ncol = 1, scales = "free_y",
             strip.position = "right") +
  scale_fill_manual(values = CELL_COLORS) +
  scale_y_continuous(breaks = NULL) +
  scale_x_discrete(labels = function(x) gsub("_", " ", x)) +
  labs(x = NULL, y = "Normalised expression") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    strip.background = element_blank(),
    strip.text.y     = element_text(size = 7, face = "italic", angle = 0, hjust = 0),
    axis.text.x      = element_text(size = 7, angle = 45, hjust = 1, vjust = 1, lineheight = 1),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_text(size = 7),
    legend.position  = "none",
    panel.grid       = element_blank(),
    panel.spacing.y  = unit(0.3, "mm"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ── Save: 135×95 mm ──────────────────────────────────────────────────────────
w_mm <- 95; h_mm <- 135
w <- w_mm / 25.4; h <- h_mm / 25.4
ragg::agg_png(file.path(OUT, "fig1B_stacked_violin_celltypes.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p1B_v); dev.off()
grDevices::cairo_pdf(file.path(OUT, "fig1B_stacked_violin_celltypes.pdf"),
                     width = w, height = h, family = "Arial")
print(p1B_v); dev.off()
log_msg("Saved: final/fig1/fig1B_stacked_violin_celltypes.{png,pdf} [135×95mm]")
log_msg("=== fig1B_final.R complete ===")
