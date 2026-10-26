#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1c_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1C_source.R
# ==========================================================
# fig1C_final.R — 基于 fig1c_composition_bar.R 修改
# 修改项：输出路径→final/fig1/、画布60×105mm、字体统一（标题9pt居中无副标题、坐标轴7-8pt、legend7pt）
# 其他代码严格保持原脚本不变
#
# Input:  results/01_scRNAseq/annotation/seurat_annotated.rds
# Output: final/fig1/fig1C_alluvial_composition.{png,pdf}

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(ragg)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "final/fig1"

# ── Canonical cell-type palette ──────────────────────────────────────────────
CELL_COLORS <- c(
  "B_cell"          = "#2171B5",
  "Plasma_cell"     = "#6BAED6",
  "CD8_T_Cytotoxic" = "#D24B40",
  "CD4_T_Helper"    = "#E28E2C",
  "Treg"            = "#F0C040",
  "Cycling"         = "#9E7BB5",
  "Macrophage"      = "#8C4B2F",
  "DC_pDC"          = "#A1C45A",
  "Mast_cell"       = "#FDAE6B",
  "CAF"             = "#C4A882",
  "Endothelial"     = "#74C476",
  "Pericyte"        = "#BDBDBD",
  "Tumor_cell"      = "#BCBDDC"
)
ct_order <- names(CELL_COLORS)

# ── Load and compute proportions ──────────────────────────────────────────────
log_msg("Loading Seurat object...")
seu <- readRDS("results/01_scRNAseq/annotation/seurat_annotated.rds")
seu$cell_type <- factor(seu$cell_type, levels = ct_order)
log_msg(sprintf("  %d cells | %d patients", ncol(seu), length(unique(seu$patient))))

prop_df <- seu@meta.data %>%
  group_by(patient, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient) %>%
  mutate(pct = n / sum(n) * 100) %>%
  ungroup()

b_order <- prop_df %>%
  filter(cell_type == "B_cell") %>%
  arrange(desc(pct)) %>%
  pull(patient)

prop_df$patient   <- factor(prop_df$patient, levels = b_order)
prop_df$cell_type <- factor(prop_df$cell_type, levels = rev(ct_order))

# ── Stacked bar ───────────────────────────────────────────────────────────────
p1C <- ggplot(prop_df, aes(patient, pct, fill = cell_type)) +
  geom_col(width = 0.85) +
  scale_fill_manual(
    values = CELL_COLORS,
    guide  = guide_legend(reverse = TRUE, ncol = 1,
                          override.aes = list(size = 3))
  ) +
  scale_y_continuous(
    labels = scales::percent_format(scale = 1),
    expand = c(0, 0)
  ) +
  labs(
    title = "Cell-type composition across 18 HNSCC patients",
    x     = "Patient",
    y     = "Proportion (%)",
    fill  = "Cell type"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line        = element_line(linewidth = 0.35),
    axis.ticks       = element_line(linewidth = 0.35),
    axis.title       = element_text(size = 7),
    axis.text.y      = element_text(size = 7),
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 7),
    legend.title     = element_text(size = 7),
    legend.text      = element_text(size = 7),
    legend.key.size  = unit(3, "mm"),
    plot.title       = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle    = element_blank(),
    panel.grid       = element_blank(),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ── Save: 60×105 mm ──────────────────────────────────────────────────────────
w_mm <- 105; h_mm <- 60
w <- w_mm / 25.4; h <- h_mm / 25.4

log_msg("Saving fig1C...")
ragg::agg_png(file.path(OUT, "fig1C_alluvial_composition.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p1C); dev.off()

grDevices::cairo_pdf(file.path(OUT, "fig1C_alluvial_composition.pdf"),
                     width = w, height = h, family = "Arial")
print(p1C); dev.off()

log_msg("Saved: final/fig1/fig1C_alluvial_composition.{png,pdf} [60×105mm]")
log_msg("=== fig1C_final.R complete ===")
