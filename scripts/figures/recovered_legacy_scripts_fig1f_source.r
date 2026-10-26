#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1f_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1F_source.R
# ==========================================================
# fig1F_final.R — 基于 fig1bf_violin_hv.R 修改（仅 Fig1F vertical 版本）
# 修改项：输出路径→final/fig1/、画布115×70mm、字体统一（标题9pt居中无副标题、坐标轴7-8pt、legend7pt、基因名7pt）
# 其他代码严格保持原脚本不变
#
# Input:  results/subcluster/b_subclusters.rds
# Output: final/fig1/fig1F_stacked_violin_bcell_states.{png,pdf}

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr); library(ragg)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "final/fig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── B-cell subtype palette ────────────────────────────────────────────────────
B_COLORS <- c(
  "Naive B"        = "#1560A8",
  "Transitional B" = "#4D9DD4",
  "GC-like B"      = "#D97E12",
  "Memory B"       = "#38A038",
  "Plasma cell"    = "#C8351A"
)

# ── Fig1F — B-cell subtype stacked violin (vertical) ──────────────────────────
log_msg("Loading B-cell subclusters...")
b_seu <- readRDS("results/subcluster/b_subclusters.rds")
b_seu$b_subtype <- factor(b_seu$b_subtype, levels = names(B_COLORS))

b_markers <- c(
  "TCL1A","IGHD","IGHM",
  "IL4R","FCER2",
  "BCL6","AICDA","RGS13",
  "CD27","TNFRSF13B",
  "JCHAIN","MZB1","XBP1","SDC1",
  "MS4A1","CD79A","CXCR5"
)
avail_b <- intersect(b_markers, rownames(b_seu))
log_msg(sprintf("Available B-cell markers: %d / %d", length(avail_b), length(b_markers)))

log_msg("Extracting expression data for Fig1F...")
b_mat <- as.matrix(GetAssayData(b_seu, layer = "data")[avail_b, ])
b_expr <- as.data.frame(t(b_mat))
b_expr$b_subtype <- as.character(b_seu$b_subtype)

b_long <- b_expr %>%
  pivot_longer(-b_subtype, names_to = "gene", values_to = "expr") %>%
  mutate(
    gene      = factor(gene, levels = avail_b),
    b_subtype = factor(b_subtype, levels = names(B_COLORS))
  )

# ── Vertical: genes = facet rows, cell types = x-axis ─────────────────────────
log_msg("Fig1F: vertical version...")
p1F_v <- ggplot(b_long, aes(x = b_subtype, y = expr, fill = b_subtype)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.15) +
  facet_wrap(~gene, ncol = 1, scales = "free_y",
             strip.position = "right") +
  scale_fill_manual(values = B_COLORS) +
  scale_y_continuous(breaks = NULL) +
  labs(x = NULL, y = "Normalised expression") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line        = element_line(linewidth = 0.3),
    axis.ticks       = element_line(linewidth = 0.3),
    strip.background = element_blank(),
    strip.text.y     = element_text(size = 7, face = "italic", angle = 0, hjust = 0),
    axis.text.x      = element_text(size = 7, angle = 30, hjust = 1),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.title.y     = element_text(size = 7),
    legend.position  = "none",
    panel.grid       = element_blank(),
    panel.spacing.y  = unit(0.3, "mm"),
    plot.background  = element_rect(fill = "white", colour = NA),
    panel.background = element_rect(fill = "white", colour = NA)
  )

# ── Save: 115×70 mm ──────────────────────────────────────────────────────────
w_mm <- 70; h_mm <- 115
w <- w_mm / 25.4; h <- h_mm / 25.4
ragg::agg_png(file.path(OUT, "fig1F_stacked_violin_bcell_states.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p1F_v); dev.off()
grDevices::cairo_pdf(file.path(OUT, "fig1F_stacked_violin_bcell_states.pdf"),
                     width = w, height = h, family = "Arial")
print(p1F_v); dev.off()
log_msg("Saved: final/fig1/fig1F_stacked_violin_bcell_states.{png,pdf} [115×70mm]")
log_msg("=== fig1F_final.R complete ===")
