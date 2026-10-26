#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig2b_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig2B_source.R
# ==========================================================
# =============================================================================
# fig2B_wgcna_dendrogram.R — Fig2B: Gene clustering dendrogram + module colors
# 原脚本: scripts/module02/02c_wgcna.R Step 4（plotDendroAndColors）
# 用户新增需求: 140×75mm；标题9pt居中；最小字号≥7pt；只输出PDF
# 输入: results/02_WGCNA/network/wgcna_results.rds
# 输出: final/fig2/fig2B_wgcna_dendrogram.pdf
# 预计内存: ~5 GB
# =============================================================================

suppressPackageStartupMessages({ library(WGCNA) })
options(stringsAsFactors = FALSE)

source(file.path("scripts", "theme_fig2.R"))

log_msg("Loading WGCNA results...")
wg <- readRDS(file.path(BASE, "results/02_WGCNA/network/wgcna_results.rds"))

net          <- wg$net
module_colors <- net$colors
dendro        <- net$dendrograms[[1]]
block_genes   <- net$blockGenes[[1]]

# hclust 对象保存时可能丢失 labels，用 gene_names 补全（仅用于内部索引，绘图 dendroLabels=FALSE）
if (length(dendro$labels) == 0 && !is.null(wg$gene_names)) {
  dendro$labels <- wg$gene_names[block_genes]
}
cat(sprintf("dendrogram leaves: %d, colors: %d\n", length(dendro$labels), length(module_colors[block_genes])))

# ── 绘图（逻辑与原脚本一致；字号按新规范: 标题9pt=cex.main 0.75）───────────
w_mm <- 140; h_mm <- 75
f <- file.path(OUT, "fig2B_wgcna_dendrogram.pdf")
grDevices::cairo_pdf(f, width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
par(mar = c(1.5, 4, 2, 1), cex.main = cex_pt(9), cex.lab = cex_pt(7.5),
    cex.axis = cex_pt(7.5), family = "Arial")
plotDendroAndColors(
  dendro,
  module_colors[block_genes],
  "Module colors",
  dendroLabels = FALSE, hang = 0.03,
  addGuide = TRUE, guideHang = 0.05,
  main = "Gene Clustering Dendrogram and Module Colors"
)
dev.off()
log_msg(sprintf("Saved: %s (%.0f x %.0f mm)", f, w_mm, h_mm))
log_msg("=== fig2B complete ===")
