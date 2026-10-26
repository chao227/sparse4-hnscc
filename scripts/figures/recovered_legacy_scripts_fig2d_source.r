#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig2d_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig2D_source.R
# ==========================================================
# =============================================================================
# fig2D_yellow_hub_network.R — Fig2D: Top 13 yellow-module hub genes (kME bars)
# 原脚本: scripts/manuscript_fix/fix_fig2D_v2.R（仅改输出尺寸/字体/仅PDF）
# 用户新增需求: 70×50mm；标题9pt居中无副标题；坐标轴7.5pt；基因名7.5pt；
#              图内数值≥7pt；只输出PDF
# 输入: results/02_WGCNA/hub_genes/hub_genes_all.csv
# 输出: final/fig2/fig2D_yellow_hub_network.pdf
# 预计内存: <1 GB
# =============================================================================

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr)
})

source(file.path("scripts", "theme_fig2.R"))

log_msg("Building Fig2D hub-gene kME bar chart...")

hub <- read.csv(file.path(BASE, "results/02_WGCNA/hub_genes/hub_genes_all.csv")) %>%
  filter(module == "yellow") %>%
  arrange(desc(kME)) %>%
  head(13) %>%
  mutate(
    gene    = factor(gene, levels = rev(gene)),
    is_hub  = gene == "CD79A",
    kME_pct = (kME - min(kME)) / (max(kME) - min(kME))
  )

kmin <- min(hub$kME)
kmax <- max(hub$kME)

pal_gradient <- colorRampPalette(c("#A8DADC","#457B9D","#1D3557"))(100)

# 绘图逻辑与原脚本一致；主题/字号按新规范
p2d <- ggplot(hub, aes(x = kME, y = gene)) +
  geom_col(aes(x = kmax + 0.003), fill = "#F0F4F8", width = 0.65) +
  geom_col(aes(fill = kME_pct), width = 0.65, color = NA) +
  geom_col(data = filter(hub, is_hub),
           aes(x = kME), fill = "#E63946", width = 0.65, color = NA) +
  geom_text(aes(x = kME - 0.0004, label = sprintf("%.3f", kME)),
            hjust = 1, size = pt2mm(7), fontface = "bold", color = "white") +
  scale_fill_gradientn(colors = pal_gradient, guide = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.02)),
    breaks = seq(0.86, 0.94, by = 0.01),
    labels = function(x) sprintf("%.2f", x)
  ) +
  # coord_cartesian zooms without censoring bars that start at x=0
  # 下限 0.86（原 0.885-0.004）：避免条内白色数值标签左端与 y 轴基因名重叠
  coord_cartesian(xlim = c(0.86, kmax + 0.006)) +
  labs(
    title    = "Top 13 Yellow Module Hub Genes",
    x = "Module Membership (kME)", y = NULL
  ) +
  theme_fig2() +
  theme(
    axis.text.y = element_text(size = 7.5, face = "italic",
                               color = ifelse(levels(hub$gene) == "CD79A",
                                              "#E63946", "grey25")),
    axis.text.x = element_text(size = 7.5, color = "grey40"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#E8ECF0", linewidth = 0.5),
    axis.ticks         = element_blank()
  )

save_pdf(p2d, "fig2D_yellow_hub_network", 70, 50)
log_msg("=== fig2D complete ===")
