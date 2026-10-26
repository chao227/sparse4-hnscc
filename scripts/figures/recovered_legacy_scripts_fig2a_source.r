#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig2a_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig2A_source.R
# ==========================================================
# =============================================================================
# fig2A_wgcna_soft_threshold.R — Fig2A: WGCNA soft-threshold plot
# 原脚本: scripts/revision/fig2a_soft_threshold_hv.R（仅改输出尺寸/字体/仅PDF）
# 用户新增需求: 竖版单列（上下堆叠）50×75mm；标题9pt居中无副标题；
#              坐标轴7.5pt；图内文字≥7pt；只输出PDF
# 输入: results/02_WGCNA/vst_top5000.rds, results/02_WGCNA/network/wgcna_results.rds
# 输出: final/fig2/fig2A_wgcna_soft_threshold.pdf
# 预计内存: ~10 GB
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA); library(ggplot2); library(patchwork); library(dplyr)
})
options(stringsAsFactors = FALSE)
enableWGCNAThreads(nThreads = min(8, parallel::detectCores() - 2))

source(file.path("scripts", "theme_fig2.R"))

# ── 数据加载（与原脚本完全一致）────────────────────────────────────────────
log_msg("Loading VST matrix (5000 genes × 520 samples)...")
vst <- readRDS(file.path(BASE, "results/02_WGCNA/vst_top5000.rds"))
datExpr <- t(vst)
log_msg(sprintf("datExpr: %d samples × %d genes", nrow(datExpr), ncol(datExpr)))

wgcna_res <- readRDS(file.path(BASE, "results/02_WGCNA/network/wgcna_results.rds"))
best_power <- wgcna_res$best_power
log_msg(sprintf("Known best soft-threshold power: %d", best_power))

# ── pickSoftThreshold（与原脚本一致）───────────────────────────────────────
powers <- c(1:10, 12, 14, 16, 18, 20)
log_msg("Running pickSoftThreshold ...")
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 0,
                         networkType = "unsigned", corFnc = "cor")
sft_df <- sft$fitIndices
log_msg(sprintf("Best power from pickSoftThreshold: %s", sft$powerEstimate))

# ── 主题（新增需求: 标题9pt居中/无副标题/坐标轴7.5pt/最小≥7pt）─────────────
theme_sft <- function(bs = 7.5) {
  theme_classic(base_size = bs, base_family = "Arial") +
    theme(
      axis.line   = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks  = element_line(linewidth = 0.35, colour = "black"),
      axis.title  = element_text(size = 7.5),
      axis.text   = element_text(size = 7.5),
      panel.grid.major = element_line(colour = "grey93", linewidth = 0.2),
      panel.grid.minor = element_blank(),
      plot.background  = element_rect(fill = "white", colour = NA),
      panel.background = element_rect(fill = "white", colour = NA),
      plot.title  = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle = element_blank()
    )
}

# ── Panel 1: Scale-free topology R²（绘图逻辑与原脚本一致，仅字号≥7pt）─────
# 50mm 宽竖版下 15 个 power 标签物理放不下（原 183mm 横版可以），
# 仅标注关键 power 子集（c(1,3,5,8,12,16,20)）避免粘连
label_powers <- sft_df %>% filter(Power %in% c(1, 3, 5, 8, 12, 16, 20))
p_r2 <- ggplot(sft_df, aes(x = Power, y = SFT.R.sq)) +
  geom_hline(yintercept = 0.85, colour = "#CB181D", linewidth = 0.4, linetype = "dashed") +
  geom_vline(xintercept = best_power, colour = "#238B45", linewidth = 0.4, linetype = "dashed") +
  geom_line(colour = "#3182BD", linewidth = 0.5) +
  geom_point(colour = "#3182BD", size = 1.8, shape = 16) +
  geom_label(data = label_powers, aes(label = Power), size = pt2mm(7), family = "Arial",
             label.size = 0.15, label.padding = unit(1.2, "pt"),
             colour = "#252525", fill = "white") +
  annotate("text", x = best_power + 0.6, y = 0.05, label = sprintf("Power = %d", best_power),
           hjust = 0, size = pt2mm(7), family = "Arial", colour = "#238B45", fontface = "italic") +
  annotate("text", x = 1.0, y = 0.50,
           label = "R² threshold = 0.85", hjust = 0,
           size = pt2mm(7), family = "Arial", colour = "#CB181D", fontface = "italic") +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20),
                     expand = expansion(mult = 0.03)) +
  scale_y_continuous(breaks = c(0, 0.25, 0.50, 0.75, 1.00),
                     limits = c(0, 1.05), expand = expansion(mult = 0.02)) +
  labs(title = "Scale-free Topology Fit",
       x     = "Soft threshold (power)",
       y     = expression(R^2)) +
  theme_sft()

# ── Panel 2: Mean connectivity ──────────────────────────────────────────────
p_mk <- ggplot(sft_df, aes(x = Power, y = mean.k.)) +
  geom_vline(xintercept = best_power, colour = "#238B45", linewidth = 0.4, linetype = "dashed") +
  geom_line(colour = "#CB181D", linewidth = 0.5) +
  geom_point(colour = "#CB181D", size = 1.8, shape = 16) +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 20),
                     expand = expansion(mult = 0.03)) +
  scale_y_continuous(expand = expansion(mult = 0.05)) +
  labs(title = "Mean Connectivity",
       x     = "Soft threshold (power)",
       y     = "Mean connectivity") +
  theme_sft()

# ── 组合: 竖版单列（用户确认 50×75mm 用竖版）───────────────────────────────
p_v <- p_r2 / p_mk + plot_layout(nrow = 2, heights = c(1, 1))
save_pdf(p_v, "fig2A_wgcna_soft_threshold", 50, 75)

log_msg("=== fig2A complete ===")
