#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_sfig1a_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/SFig1A_source.R
# ==========================================================
# sfig1A_final.R — 基于 scripts/module01/01a_qc_per_sample.R（Step 5 的 p_qc）修改
# 修改项（仅限本次新增需求，其余代码严格保持原脚本）：
#   1) 输出路径 → final/sfig1/
#   2) 画布 195×65 mm（原 16×5 in）
#   3) 标题 9pt 居中（hjust=0.5）、无副标题；n 值按 REASONIX.md 要求并入标题
#   4) 坐标轴 7-8pt（axis.title=7.5, axis.text=7）
#   5) figure legend 7pt
#   6) x 轴 33 个样本标签改为 90° 垂直 7pt（原 45° 在 195mm 宽下必然重叠，用户已确认）
#   7) 仅输出 PDF（用户要求）
# 保持原样：violin+boxplot 几何、CD45p/CD45n 配色、标题文字、x/y 轴名
#
# Input:  results/01_scRNAseq/QC/qc_meta_all_cells.csv（01a Step6 保存的 all_meta，102,117 细胞）
# Output: final/sfig1/sfig1A_qc_gene_count.pdf
# 预计内存: <2 GB

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))

OUT <- "final/sfig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── 数据加载（与原脚本 all_meta 完全一致的数据源）────────────────────────────
log_msg("Loading qc_meta_all_cells.csv ...")
all_meta <- read.csv("results/01_scRNAseq/QC/qc_meta_all_cells.csv",
                     row.names = 1, stringsAsFactors = FALSE)
log_msg(sprintf("  %d cells | %d samples | fractions: %s",
                nrow(all_meta), length(unique(all_meta$sample)),
                paste(unique(all_meta$fraction), collapse = ", ")))

# ── 绘图（原脚本 p_qc 代码，仅字体/角度/标题按新需求调整）────────────────────
p_qc <- ggplot(all_meta, aes(x = sample, y = nFeature_RNA, fill = fraction)) +
  geom_violin(scale = "width") +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  theme_bw(base_size = 7, base_family = "Arial") +
  theme(
    axis.text.x    = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
    axis.text.y    = element_text(size = 7),
    axis.title     = element_text(size = 7.5),
    legend.title   = element_text(size = 7),
    legend.text    = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.title     = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle  = element_blank(),
    panel.grid     = element_blank()
  ) +
  labs(title = sprintf("Genes per Cell after QC (n = %s cells)",
                       format(nrow(all_meta), big.mark = ",")),
       x = "Sample", y = "Number of Genes") +
  scale_fill_manual(values = c("CD45p" = "#E41A1C", "CD45n" = "#377EB8"))

# ── 保存：195×65 mm，仅 PDF ───────────────────────────────────────────────────
w_mm <- 195; h_mm <- 65
grDevices::cairo_pdf(file.path(OUT, "sfig1A_qc_gene_count.pdf"),
                     width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
print(p_qc); dev.off()
log_msg("Saved: final/sfig1/sfig1A_qc_gene_count.pdf [195×65 mm]")
log_msg("=== sfig1A_final.R complete ===")
