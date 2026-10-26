#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_sfig1d_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/SFig1D_source.R
# ==========================================================
# sfig1D_final.R — 基于 scripts/revision/fix_supp_figures_v2.R（FIX 4）修改
# 修改项（仅限本次新增需求，其余代码严格保持原脚本）：
#   1) 输出路径 → final/sfig1/
#   2) 画布 60×65 mm（原 110×95 mm）
#   3) 标题 9pt 居中（hjust=0.5）、无副标题；n 值按 REASONIX.md 要求并入标题
#      （60mm 页宽下右侧图例挤压标题导致裁剪，图例移到底部横排后标题保留完整信息；
#        原副标题 "5 transcriptional states" 与当前 rds 实际 3 类程序不符，按数据一致性规范弃用）
#   4) 坐标轴 7pt（axis.title=7, axis.text=7 保留刻度）
#   5) figure legend 7pt（原 legend.title 6.7 / legend.text 6.3）
#   6) 仅输出 PDF（用户要求）
# 保持原样：TUMOR_COLORS 配色、Program 图例、pt=0.5/alpha=0.85、ggrastr 光栅化、
#          UMAP 1/UMAP 2 轴名、图例位置 right
#
# Input:  results/subcluster/tumor_subclusters.rds（13,742 细胞）
# Output: final/sfig1/sfig1D_tumor_subtype_umap.pdf
# 预计内存: ~8 GB

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrastr)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))

OUT <- "final/sfig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── 加载肿瘤亚群对象 ──────────────────────────────────────────────────────────
log_msg("Loading tumor_subclusters.rds ...")
seu <- readRDS("results/subcluster/tumor_subclusters.rds")
n_tc <- ncol(seu)
log_msg(sprintf("  %d tumor cells | subtypes: %s",
                n_tc, paste(unique(seu$tumor_subtype), collapse = ", ")))

pt <- if (n_tc < 30000) 0.5 else 0.3   # 原脚本档位

# ── 配色（原脚本，未出现的类别自动不显示）────────────────────────────────────
TUMOR_COLORS <- c(
  "Differentiated" = "#E28E2C", "Basal" = "#D24B40",
  "pEMT" = "#9E7BB5", "Cycling" = "#74C476", "Stress" = "#BDBDBD"
)

# 标题：60mm 页宽下右侧图例会挤压标题可用宽度导致裁剪，故将图例移到底部横排
# （排版调整，图例内容/字号不变），标题保留完整信息（主诉+n值）；
# 原副标题 "5 transcriptional states" 与当前 rds 实际 3 类程序不符，按数据一致性规范弃用
p <- DimPlot(seu, reduction = "umap", group.by = "tumor_subtype",
             cols = TUMOR_COLORS, pt.size = pt, raster = FALSE) +
  labs(title = sprintf("Tumor cell programs (n = %s)", format(n_tc, big.mark = ",")),
       x = "UMAP 1", y = "UMAP 2") +
  guides(colour = guide_legend(override.aes = list(size = 3), title = "Program", nrow = 2)) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line      = element_line(linewidth = 0.35),
    axis.ticks     = element_line(linewidth = 0.35),
    axis.text      = element_text(size = 7),
    axis.title     = element_text(size = 7),
    legend.title   = element_text(size = 7),
    legend.text    = element_text(size = 7),
    legend.key.size = unit(2, "mm"),
    legend.spacing.x = unit(0.5, "mm"),
    legend.spacing  = unit(2, "pt"),
    legend.margin  = margin(0, 0, 0, 0),
    plot.title     = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle  = element_blank(),
    panel.grid     = element_blank(),
    legend.position = "bottom",
    legend.direction = "horizontal"
  )
p$layers[[1]]$aes_params$alpha <- 0.85   # 原脚本 alpha

# 光栅化点云层（REASONIX.md 强制）
p_rast <- ggrastr::rasterize(p, layers = "Point", dpi = 600)
rm(seu); gc()

# ── 保存：60×65 mm，仅 PDF ────────────────────────────────────────────────────
w_mm <- 60; h_mm <- 65
grDevices::cairo_pdf(file.path(OUT, "sfig1D_tumor_subtype_umap.pdf"),
                     width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
print(p_rast); dev.off()
log_msg("Saved: final/sfig1/sfig1D_tumor_subtype_umap.pdf [60×65 mm, rasterized]")
log_msg("=== sfig1D_final.R complete ===")
