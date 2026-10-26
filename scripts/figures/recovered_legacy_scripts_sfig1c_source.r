#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_sfig1c_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/SFig1C_source.R
# ==========================================================
# sfig1C_final.R — 基于 scripts/revision/rasterize_supp_umaps.R（sfig1F: post-Harmony UMAP）修改
# 修改项（仅限本次新增需求，其余代码严格保持原脚本）：
#   1) 输出路径 → final/sfig1/
#   2) 画布 60×65 mm（原 140×115 mm）
#   3) 标题 9pt 居中（hjust=0.5）、无副标题；n 值按 REASONIX.md 要求并入标题
#      （60mm 页宽下 9pt 标题仅容纳约33字符，标题精简为 "UMAP post-Harmony (n = ...)"，
#        "colored by fraction" 由图例 Fraction 承载）
#   4) 坐标轴 7pt（axis.title=7，axis.text 保持原脚本 blank）
#   5) figure legend 7pt（原 legend.title 8 / legend.text 7.5）
#   6) 仅输出 PDF（用户要求）
# 保持原样：FRAC_COLORS 配色、fraction 着色列、pt/alpha 档位（102,117 → pt=0.1/alpha=0.70）、
#          ggrastr 光栅化、UMAP 1/UMAP 2 轴名
#
# Input:  results/01_scRNAseq/annotation/seurat_annotated.rds（102,117 细胞）
# Output: final/sfig1/sfig1C_umap_post_harmony.pdf
# 预计内存: ~20 GB（2.9GB Seurat 对象 + 光栅化）

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(ggrastr)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))

OUT <- "final/sfig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# ── 加载主 Seurat 对象 ────────────────────────────────────────────────────────
log_msg("Loading seurat_annotated.rds ...")
seu <- readRDS("results/01_scRNAseq/annotation/seurat_annotated.rds")
n <- ncol(seu)
log_msg(sprintf("  %d cells | reductions: %s", n, paste(Reductions(seu), collapse = ", ")))

# ── 细胞数自适应点参数（原脚本档位，同时满足 REASONIX.md 防蒙雾要求）──────────
if (n < 30000)      { pt <- 0.5; al <- 0.90
} else if (n < 100000) { pt <- 0.3; al <- 0.82
} else if (n < 300000) { pt <- 0.1; al <- 0.70
} else               { pt <- 0.05; al <- 0.55 }

# ── fraction 列与配色（原脚本）────────────────────────────────────────────────
frac_col <- intersect(c("fraction", "orig.ident", "sample"), colnames(seu@meta.data))[1]
log_msg(sprintf("  Fraction column: %s | levels: %s",
                frac_col, paste(unique(seu@meta.data[[frac_col]]), collapse = ", ")))

FRAC_COLORS <- c("CD45p" = "#41B6C4", "CD45n" = "#FC8D59",
                 "PBL"   = "#78C679", "CD45"   = "#41B6C4")

# ── post-Harmony UMAP（原脚本 p1F，仅字体/标题按新需求调整）──────────────────
p <- DimPlot(seu, reduction = "umap", group.by = frac_col,
             pt.size = pt, raster = FALSE) +
  scale_colour_manual(values = FRAC_COLORS, name = "Fraction",
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1))) +
# 标题：60mm 页宽下右侧图例会挤压标题可用宽度导致裁剪，故将图例移到底部横排
# （排版调整，图例内容/字号不变），标题保留完整信息（主诉+n值）
  labs(title = sprintf("UMAP post-Harmony (n = %s)", format(n, big.mark = ",")),
       x = "UMAP 1", y = "UMAP 2") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line      = element_line(linewidth = 0.35),
    axis.ticks     = element_blank(),
    axis.text      = element_blank(),
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
p$layers[[1]]$aes_params$alpha <- al

# 光栅化点云层（REASONIX.md 强制：PDF 点云必须 ggrastr 光栅化）
p_rast <- ggrastr::rasterize(p, layers = "Point", dpi = 600)
rm(seu); gc()

# ── 保存：60×65 mm，仅 PDF ────────────────────────────────────────────────────
w_mm <- 60; h_mm <- 65
grDevices::cairo_pdf(file.path(OUT, "sfig1C_umap_post_harmony.pdf"),
                     width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
print(p_rast); dev.off()
log_msg("Saved: final/sfig1/sfig1C_umap_post_harmony.pdf [60×65 mm, rasterized]")
log_msg("=== sfig1C_final.R complete ===")
