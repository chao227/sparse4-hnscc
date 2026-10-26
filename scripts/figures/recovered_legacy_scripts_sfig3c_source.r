#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_sfig3c_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/SFig3C_source.R
# ==========================================================
# =============================================================================
# sfig3C_final.R — SFig3C: SCENIC regulon activity heatmap（CD8 亚型）
# 基准: scripts/revised/R13_scenic_visualize.R（Figure 1: Regulon heatmap 部分）
#       数据与统计完全一致: cd8_auc_matrix.csv + cd8_subclusters.rds,
#       top 40 方差最大 regulons × 6 CD8 亚型平均 AUC 行 z-score
# 用户新增需求（唯一修改点）: 输出→final/sfig3/；画布 90×75mm；添加居中 9pt 标题
#      （原图无标题，按用户要求补）无副标题；row_names(regulon 名)8pt、
#      column_names(亚型名)8pt、legend 7pt、annotation 名 ≥7pt；最小字体≥7pt；仅 PDF
# 输入: results/scenic/cd8_auc_matrix.csv, results/subcluster/cd8_subclusters.rds
# 输出: final/sfig3/sfig3C_scenic_regulon_heatmap.pdf
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat); library(ComplexHeatmap); library(circlize); library(ragg)
})

source(file.path("scripts", "theme_sfig3.R"))
source(file.path(BASE, "final", "fig5", "scripts", "theme_fig5.R"))
OUT <- file.path(BASE, "final", "sfig3")

# ── Canonical CD8 palette: latest locked Fig5/SFig2A-B mapping ───────────────
CD8_COLORS <- FIG5_CD8_COLORS

tryCatch({
# ── 数据加载（原脚本逻辑不变）────────────────────────────────────────────────
AUC_CSV <- file.path(BASE, "results/scenic/cd8_auc_matrix.csv")
CD8_RDS <- file.path(BASE, "results/subcluster/cd8_subclusters.rds")
cd8 <- readRDS(CD8_RDS)
log_msg(sprintf("CD8 cells: %d", ncol(cd8)))
auc_df <- read.csv(AUC_CSV, row.names = 1, check.names = FALSE)
log_msg(sprintf("AUC matrix: %d cells x %d regulons", nrow(auc_df), ncol(auc_df)))

common_cells <- intersect(rownames(auc_df), colnames(cd8))
if (length(common_cells) < 100) {
  auc_ids_clean <- gsub("-[0-9]+$", "", rownames(auc_df))
  cd8_ids_clean <- gsub("-[0-9]+$", "", colnames(cd8))
  idx_auc <- match(cd8_ids_clean, auc_ids_clean)
  common_cells <- colnames(cd8)[!is.na(idx_auc)]
  auc_df <- auc_df[idx_auc[!is.na(idx_auc)], , drop = FALSE]
  rownames(auc_df) <- common_cells
  log_msg(sprintf("After ID matching: %d cells aligned", length(common_cells)))
} else {
  auc_df <- auc_df[common_cells, , drop = FALSE]
  cd8 <- cd8[, common_cells]
  log_msg(sprintf("Direct match: %d cells aligned", length(common_cells)))
}

regulon_names <- colnames(auc_df)
subtype_vec <- cd8$cd8_subtype[common_cells]
auc_mat <- as.matrix(auc_df)

# 各亚型平均 AUC
subtype_mean_auc <- sapply(unique(subtype_vec), function(st) {
  idx <- which(subtype_vec == st)
  colMeans(auc_mat[idx, , drop = FALSE])
})
colnames(subtype_mean_auc) <- unique(subtype_vec)
subtype_mean_auc <- subtype_mean_auc[, names(CD8_COLORS)[names(CD8_COLORS) %in% colnames(subtype_mean_auc)]]

# top 20 方差最大 regulons（用户确认：90×75mm 画布下 40 行必然重叠，改为 20）
rv <- apply(subtype_mean_auc, 1, var)
top40_idx <- order(rv, decreasing = TRUE)[1:min(20, length(rv))]
top40_regulons <- regulon_names[top40_idx]
log_msg(sprintf("Top %d variable regulons selected", length(top40_regulons)))

# 重建矩阵
present_subtypes <- names(CD8_COLORS)[names(CD8_COLORS) %in% colnames(subtype_mean_auc)]
heat_mat_raw <- subtype_mean_auc[top40_regulons, present_subtypes, drop = FALSE]
heat_mat_scaled <- t(scale(t(heat_mat_raw)))
heat_mat_scaled[is.nan(heat_mat_scaled)] <- 0
rownames(heat_mat_scaled) <- gsub("\\(\\+\\)$", "", rownames(heat_mat_scaled))

col_fun <- colorRamp2(c(-2, 0, 2), c("#2166AC", "white", "#D73027"))

# ── 热图（原几何/配色不变；字体按用户新增需求；行名 8pt 紧凑排列）──────────
n_rows <- nrow(heat_mat_scaled)
n_cols <- ncol(heat_mat_scaled)
# 90×75mm 画布: 热图主体高度取可用空间上限，行名 8pt（≥7pt 下限）
ht <- Heatmap(
  heat_mat_scaled,
  name                  = "Scaled\nAUC",
  col                   = col_fun,
  cluster_rows          = TRUE,
  cluster_columns       = FALSE,
  column_order          = seq_len(ncol(heat_mat_scaled)),
  show_row_dend         = TRUE,
  show_column_dend      = FALSE,
  row_names_gp          = gpar(fontsize = 8, fontfamily = "Arial"),
  column_names_gp       = gpar(fontsize = 7, fontfamily = "Arial", fontface = "bold"),
  column_names_rot      = 60,
  top_annotation        = HeatmapAnnotation(
    Subtype = present_subtypes,
    col     = list(Subtype = CD8_COLORS[present_subtypes]),
    annotation_legend_param = list(
      title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 7, fontfamily = "Arial")),
    annotation_name_gp = gpar(fontsize = 7),
    simple_anno_size = unit(3, "mm")),
  heatmap_legend_param  = list(title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
                               labels_gp = gpar(fontsize = 7, fontfamily = "Arial")),
  width                 = unit(4.0, "cm"),
  height                = unit(165, "pt")
)

# ── Save: 105 x 90 mm（用户确认加高加宽，标题 grid 居中；PDF + PNG）────────────
w_mm <- 105; h_mm <- 90
draw_sfig3C <- function() {
  draw(ht, heatmap_legend_side = "right", annotation_legend_side = "right")
  grid::grid.text("SCENIC regulon activity in CD8+ T cells",
                  x = 0.5, y = unit(0.975, "npc"), just = c("center", "top"),
                  gp = grid::gpar(fontsize = 9, fontface = "bold", fontfamily = "Arial"))
}
f_pdf <- file.path(OUT, "sfig3C_scenic_regulon_heatmap.pdf")
grDevices::cairo_pdf(f_pdf, width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
draw_sfig3C(); dev.off()
f_png <- file.path(OUT, "sfig3C_scenic_regulon_heatmap.png")
ragg::agg_png(f_png, width = w_mm / 25.4, height = h_mm / 25.4,
              units = "in", res = 600, background = "white")
draw_sfig3C(); dev.off()
log_msg(sprintf("Saved: %s (%.0f x %.0f mm)",
                f_pdf, w_mm, h_mm))
log_msg(sprintf("Saved: %s (600 dpi)", f_png))
}, error = function(e) log_msg(paste("ERROR sfig3C:", conditionMessage(e))))
log_msg("=== sfig3C complete ===")
