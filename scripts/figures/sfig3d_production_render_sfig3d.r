#!/usr/bin/env Rscript

# ==========================================================
# Script   : sfig3d_production_render_sfig3d.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_repairs/SFig3D_production/scripts/render_SFig3D.R
# ==========================================================

suppressPackageStartupMessages({
  library(ComplexHeatmap)
  library(circlize)
  library(grid)
  library(ragg)
})

options(stringsAsFactors = FALSE, bitmapType = "cairo", warn = 1)

args <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args, value = TRUE)
stopifnot(length(script_arg) == 1L)
script_file <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
pkg_dir <- normalizePath(file.path(dirname(script_file), ".."), mustWork = TRUE)
input_dir <- file.path(pkg_dir, "inputs")
source_dir <- file.path(pkg_dir, "source_data")
output_dir <- file.path(pkg_dir, "outputs")
log_dir <- file.path(pkg_dir, "logs")
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

log_con <- file(file.path(log_dir, "render_SFig3D.log"), open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")
on.exit({
  while (sink.number(type = "message") > 0) sink(type = "message")
  while (sink.number(type = "output") > 0) sink(type = "output")
  close(log_con)
}, add = TRUE)

cat(sprintf("started_at\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
cat(sprintf("package_dir\t%s\n", pkg_dir))

vst <- readRDS(file.path(input_dir, "vst_top25mad.rds"))
master <- read.csv(file.path(input_dir, "master_survival_IES_table.csv"),
                   check.names = FALSE)
trait <- read.csv(file.path(input_dir, "trait_matrix.csv"),
                  check.names = FALSE)

required_genes <- c("CD79A", "IKZF1")
stopifnot(all(required_genes %in% rownames(vst)))
stopifnot("sample_id" %in% names(master), "sample_id" %in% names(trait))
stopifnot(!anyDuplicated(master$sample_id))

vst_id15 <- substr(colnames(vst), 1L, 15L)
trait_id15 <- substr(trait$sample_id, 1L, 15L)
vst_idx <- match(master$sample_id, vst_id15)
trait_idx <- match(master$sample_id, trait_id15)
stopifnot(!anyNA(vst_idx), !anyNA(trait_idx))
stopifnot(!anyDuplicated(vst_id15[vst_idx]), !anyDuplicated(trait_id15[trait_idx]))

score_cols <- c(
  "B_cell_naive", "B_cell_memory", "Plasma_cell", "T_CD8_cytotoxic",
  "T_CD4_naive", "T_CD4_memory_resting", "T_regulatory", "NK_resting",
  "NK_activated", "Monocyte", "Macrophage_M1", "Macrophage_M2",
  "DC_resting", "DC_activated", "Mast_cell_resting", "Neutrophil",
  "NK_T_cell"
)
stopifnot(length(score_cols) == 17L, all(score_cols %in% names(trait)))

cd79a <- as.numeric(vst["CD79A", vst_idx])
ikzf1 <- as.numeric(vst["IKZF1", vst_idx])
stopifnot(all(is.finite(cd79a)), all(is.finite(ikzf1)))
med_cd79a <- median(cd79a)
med_ikzf1 <- median(ikzf1)

quad_chr <- ifelse(cd79a >= med_cd79a & ikzf1 >= med_ikzf1, "LymphHot",
            ifelse(cd79a >= med_cd79a & ikzf1 < med_ikzf1, "B-skewed",
            ifelse(cd79a < med_cd79a & ikzf1 >= med_ikzf1, "T-skewed", "Cold")))
quad <- factor(quad_chr, levels = c("LymphHot", "B-skewed", "T-skewed", "Cold"))
expected_counts <- c(LymphHot = 171L, `B-skewed` = 67L,
                     `T-skewed` = 67L, Cold = 171L)
observed_counts <- table(quad)
stopifnot(identical(as.integer(observed_counts), as.integer(expected_counts)))

ssg <- trait[trait_idx, score_cols, drop = FALSE]
stopifnot(nrow(ssg) == 476L)
stopifnot(all(vapply(ssg, is.numeric, logical(1))))
stopifnot(all(is.finite(as.matrix(ssg))))

membership <- data.frame(
  sample_id = master$sample_id,
  source_vst_sample_id = colnames(vst)[vst_idx],
  source_trait_sample_id = trait$sample_id[trait_idx],
  CD79A = cd79a,
  IKZF1 = ikzf1,
  CD79A_median = med_cd79a,
  IKZF1_median = med_ikzf1,
  quadrant = as.character(quad),
  stringsAsFactors = FALSE
)
write.table(membership, file.path(source_dir, "SFig3D_patient_mapping_and_quadrants.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

kw_p <- vapply(score_cols, function(sc) {
  stats::kruskal.test(ssg[[sc]] ~ quad)$p.value
}, numeric(1))
kw_df <- data.frame(score = score_cols, p = unname(kw_p), stringsAsFactors = FALSE)
kw_df <- kw_df[order(kw_df$p, kw_df$score), , drop = FALSE]
kw_df$BH_q <- p.adjust(kw_df$p, method = "BH")
kw_df$nominal_p_rank <- seq_len(nrow(kw_df))
kw_df$selected_top15 <- kw_df$nominal_p_rank <= 15L
write.table(kw_df, file.path(source_dir, "SFig3D_Kruskal_Wallis_BH17.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
top15 <- kw_df$score[kw_df$selected_top15]
stopifnot(length(top15) == 15L, all(kw_df$BH_q[kw_df$selected_top15] < 0.05))

mean_by_quad <- sapply(levels(quad), function(q) {
  colMeans(ssg[quad == q, score_cols, drop = FALSE])
})
mat_mean <- t(mean_by_quad)
mat_scaled <- t(scale(t(mat_mean)))
stopifnot(all(is.finite(mat_scaled)))
mat_plot <- mat_scaled[, top15, drop = FALSE]

write.table(data.frame(quadrant = rownames(mat_mean), mat_mean, check.names = FALSE),
            file.path(source_dir, "SFig3D_group_mean_ssGSEA_all17.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(data.frame(quadrant = rownames(mat_scaled), mat_scaled, check.names = FALSE),
            file.path(source_dir, "SFig3D_row_zscore_all17.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(data.frame(quadrant = rownames(mat_plot), mat_plot, check.names = FALSE),
            file.path(source_dir, "SFig3D_plotted_matrix_top15.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

quad_cols <- c(
  LymphHot = "#D73027", `B-skewed` = "#FC8D59",
  `T-skewed` = "#74ADD1", Cold = "#4575B4"
)
col_fun <- colorRamp2(c(-2, -1, 0, 1, 2),
                      c("#4575B4", "#74ADD1", "white", "#FC8D59", "#D73027"))

text_contract <- data.frame(
  element = c("column_title", "row_names", "column_names", "annotation_name",
              "heatmap_legend_title", "heatmap_legend_labels"),
  font_family = "Arial",
  font_size_pt = c(9, 8, 7, 7.5, 7, 7),
  stringsAsFactors = FALSE
)
stopifnot(min(text_contract$font_size_pt) >= 7)
write.table(text_contract, file.path(source_dir, "SFig3D_typography_contract.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

make_heatmap <- function() {
  Heatmap(
    mat_plot,
    name = "z-score",
    col = col_fun,
    row_names_gp = gpar(fontsize = 8, fontfamily = "Arial"),
    column_names_gp = gpar(fontsize = 7, fontfamily = "Arial"),
    column_names_rot = 45,
    column_title = "ssGSEA immune cell scores across\nLymphHot quadrants (top 15; 15 FDR<0.05)",
    column_title_gp = gpar(fontsize = 9, fontface = "bold", fontfamily = "Arial"),
    cluster_rows = FALSE,
    cluster_columns = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "complete",
    border = TRUE,
    rect_gp = gpar(col = "white", lwd = 0.5),
    left_annotation = rowAnnotation(
      Subtype = anno_simple(
        rownames(mat_plot),
        col = quad_cols[rownames(mat_plot)],
        border = FALSE
      ),
      annotation_name_side = "top",
      annotation_name_gp = gpar(fontsize = 7.5, fontfamily = "Arial")
    ),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 7, fontfamily = "Arial"),
      labels_gp = gpar(fontsize = 7, fontfamily = "Arial")
    ),
    width = unit(6.4, "cm"),
    height = unit(3.7, "cm")
  )
}

draw_heatmap <- function(ht) {
  withCallingHandlers(
    draw(
      ht,
      heatmap_legend_side = "right",
      annotation_legend_side = "right",
      padding = unit(c(2, 2, 2, 2), "mm")
    ),
    warning = function(w) {
      if (grepl("PostScript.*Arial", conditionMessage(w))) {
        invokeRestart("muffleWarning")
      }
    }
  )
}

w_mm <- 100
h_mm <- 78
pdf_file <- file.path(output_dir, "SFig3D_ssGSEA_quadrant_heatmap.pdf")
png_file <- file.path(output_dir, "SFig3D_ssGSEA_quadrant_heatmap.png")

grDevices::cairo_pdf(pdf_file, width = w_mm / 25.4, height = h_mm / 25.4,
                     family = "Arial", onefile = FALSE)
ht_pdf <- draw_heatmap(make_heatmap())
pdf_order <- column_order(ht_pdf)
dev.off()

ragg::agg_png(png_file, width = w_mm / 25.4, height = h_mm / 25.4,
              units = "in", res = 600, background = "white", scaling = 1)
ht_png <- draw_heatmap(make_heatmap())
png_order <- column_order(ht_png)
dev.off()

stopifnot(identical(as.integer(pdf_order), as.integer(png_order)))
plot_order <- data.frame(
  position = seq_along(pdf_order),
  source_column_index = as.integer(pdf_order),
  score = colnames(mat_plot)[as.integer(pdf_order)],
  stringsAsFactors = FALSE
)
write.table(plot_order, file.path(source_dir, "SFig3D_rendered_column_order.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

geometry <- data.frame(
  format = c("PDF", "PNG"),
  width_mm = w_mm,
  height_mm = h_mm,
  dpi = c(NA, 600),
  stringsAsFactors = FALSE
)
write.table(geometry, file.path(source_dir, "SFig3D_export_geometry.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

writeLines(capture.output(sessionInfo()), file.path(log_dir, "sessionInfo.txt"))
cat(sprintf("quadrant_counts\t%s\n", paste(as.integer(observed_counts), collapse = ",")))
cat(sprintf("medians\tCD79A=%.17g;IKZF1=%.17g\n", med_cd79a, med_ikzf1))
cat(sprintf("BH_lt_0.05_all17\t%d\n", sum(kw_df$BH_q < 0.05)))
cat(sprintf("BH_lt_0.05_plotted15\t%d\n",
            sum(kw_df$BH_q[kw_df$selected_top15] < 0.05)))
cat(sprintf("top15\t%s\n", paste(top15, collapse = ",")))
cat(sprintf("column_order\t%s\n", paste(plot_order$score, collapse = ",")))
cat(sprintf("pdf\t%s\n", pdf_file))
cat(sprintf("png\t%s\n", png_file))
cat(sprintf("completed_at\t%s\n", format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z")))
