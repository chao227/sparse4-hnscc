#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_sfig4bcd_min7pt.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_sfig4bcd_min7pt.R
# ==========================================================

# Presentation-only successor for SFig.4B-D. Existing barcode-matched raster
# tiles and frozen source summaries are retained; only panel typography and
# layout are revised. All visible semantic text is at least 7.2 pt.

suppressPackageStartupMessages({
  library(png)
  library(ragg)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
cache <- file.path(root, "replacement_spatial_figures_20261019/source_data/render_qc/tile_cache")
source_root <- file.path(root, "replacement_spatial_figures_20261019/source_data")
out <- file.path(root, "yellow_module_spatial_20261019/render/figures/SFig4")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(file.path(source_root, "GSE208253_sample_manifest.tsv"), stringsAsFactors = FALSE)
detection <- read.delim(file.path(source_root, "GSE208253_IKZF1_direct_spot_detection.tsv"), stringsAsFactors = FALSE)
rctd_limits <- read.delim(file.path(source_root, "render_qc/RCTD_display_limits.tsv"), stringsAsFactors = FALSE)
global_limits <- read.delim(file.path(source_root, "render_qc/global_display_limits.tsv"), stringsAsFactors = FALSE)

picks <- c("GSM6339638", "GSM6339633", "GSM6339632")
pick_rows <- manifest[match(picks, manifest$GSM), , drop = FALSE]
stopifnot(!anyNA(pick_rows$GSM), nrow(detection) == 12L)

read_tile <- function(gsm, suffix) {
  path <- file.path(cache, paste0(gsm, "_", suffix, ".png"))
  stopifnot(file.exists(path))
  readPNG(path)
}
draw_tile <- function(gsm, suffix, x0, y0, x1, y1) {
  rasterImage(read_tile(gsm, suffix), x0, y0, x1, y1, interpolate = FALSE)
  rect(x0, y0, x1, y1, border = "grey55", lwd = 0.35)
}
draw_key <- function(x0, x1, y0, y1, limits, colour) {
  pal <- colorRampPalette(c("#F7FBFF", colour))(100)
  for (i in seq_along(pal)) rect(x0 + (i - 1) / 100 * (x1 - x0), y0,
                                  x0 + i / 100 * (x1 - x0), y1,
                                  col = pal[i], border = NA)
  text(x0, y0 - 0.012, formatC(limits[1], format = "f", digits = 2), adj = c(0, 1), cex = 0.65)
  text(x1, y0 - 0.012, formatC(limits[2], format = "f", digits = 2), adj = c(1, 1), cex = 0.65)
}
save_pair <- function(fun, stem, width_mm, height_mm) {
  ragg::agg_png(paste0(stem, ".png"), width = width_mm / 25.4, height = height_mm / 25.4,
                units = "in", res = 600, background = "white")
  fun(); dev.off()
  cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
            family = "Arial", onefile = TRUE)
  fun(); dev.off()
}

cell_suffix <- c("B_cell", "CD8_T_Cytotoxic", "CD4_T_Helper", "Tumor_cell")
cell_labels <- c("B cell", "CD8 T cytotoxic", "CD4 T helper", "Tumour cell")
cell_colours <- c("#3182BD", "#756BB1", "#41AB5D", "#E6550D")

draw_b <- function() {
  par(family = "Arial", mar = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.975, "Representative RCTD-based immune-cell spatial distributions", font = 2, cex = 0.82)
  xs <- seq(0.13, 0.985, length.out = 6)
  headers <- c("H&E", cell_labels)
  for (j in 1:5) text(mean(xs[j:(j + 1)]), 0.91, headers[j], font = 2, cex = 0.65)
  yy <- c(0.79, 0.54, 0.29)
  for (i in seq_along(picks)) {
    gsm <- picks[i]
    draw_tile(gsm, "he", xs[1], yy[i] - 0.15, xs[2] - 0.005, yy[i])
    for (j in seq_along(cell_suffix)) {
      draw_tile(gsm, cell_suffix[j], xs[j + 1], yy[i] - 0.15, xs[j + 2] - 0.005, yy[i])
    }
    label <- paste0(pick_rows$patient_id[i], "-", pick_rows$section_id[i])
    text(0.012, yy[i] - 0.075, label, font = 2, cex = 0.65, adj = 0)
  }
  for (j in seq_along(cell_suffix)) {
    limits <- as.numeric(rctd_limits[rctd_limits$feature == cell_suffix[j], c("lower", "upper")])
    draw_key(xs[j + 1] + 0.012, xs[j + 2] - 0.025, 0.050, 0.058, limits, cell_colours[j])
  }
}

detection$label <- paste0(manifest$patient_id[match(detection$GSM, manifest$GSM)], "-", detection$section_id)
ik_limits <- as.numeric(global_limits[global_limits$feature == "IKZF1 direct log1p expression", c("lower", "upper")])
draw_c <- function() {
  par(family = "Arial", mar = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.978, "IKZF1 expression distribution across spatial sections", font = 2, cex = 0.80)
  par(fig = c(0.055, 0.46, 0.20, 0.86), new = TRUE, mar = c(3.3, 4.3, 0.4, 0.3), xpd = FALSE, family = "Arial")
  plot(100 * detection$detection_rate, seq_len(nrow(detection)), type = "n", xlim = c(0, 28), ylim = c(0.5, 12.5),
       yaxt = "n", xlab = "IKZF1 detection fraction (%)", ylab = "Patient-section", cex.axis = 0.65, cex.lab = 0.65)
  axis(2, at = 1:12, labels = detection$label, las = 1, cex.axis = 0.65)
  abline(v = 5, lty = 2, col = "#C26D3A", lwd = 0.6)
  points(100 * detection$detection_rate, 1:12, pch = 21, bg = "#5B83A6", col = "white", cex = 0.72)
  par(fig = c(0, 1, 0, 1), new = TRUE, mar = c(0, 0, 0, 0), xpd = NA, family = "Arial")
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.25, 0.885, "2.2-24.8% detected\n4/12 sections below 5%", font = 2, cex = 0.65, col = "#A65A2C")
  xs <- c(0.52, 0.755, 0.99)
  text(c(0.637, 0.872), 0.90, c("H&E", "IKZF1 expression"), font = 2, cex = 0.65)
  yy <- c(0.79, 0.54, 0.29)
  for (i in seq_along(picks)) {
    draw_tile(picks[i], "he", xs[1], yy[i] - 0.15, xs[2] - 0.006, yy[i])
    draw_tile(picks[i], "ik", xs[2], yy[i] - 0.15, xs[3], yy[i])
    text(xs[1] - 0.007, yy[i] - 0.075, paste0(pick_rows$patient_id[i], "-", pick_rows$section_id[i]),
         font = 2, cex = 0.65, adj = 1)
  }
  draw_key(0.77, 0.95, 0.050, 0.058, ik_limits, "#FB6A4A")
}

draw_d <- function() {
  par(family = "Arial", mar = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.96, "Spatial image-coordinate alignment QC", font = 2, cex = 0.88)
  qc_gsm <- c("GSM6339632", "GSM6339638")
  for (j in 1:2) {
    x0 <- if (j == 1) 0.04 else 0.52
    x1 <- if (j == 1) 0.48 else 0.96
    draw_tile(qc_gsm[j], "raw", x0, 0.17, x1, 0.86)
    row <- manifest[manifest$GSM == qc_gsm[j], , drop = FALSE]
    text(mean(c(x0, x1)), 0.11, paste0(row$patient_id, "-", row$section_id), font = 2, cex = 0.65)
  }
}

save_pair(draw_b, file.path(out, "SFig4B_RCTD_immune_cell_spatial_distribution"), 130, 85)
save_pair(draw_c, file.path(out, "SFig4C_IKZF1_expression_distribution"), 95, 90)
save_pair(draw_d, file.path(out, "SFig4D_coordinate_alignment_QC"), 183, 90)
cat("SFig4B-D typography rerender complete\n")
