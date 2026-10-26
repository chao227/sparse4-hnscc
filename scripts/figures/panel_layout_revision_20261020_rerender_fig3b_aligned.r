#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_fig3b_aligned.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_fig3b_aligned.R
# ==========================================================

# Presentation-only successor for Fig.3B. It uses the frozen representative
# sections and Yellow/TLS10 values, but draws H&E and both spot maps from one
# shared official-coordinate crop per section so their tissue shapes align.

suppressPackageStartupMessages({
  library(png)
  library(jsonlite)
  library(ragg)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
input_root <- file.path(root, "inputs/spatial")
spot_path <- file.path(root, "yellow_module_spatial_20261019/production/source_data/YELLOW_SPATIAL_SPOT_LEVEL.tsv.gz")
manifest_path <- file.path(root, "replacement_spatial_figures_20261019/source_data/GSE208253_sample_manifest.tsv")
official_coord_path <- file.path(input_root, "official/GSE208253_official_spatial_coordinates.tsv")
out1 <- file.path(root, "yellow_module_spatial_20261019/render/figures/Fig3/Fig3B_yellow_representative_spatial_maps")
out2 <- file.path(root, "fig3_split_20261020/figures/Fig3/Fig3B_yellow_representative_spatial_maps")
tmp <- file.path(root, "panel_layout_revision_20261020/qa/fig3b_tiles")
dir.create(tmp, recursive = TRUE, showWarnings = FALSE)

spot_legacy <- read.delim(gzfile(spot_path), check.names = FALSE, stringsAsFactors = FALSE)
official <- read.delim(official_coord_path, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(!anyDuplicated(official[c("GSM", "barcode")]))
coord_key <- paste(official$GSM, official$barcode, sep = "\r")
spot_key <- paste(spot_legacy$GSM, spot_legacy$barcode, sep = "\r")
coord_index <- match(spot_key, coord_key)
stopifnot(!anyNA(coord_index))
spot <- spot_legacy
spot$pixel_col_fullres <- official$pixel_col_fullres[coord_index]
spot$pixel_row_fullres <- official$pixel_row_fullres[coord_index]
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
representative <- c("GSM6339638", "GSM6339633", "GSM6339632")
stopifnot(nrow(spot) == 23882L, all(representative %in% spot$GSM), all(representative %in% manifest$GSM))

read_gz_png <- function(path) {
  con <- gzfile(path, "rb"); on.exit(close(con), add = TRUE)
  readPNG(readBin(con, "raw", n = 1e9))
}
read_gz_json <- function(path) {
  con <- gzfile(path, "rt"); on.exit(close(con), add = TRUE)
  fromJSON(paste(readLines(con, warn = FALSE), collapse = "\n"))
}
make_crop <- function(img, data, scale, pad = 0.035) {
  x <- data$pixel_col_fullres * scale$tissue_hires_scalef
  y <- data$pixel_row_fullres * scale$tissue_hires_scalef
  xr <- range(x); yr <- range(y)
  px <- max(10, diff(xr) * pad); py <- max(10, diff(yr) * pad)
  list(x0 = max(0, xr[1] - px), x1 = min(dim(img)[2], xr[2] + px),
       y0 = max(0, yr[1] - py), y1 = min(dim(img)[1], yr[2] + py),
       w = dim(img)[2], h = dim(img)[1])
}

yellow_limits <- as.numeric(quantile(spot$yellow_activity_kME, c(.01, .99), na.rm = TRUE))
tls_limits <- as.numeric(quantile(spot$TLS10_score, c(.01, .99), na.rm = TRUE))
pal_yellow <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256)
pal_tls <- colorRampPalette(c("#F1F8F2", "#A8DDB5", "#43A86B", "#005A32"))(256)

sections <- lapply(representative, function(gsm) {
  row <- manifest[match(gsm, manifest$GSM), , drop = FALSE]
  section <- spot[spot$GSM == gsm, , drop = FALSE]
  stem <- paste0(gsm, "_", tolower(row$section_id))
  image_path <- file.path(input_root, "histology", paste0(stem, "_tissue_hires_image.png.gz"))
  scale_path <- file.path(input_root, "histology", paste0(stem, "_scalefactors_json.json.gz"))
  stopifnot(file.exists(image_path), file.exists(scale_path))
  image <- read_gz_png(image_path)
  scale <- read_gz_json(scale_path)
  list(gsm = gsm, patient = row$patient_id, section_id = row$section_id,
       data = section, image = image, scale = scale,
       crop = make_crop(image, section, scale))
})

render_tile <- function(z, kind, values = NULL, limits = NULL, palette = NULL) {
  path <- file.path(tmp, paste0(z$gsm, "_", kind, ".png"))
  ragg::agg_png(path, width = 1800, height = 1400, units = "px", res = 600, background = "white")
  par(mar = c(0, 0, 0, 0), xpd = NA)
  cr <- z$crop
  plot.new()
  plot.window(xlim = c(cr$x0, cr$x1), ylim = c(cr$y1, cr$y0), asp = 1, xaxs = "i", yaxs = "i")
  if (kind == "he") {
    rasterImage(z$image, 0, cr$h, cr$w, 0, interpolate = FALSE)
  } else {
    # Keep the identical H&E field behind each feature map. This makes the
    # barcode-coordinate registration directly visible rather than asking the
    # reader to compare an isolated spot silhouette with a separate H&E crop.
    bg <- z$image
    bg[, , 1:3] <- pmin(1, pmax(0, bg[, , 1:3] * 0.28 + 0.68))
    rasterImage(bg, 0, cr$h, cr$w, 0, interpolate = FALSE)
    x <- z$data$pixel_col_fullres * z$scale$tissue_hires_scalef
    y <- z$data$pixel_row_fullres * z$scale$tissue_hires_scalef
    idx <- pmax(1, pmin(length(palette), round((values - limits[1]) / diff(limits) * (length(palette) - 1)) + 1))
    diameter <- z$scale$spot_diameter_fullres * z$scale$tissue_hires_scalef
    symbols(x, y, circles = rep(diameter / 2, length(x)), inches = FALSE,
            bg = grDevices::adjustcolor(palette[idx], alpha.f = 0.86), fg = NA, add = TRUE)
  }
  box(col = "grey45", lwd = 0.45)
  dev.off()
  path
}

tiles <- lapply(sections, function(z) {
  c(he = render_tile(z, "he"),
    yellow = render_tile(z, "yellow", z$data$yellow_activity_kME, yellow_limits, pal_yellow),
    tls = render_tile(z, "tls", z$data$TLS10_score, tls_limits, pal_tls))
})

draw_tile <- function(path, x0, y0, x1, y1) rasterImage(readPNG(path), x0, y0, x1, y1, interpolate = FALSE)
draw_key <- function(x0, y0, width, height, limits, palette, label) {
  n <- length(palette)
  for (i in seq_len(n)) {
    rect(x0 + (i - 1) / n * width, y0, x0 + i / n * width, y0 + height,
         col = palette[i], border = NA)
  }
  text(x0, y0 - .008, formatC(limits[1], format = "f", digits = 2), adj = c(0, 1), cex = .62)
  text(x0 + width, y0 - .008, formatC(limits[2], format = "f", digits = 2), adj = c(1, 1), cex = .62)
  text(x0 + width / 2, y0 + height + .010, label, cex = .62, font = 2)
}

draw_panel <- function() {
  par(family = "Arial", mar = c(0, 0, 0, 0), xpd = NA)
  plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(.5, .982, "Representative multi-section spatial maps", font = 2, cex = .82)
  xs <- c(.045, .355, .665); width <- .285
  text(xs + width / 2, .925, c("H&E", "Yellow activity", "TLS10 program"), font = 2, cex = .62)
  yy <- c(.845, .605, .365); height <- .185
  for (i in seq_along(sections)) {
    z <- sections[[i]]
    for (j in seq_len(3)) draw_tile(tiles[[i]][j], xs[j], yy[i] - height, xs[j] + width, yy[i])
    text(.010, yy[i] - height / 2, paste0(z$patient, "-", z$section_id),
         font = 2, cex = .62, adj = c(0, .5))
  }
  # Compact keys occupy approximately one tenth of the previous legend area.
  draw_key(.665, .045, .115, .010, yellow_limits, pal_yellow, "Yellow activity")
  draw_key(.825, .045, .115, .010, tls_limits, pal_tls, "TLS10")
}

save_pair <- function(stem) {
  ragg::agg_png(paste0(stem, ".png"), width = 140 / 25.4, height = 110 / 25.4,
                units = "in", res = 600, background = "white")
  draw_panel(); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = 140 / 25.4, height = 110 / 25.4,
                       family = "Arial", onefile = TRUE)
  draw_panel(); dev.off()
}

save_pair(out1)
file.copy(paste0(out1, ".pdf"), paste0(out2, ".pdf"), overwrite = TRUE)
file.copy(paste0(out1, ".png"), paste0(out2, ".png"), overwrite = TRUE)

write.table(
  data.frame(
    GSM = representative,
    patient = vapply(sections, `[[`, character(1), "patient"),
    section_id = vapply(sections, `[[`, character(1), "section_id"),
    n_spots = vapply(sections, function(z) nrow(z$data), integer(1)),
    shared_crop = TRUE,
    rotation = "NONE", mirror = "NONE", manual_shift = "NONE"
  ),
  file.path(root, "panel_layout_revision_20261020/reports/FIG3B_ALIGNMENT_QA.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(
    check = c("all_QC_spots_barcode_matched", "legacy_col_equals_official_row", "legacy_row_equals_official_col"),
    value = c(
      !anyNA(coord_index),
      isTRUE(all.equal(spot_legacy$pixel_col_fullres, spot$pixel_row_fullres, tolerance = 0)),
      isTRUE(all.equal(spot_legacy$pixel_row_fullres, spot$pixel_col_fullres, tolerance = 0))
    ),
    interpretation = c(
      "Official coordinates recovered for every plotted QC-retained spot",
      "Confirms the legacy display-table column label swap",
      "Confirms the legacy display-table column label swap"
    )
  ),
  file.path(root, "panel_layout_revision_20261020/reports/FIG3B_COORDINATE_RECONCILIATION.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("Fig3B aligned rerender complete\n")
