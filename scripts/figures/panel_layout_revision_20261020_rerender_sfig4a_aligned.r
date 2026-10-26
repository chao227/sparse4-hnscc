#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_sfig4a_aligned.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_sfig4a_aligned.R
# ==========================================================

# Presentation-only successor for SFig.4A. Frozen Yellow/TLS10 values are
# retained. Display coordinates are recovered by GSM + barcode from the
# authoritative official 10x coordinate table, and H&E plus both spot maps use
# the same crop, orientation and aspect ratio within each section.

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
output_stem <- file.path(root, "yellow_module_spatial_20261019/render/figures/SFig4/SFig4A_yellow_complete_spatial_landscape")
qa_root <- file.path(root, "panel_layout_revision_20261020")
tile_root <- file.path(qa_root, "qa/sfig4a_tiles")
dir.create(tile_root, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(qa_root, "reports"), recursive = TRUE, showWarnings = FALSE)

required <- c(spot_path, manifest_path, official_coord_path)
stopifnot(all(file.exists(required)))

spot_legacy <- read.delim(gzfile(spot_path), check.names = FALSE, stringsAsFactors = FALSE)
official <- read.delim(official_coord_path, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)

required_spot <- c("GSM", "barcode", "patient_id", "section_id", "yellow_activity_kME", "TLS10_score")
required_official <- c("GSM", "barcode", "pixel_col_fullres", "pixel_row_fullres")
stopifnot(all(required_spot %in% names(spot_legacy)), all(required_official %in% names(official)))
stopifnot(!anyDuplicated(official[c("GSM", "barcode")]))

# Never use the two legacy pixel columns for rendering. They are retained only
# below for a reconciliation audit that documents the historical label swap.
official_key <- paste(official$GSM, official$barcode, sep = "\r")
spot_key <- paste(spot_legacy$GSM, spot_legacy$barcode, sep = "\r")
official_index <- match(spot_key, official_key)
stopifnot(!anyNA(official_index))
spot <- spot_legacy
spot$pixel_col_fullres <- official$pixel_col_fullres[official_index]
spot$pixel_row_fullres <- official$pixel_row_fullres[official_index]

manifest <- manifest[manifest$eligible == "YES", , drop = FALSE]
manifest <- manifest[order(as.integer(sub("^GSM63396", "", manifest$GSM))), , drop = FALSE]
stopifnot(nrow(manifest) == 12L, setequal(manifest$GSM, unique(spot$GSM)))

read_gz_png <- function(path) {
  stopifnot(file.exists(path))
  connection <- gzfile(path, "rb")
  on.exit(close(connection), add = TRUE)
  readPNG(readBin(connection, "raw", n = 1e9))
}

read_gz_json <- function(path) {
  stopifnot(file.exists(path))
  connection <- gzfile(path, "rt")
  on.exit(close(connection), add = TRUE)
  fromJSON(paste(readLines(connection, warn = FALSE), collapse = "\n"))
}

make_crop <- function(image, section, scalefactor, pad = 0.035) {
  x <- section$pixel_col_fullres * scalefactor$tissue_hires_scalef
  y <- section$pixel_row_fullres * scalefactor$tissue_hires_scalef
  x_range <- range(x)
  y_range <- range(y)
  x_pad <- max(10, diff(x_range) * pad)
  y_pad <- max(10, diff(y_range) * pad)
  list(
    x0 = max(0, x_range[1] - x_pad),
    x1 = min(dim(image)[2], x_range[2] + x_pad),
    y0 = max(0, y_range[1] - y_pad),
    y1 = min(dim(image)[1], y_range[2] + y_pad),
    width = dim(image)[2],
    height = dim(image)[1]
  )
}

yellow_limits <- as.numeric(quantile(spot$yellow_activity_kME, c(0.01, 0.99), na.rm = TRUE))
tls_limits <- as.numeric(quantile(spot$TLS10_score, c(0.01, 0.99), na.rm = TRUE))
yellow_palette <- colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(256)
tls_palette <- colorRampPalette(c("#F1F8F2", "#A8DDB5", "#43A86B", "#005A32"))(256)

sections <- lapply(seq_len(nrow(manifest)), function(index) {
  manifest_row <- manifest[index, , drop = FALSE]
  section <- spot[spot$GSM == manifest_row$GSM, , drop = FALSE]
  stem <- paste0(manifest_row$GSM, "_", tolower(manifest_row$section_id))
  image_path <- file.path(input_root, "histology", paste0(stem, "_tissue_hires_image.png.gz"))
  scale_path <- file.path(input_root, "histology", paste0(stem, "_scalefactors_json.json.gz"))
  image <- read_gz_png(image_path)
  scale <- read_gz_json(scale_path)
  list(
    gsm = manifest_row$GSM,
    patient = manifest_row$patient_id,
    section_id = manifest_row$section_id,
    data = section,
    image = image,
    scale = scale,
    crop = make_crop(image, section, scale)
  )
})

render_tile <- function(section, kind, values = NULL, limits = NULL, palette = NULL) {
  path <- file.path(tile_root, paste0(section$gsm, "_", kind, ".png"))
  # Match the physical aspect ratio of the final 14.8 x 20.5 mm tile so the
  # cached raster is not stretched during panel assembly.
  ragg::agg_png(path, width = 720, height = 1000, units = "px", res = 600, background = "white")
  par(mar = c(0, 0, 0, 0), xpd = NA)
  crop <- section$crop
  plot.new()
  plot.window(xlim = c(crop$x0, crop$x1), ylim = c(crop$y1, crop$y0), asp = 1, xaxs = "i", yaxs = "i")
  if (kind == "he") {
    rasterImage(section$image, 0, crop$height, crop$width, 0, interpolate = FALSE)
  } else {
    background <- section$image
    background[, , 1:3] <- pmin(1, pmax(0, background[, , 1:3] * 0.28 + 0.68))
    rasterImage(background, 0, crop$height, crop$width, 0, interpolate = FALSE)
    scaled <- (values - limits[1]) / diff(limits)
    palette_index <- pmax(1, pmin(length(palette), round(scaled * (length(palette) - 1)) + 1))
    x <- section$data$pixel_col_fullres * section$scale$tissue_hires_scalef
    y <- section$data$pixel_row_fullres * section$scale$tissue_hires_scalef
    diameter <- section$scale$spot_diameter_fullres * section$scale$tissue_hires_scalef
    symbols(
      x, y, circles = rep(diameter / 2, length(x)), inches = FALSE,
      bg = adjustcolor(palette[palette_index], alpha.f = 0.86), fg = NA, add = TRUE
    )
  }
  box(col = "grey45", lwd = 0.35)
  dev.off()
  path
}

tile_paths <- lapply(sections, function(section) {
  c(
    he = render_tile(section, "he"),
    yellow = render_tile(section, "yellow", section$data$yellow_activity_kME, yellow_limits, yellow_palette),
    tls = render_tile(section, "tls", section$data$TLS10_score, tls_limits, tls_palette)
  )
})

draw_tile <- function(path, x0, y0, x1, y1) {
  rasterImage(readPNG(path), x0, y0, x1, y1, interpolate = FALSE)
}

draw_key <- function(x0, y0, width, height, limits, palette, label) {
  n <- length(palette)
  for (index in seq_len(n)) {
    rect(
      x0 + (index - 1) / n * width, y0,
      x0 + index / n * width, y0 + height,
      col = palette[index], border = NA
    )
  }
  text(x0, y0 - 0.004, formatC(limits[1], format = "f", digits = 2), adj = c(0, 1), cex = 0.65)
  text(x0 + width, y0 - 0.004, formatC(limits[2], format = "f", digits = 2), adj = c(1, 1), cex = 0.65)
  text(x0 + width / 2, y0 + height + 0.004, label, cex = 0.65, font = 2)
}

draw_panel <- function() {
  par(family = "Arial", mar = c(0, 0, 0, 0), xpd = NA)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  text(0.5, 0.982, "Yellow-module spatial activity and TLS10 transcriptional program across 12 sections", font = 2, cex = 0.72)

  column_left <- seq(0.012, 0.988, length.out = 5)
  row_top <- c(0.915, 0.625, 0.335)
  tile_height <- 0.205
  for (index in seq_along(sections)) {
    column <- (index - 1) %% 4 + 1
    row <- ceiling(index / 4)
    block_x0 <- column_left[column]
    block_x1 <- column_left[column + 1] - 0.008
    tile_width <- (block_x1 - block_x0 - 0.008) / 3
    tile_y1 <- row_top[row]
    tile_y0 <- tile_y1 - tile_height
    for (map_index in seq_len(3)) {
      x0 <- block_x0 + (map_index - 1) * (tile_width + 0.004)
      draw_tile(tile_paths[[index]][map_index], x0, tile_y0, x0 + tile_width, tile_y1)
      text(
        x0 + tile_width / 2, tile_y1 + 0.010,
        if (map_index == 1) paste0(sections[[index]]$patient, "-", sections[[index]]$section_id) else c("Yellow", "TLS10")[map_index - 1],
        font = 2, cex = 0.65
      )
    }
  }

  # Compact keys: approximately one tenth of the previous key-region height,
  # with the scale definitions retained in the Supplementary Figure legend.
  draw_key(0.655, 0.020, 0.120, 0.006, yellow_limits, yellow_palette, "Yellow activity")
  draw_key(0.825, 0.020, 0.120, 0.006, tls_limits, tls_palette, "TLS10")
}

save_pair <- function(stem) {
  ragg::agg_png(
    paste0(stem, ".png"), width = 200 / 25.4, height = 100 / 25.4,
    units = "in", res = 600, background = "white"
  )
  draw_panel()
  dev.off()
  grDevices::cairo_pdf(
    paste0(stem, ".pdf"), width = 200 / 25.4, height = 100 / 25.4,
    family = "Arial", onefile = TRUE
  )
  draw_panel()
  dev.off()
}

save_pair(output_stem)

alignment_qa <- do.call(rbind, lapply(sections, function(section) {
  official_subset <- official[official$GSM == section$gsm, , drop = FALSE]
  data.frame(
    GSM = section$gsm,
    patient = section$patient,
    section_id = section$section_id,
    n_plotted_spots = nrow(section$data),
    n_official_in_tissue_spots = nrow(official_subset),
    all_QC_spots_barcode_matched = all(section$data$barcode %in% official_subset$barcode),
    shared_crop_HE_Yellow_TLS10 = TRUE,
    rotation = "NONE",
    mirror = "NONE",
    manual_shift = "NONE"
  )
}))
write.table(
  alignment_qa,
  file.path(qa_root, "reports/SFIG4A_ALIGNMENT_QA.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

coordinate_qa <- data.frame(
  check = c(
    "all_QC_spots_barcode_matched",
    "legacy_pixel_col_equals_official_pixel_row",
    "legacy_pixel_row_equals_official_pixel_col",
    "frozen_yellow_values_preserved",
    "frozen_TLS10_values_preserved"
  ),
  value = c(
    !anyNA(official_index),
    isTRUE(all.equal(spot_legacy$pixel_col_fullres, spot$pixel_row_fullres, tolerance = 0)),
    isTRUE(all.equal(spot_legacy$pixel_row_fullres, spot$pixel_col_fullres, tolerance = 0)),
    isTRUE(all.equal(spot_legacy$yellow_activity_kME, spot$yellow_activity_kME, tolerance = 0)),
    isTRUE(all.equal(spot_legacy$TLS10_score, spot$TLS10_score, tolerance = 0))
  ),
  interpretation = c(
    "Every plotted QC-retained spot was joined to the official 10x coordinate table by GSM plus barcode",
    "Documents the legacy display-table coordinate-label swap",
    "Documents the legacy display-table coordinate-label swap",
    "Only display coordinates were replaced; Yellow values are byte-for-byte unchanged in memory",
    "Only display coordinates were replaced; TLS10 values are byte-for-byte unchanged in memory"
  )
)
write.table(
  coordinate_qa,
  file.path(qa_root, "reports/SFIG4A_COORDINATE_RECONCILIATION.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

cat("SFig4A authoritative-coordinate aligned rerender complete\n")
