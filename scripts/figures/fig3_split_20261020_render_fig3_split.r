#!/usr/bin/env Rscript

# ==========================================================
# Script   : fig3_split_20261020_render_fig3_split.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:fig3_split_20261020/scripts/render_fig3_split.R
# ==========================================================

# Successor renderer that splits the previously combined spatial-context panel
# into independent Fig. 3F and Fig. 3G panels. No statistical model is refit.

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(png)
  library(ragg)
  library(grid)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1]))
module_root <- normalizePath(file.path(dirname(script_path), ".."))
project_root <- normalizePath(file.path(module_root, ".."))
source_module <- file.path(project_root, "yellow_module_spatial_20261019")
source_fig <- file.path(source_module, "render", "figures", "Fig3")
source_data <- file.path(source_module, "production", "source_data")
out_fig <- file.path(module_root, "figures", "Fig3")
out_src <- file.path(module_root, "source_data")
dir.create(out_fig, recursive = TRUE, showWarnings = FALSE)
dir.create(out_src, recursive = TRUE, showWarnings = FALSE)

patient <- read.delim(file.path(source_data, "YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv"), check.names = FALSE)
context <- read.delim(file.path(source_data, "YELLOW_SPATIAL_CONTEXT_SUMMARY.tsv"), check.names = FALSE)
spot <- read.delim(gzfile(file.path(source_data, "YELLOW_SPATIAL_SPOT_LEVEL.tsv.gz")), check.names = FALSE)
stopifnot(nrow(patient) == 10L, length(unique(patient$patient_id)) == 10L, nrow(spot) == 23882L)

fixed_patients <- paste0("P", 1:10)
format_fixed <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
range_from_zero <- function(x, pad = 0.18) {
  value_range <- range(c(0, x), na.rm = TRUE)
  delta <- diff(value_range)
  if (!is.finite(delta) || delta == 0) delta <- 1
  c(value_range[1] - delta * 0.08, value_range[2] + delta * pad)
}
theme_nature <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(size = 7.2, family = "Arial"),
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(size = 7.2, colour = "grey10"),
      axis.title = element_text(size = 8),
      legend.text = element_text(size = 7.2),
      legend.title = element_text(size = 7.2),
      strip.text = element_text(size = 7.2),
      panel.grid = element_blank(),
      plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5, lineheight = 0.95),
      plot.subtitle = element_text(size = 7.2, face = "bold", hjust = 0.5, lineheight = 0.95),
      plot.caption = element_text(size = 7.2),
      plot.margin = margin(3, 4, 3, 4)
    )
}
patient_dot_plot <- function(column, title, subtitle, colour) {
  values <- as.numeric(patient[[column]])
  display <- data.frame(patient_id = factor(patient$patient_id, levels = rev(fixed_patients)), value = values)
  ggplot(display, aes(value, patient_id)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_point(size = 2, colour = colour) +
    coord_cartesian(xlim = range_from_zero(values), ylim = c(0.5, 10.5)) +
    labs(title = title, subtitle = subtitle, x = "Spearman rho", y = "Patient ID") +
    theme_nature() +
    theme(plot.subtitle = element_text(colour = colour))
}
save_pair <- function(plot, stem, width_mm, height_mm) {
  ragg::agg_png(paste0(stem, ".png"), width = width_mm / 25.4, height = height_mm / 25.4,
                units = "in", res = 600, background = "white")
  print(plot)
  dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4,
                       height = height_mm / 25.4, family = "Arial", onefile = TRUE)
  print(plot)
  dev.off()
}

rctd <- context[context$analysis == "Yellow activity versus RCTD B-cell weight", ]
ikzf1 <- context[context$analysis == "Yellow activity versus IKZF1", ]
stopifnot(nrow(rctd) == 1L, nrow(ikzf1) == 1L)
ikzf1_detected <- sum(spot$IKZF1_log1p_norm > 0, na.rm = TRUE)
ikzf1_detection_rate <- ikzf1_detected / nrow(spot)

panel_f <- patient_dot_plot(
  "yellow_RCTD_B_cell_rho",
  "Yellow activity-RCTD\nB-cell weight",
  paste0(rctd$positive_patients, "/10 positive; median rho = ", format_fixed(rctd$median_patient_effect)),
  "#3E6D8E"
)
panel_g <- patient_dot_plot(
  "yellow_IKZF1_rho",
  "Yellow activity-direct\nIKZF1 expression",
  paste0(ikzf1$positive_patients, "/10 positive; median rho = ", format_fixed(ikzf1$median_patient_effect),
         "\nIKZF1 detected in ", formatC(100 * ikzf1_detection_rate, format = "f", digits = 1), "% of spots"),
  "#B45F06"
)
save_pair(panel_f, file.path(out_fig, "Fig3F_yellow_RCTD_Bcell_association"), 60, 70)
save_pair(panel_g, file.path(out_fig, "Fig3G_yellow_IKZF1_association"), 60, 70)

# Copy already audited A-E pairs byte-for-byte into the successor package.
stems <- c(
  "Fig3A_spatial_cohort_overview",
  "Fig3B_yellow_representative_spatial_maps",
  "Fig3C_yellow_TLS10_patient_association",
  "Fig3D_yellow_activity_sensitivity",
  "Fig3E_yellow_TLS10_neighbourhood"
)
for (stem in stems) for (ext in c("pdf", "png")) {
  src <- file.path(source_fig, paste0(stem, ".", ext))
  stopifnot(file.exists(src), file.copy(src, file.path(out_fig, basename(src)), overwrite = TRUE))
}

png_files <- file.path(out_fig, paste0(c(stems, "Fig3F_yellow_RCTD_Bcell_association", "Fig3G_yellow_IKZF1_association"), ".png"))
stopifnot(all(file.exists(png_files)))
grobs <- lapply(png_files, function(path) wrap_elements(full = rasterGrob(readPNG(path), interpolate = FALSE)))
composite <- (grobs[[1]] | grobs[[2]]) /
  (grobs[[3]] | grobs[[4]]) /
  (grobs[[5]] | grobs[[6]] | grobs[[7]]) +
  plot_layout(heights = c(1.04, 0.92, 0.92)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(family = "Arial", face = "bold", size = 8))
save_pair(composite, file.path(out_fig, "Figure3_yellow_7panel_FINAL_submission"), 183, 190)

write.table(
  data.frame(
    panel = c("Fig3F", "Fig3G"),
    estimand = c("patient-mean within-section Spearman rho: Yellow activity vs RCTD B-cell weight",
                 "patient-mean within-section Spearman rho: Yellow activity vs direct unsmoothed IKZF1"),
    n_patients = c(rctd$n_independent_patients, ikzf1$n_independent_patients),
    positive_patients = c(rctd$positive_patients, ikzf1$positive_patients),
    median_rho = c(rctd$median_patient_effect, ikzf1$median_patient_effect),
    additional_context = c("RCTD weight is computational composition context, not histology",
                           paste0("IKZF1 detected in ", ikzf1_detected, "/", nrow(spot), " spots (", formatC(100 * ikzf1_detection_rate, format = "f", digits = 1), "%)")),
    claim_ceiling = c("descriptive patient-level spatial composition context only",
                      "descriptive patient-level spatial expression context only")
  ),
  file.path(out_src, "Fig3FG_source_data.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)

write.table(
  data.frame(panel = c(paste0("Fig3", LETTERS[1:7]), "Figure3_composite"),
             width_mm = c(200, 120, 65, 65, 80, 60, 60, 183),
             height_mm = c(100, 100, 60, 60, 70, 70, 70, 190)),
  file.path(out_src, "FIGURE_GEOMETRY_CONTRACT.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
