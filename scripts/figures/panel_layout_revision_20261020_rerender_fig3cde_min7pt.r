#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_fig3cde_min7pt.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_fig3cde_min7pt.R
# ==========================================================

# Typography-only successor renderer for Fig. 3C-E. The statistical source
# tables, estimands, patient unit and display geometry are inherited unchanged
# from the audited Yellow-module spatial renderer.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1]))
project_root <- normalizePath(file.path(dirname(script_path), "../.."))
source_root <- file.path(project_root, "yellow_module_spatial_20261019", "production", "source_data")
yellow_out <- file.path(project_root, "yellow_module_spatial_20261019", "render", "figures", "Fig3")
split_out <- file.path(project_root, "fig3_split_20261020", "figures", "Fig3")
dir.create(yellow_out, recursive = TRUE, showWarnings = FALSE)
dir.create(split_out, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(name) {
  path <- file.path(source_root, name)
  if (!file.exists(path)) stop("Missing source table: ", path, call. = FALSE)
  read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
}
patient <- read_tsv("YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv")
moran <- read_tsv("YELLOW_SPATIAL_PATIENT_MORAN.tsv")
primary <- read_tsv("YELLOW_SPATIAL_PRIMARY_INFERENCE.tsv")
sensitivity <- read_tsv("YELLOW_SPATIAL_SENSITIVITY_SUMMARY.tsv")
shared_diagnostic <- read_tsv("YELLOW_SPATIAL_SHARED_GENE_DIAGNOSTIC.tsv")
context <- read_tsv("YELLOW_SPATIAL_CONTEXT_SUMMARY.tsv")
uncertainty_path <- file.path(project_root, "final_audit_reconciliation_20261022", "production",
                              "uncertainty", "SPATIAL_YELLOW_TLS10_INTERVALS.tsv")
if (!file.exists(uncertainty_path)) stop("Missing independently audited uncertainty source: ", uncertainty_path)
uncertainty <- read.delim(uncertainty_path, check.names = FALSE, stringsAsFactors = FALSE)

fixed_patients <- paste0("P", 1:10)
stopifnot(
  nrow(patient) == 10L,
  setequal(patient$patient_id, fixed_patients),
  nrow(primary) == 1L,
  primary$n_independent_patients[1] == 10L,
  primary$n_sections[1] == 12L,
  nrow(shared_diagnostic) == 1L,
  shared_diagnostic$protocol_status[1] == "POST_RESULT_AMENDMENT__NOT_PRESPECIFIED"
)

format_fixed <- function(x, digits = 3) formatC(x, format = "f", digits = digits)
range_from_zero <- function(x, pad = 0.18) {
  value_range <- range(c(0, x), na.rm = TRUE)
  delta <- diff(value_range)
  if (!is.finite(delta) || delta == 0) delta <- 1
  c(value_range[1] - delta * 0.08, value_range[2] + delta * pad)
}

theme_min7 <- function() {
  theme_classic(base_size = 8, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.text = element_text(size = 7.0, colour = "grey10"),
      axis.title = element_text(size = 8.0),
      panel.grid = element_blank(),
      legend.text = element_text(size = 7.0),
      legend.title = element_text(size = 7.0),
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.0, face = "bold", hjust = 0.5),
      plot.caption = element_text(size = 7.0),
      plot.margin = margin(3, 4, 3, 4)
    )
}

patient_dot_plot <- function(data, column, title, subtitle, colour, x_label = "Spearman rho") {
  values <- as.numeric(data[[column]])
  stopifnot(length(values) == 10L, all(is.finite(values)))
  display <- data.frame(
    patient_id = factor(data$patient_id, levels = rev(fixed_patients)),
    value = values
  )
  ggplot(display, aes(value, patient_id)) +
    geom_vline(xintercept = 0, linetype = 2, colour = "grey55", linewidth = 0.3) +
    geom_point(size = 2, colour = colour) +
    coord_cartesian(xlim = range_from_zero(values), ylim = c(0.5, 10.5)) +
    labs(title = title, subtitle = subtitle, x = x_label, y = "Patient ID") +
    theme_min7() +
    theme(plot.subtitle = element_text(size = 7.0, face = "bold", hjust = 0.5, colour = colour))
}

save_pair <- function(plot, stem, width_mm, height_mm) {
  for (out_dir in c(yellow_out, split_out)) {
    prefix <- file.path(out_dir, stem)
    ragg::agg_png(paste0(prefix, ".png"), width = width_mm / 25.4,
                  height = height_mm / 25.4, units = "in", res = 600,
                  background = "white")
    print(plot)
    dev.off()
    grDevices::cairo_pdf(paste0(prefix, ".pdf"), width = width_mm / 25.4,
                         height = height_mm / 25.4, family = "Arial", onefile = TRUE)
    print(plot)
    dev.off()
  }
}

# Fig. 3C: formal patient-level primary test.
panel_c_subtitle <- paste0(
  "Mean rho = ", format_fixed(uncertainty$estimate[1]),
  " (95% CI ", format_fixed(uncertainty$t_ci_low_95[1]),
  "-", format_fixed(uncertainty$t_ci_high_95[1]), ")\n",
  primary$positive_patients[1], "/10 positive; median = ",
  format_fixed(primary$median_patient_rho[1]), "; exact P = .00195"
)
panel_c <- patient_dot_plot(
  patient, "yellow_TLS10_rho", "Yellow activity-TLS10 association",
  panel_c_subtitle, "#7A3E9D"
)
save_pair(panel_c, "Fig3C_yellow_TLS10_patient_association", 65, 60)

# Fig. 3D: unchanged prespecified sensitivities and post-result overlap check.
sensitivity_long <- rbind(
  data.frame(patient_id = patient$patient_id, method = "Primary kME", rho = patient$yellow_TLS10_rho),
  data.frame(patient_id = patient$patient_id, method = "Sign-only", rho = patient$sign_only_TLS10_rho),
  data.frame(patient_id = patient$patient_id, method = "Sparse4-excluded", rho = patient$sparse4_excluded_TLS10_rho),
  data.frame(patient_id = patient$patient_id, method = "CD79A-overlap diag.", rho = patient$yellow_minus_CD79A_TLS9_minus_CD79A_rho)
)
sensitivity_long$method <- factor(
  sensitivity_long$method,
  levels = c("Primary kME", "Sign-only", "Sparse4-excluded", "CD79A-overlap diag.")
)
sensitivity_medians <- data.frame(
  method = factor(levels(sensitivity_long$method), levels = levels(sensitivity_long$method)),
  rho = c(
    primary$median_patient_rho[1],
    sensitivity$median_patient_rho[sensitivity$analysis == "sign-only all mapped Yellow genes"],
    sensitivity$median_patient_rho[sensitivity$analysis == "primary kME activity excluding frozen sparse4 genes"],
    shared_diagnostic$median_patient_rho[1]
  )
)
stopifnot(nrow(sensitivity_long) == 40L, nrow(sensitivity_medians) == 4L,
          all(is.finite(sensitivity_long$rho)), all(is.finite(sensitivity_medians$rho)))
panel_d <- ggplot(sensitivity_long, aes(method, rho, group = patient_id)) +
  geom_hline(yintercept = 0, linetype = 2, colour = "grey60", linewidth = 0.25) +
  geom_line(colour = "grey75", linewidth = 0.25, alpha = 0.8) +
  geom_point(colour = "#7A3E9D", size = 0.9, alpha = 0.75) +
  geom_crossbar(
    data = sensitivity_medians,
    aes(x = method, y = rho, ymin = rho, ymax = rho, group = method),
    width = 0.55, linewidth = 0.4, colour = "black", inherit.aes = FALSE
  ) +
  labs(
    title = "Activity robustness",
    subtitle = "CD79A overlap check is post-result\nAll methods: 10/10 positive",
    x = NULL, y = "TLS10 Spearman rho"
  ) +
  scale_x_discrete(labels = c(
    "Primary kME" = "Primary\nkME",
    "Sign-only" = "Sign-\nonly",
    "Sparse4-excluded" = "Exclude\nsparse4",
    "CD79A-overlap diag." = "Exclude\nCD79A\n(post-result)"
  )) +
  theme_min7() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 7.0),
    plot.subtitle = element_text(size = 7.0, face = "bold", hjust = 0.5,
                                 colour = "#7A3E9D", lineheight = 1.05)
  )
save_pair(panel_d, "Fig3D_yellow_activity_sensitivity", 65, 60)

# Fig. 3E: unchanged descriptive patient-level neighbourhood statistic.
neighbourhood_context <- context[
  context$analysis == "Yellow focal activity versus neighboring TLS10 rank-based bivariate Moran I", ,
  drop = FALSE
]
stopifnot(nrow(neighbourhood_context) == 1L)
panel_e <- patient_dot_plot(
  moran, "yellow_focal_TLS10_neighbor_bivariate_Moran_I",
  "TLS10-pattern neighbourhood association",
  paste0(neighbourhood_context$positive_patients, "/10 positive; median I = ",
         format_fixed(neighbourhood_context$median_patient_effect)),
  "#26734D", "Bivariate Moran's I"
)
save_pair(panel_e, "Fig3E_yellow_TLS10_neighbourhood", 80, 70)

cat("Rendered Fig3C-E with all visible text >= 7 pt; data and geometry unchanged.\n")
