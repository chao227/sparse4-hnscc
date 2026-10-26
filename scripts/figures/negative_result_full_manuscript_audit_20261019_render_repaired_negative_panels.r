#!/usr/bin/env Rscript

# ==========================================================
# Script   : negative_result_full_manuscript_audit_20261019_render_repaired_negative_panels.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:negative_result_full_manuscript_audit_20261019/scripts/render_repaired_negative_panels.R
# ==========================================================
suppressPackageStartupMessages(library(ggplot2))

args <- commandArgs(FALSE)
script_file <- sub("^--file=", "", args[grep("^--file=", args)])
repo <- normalizePath(file.path(dirname(script_file), "..", ".."))
module <- file.path(repo, "negative_result_full_manuscript_audit_20261019")
figdir <- file.path(module, "figures")
srcdir <- file.path(module, "source_data")
logdir <- file.path(module, "logs")
target_stems_env <- Sys.getenv("RENDER_TARGET_STEMS", unset = "")
target_stems <- if (nzchar(target_stems_env)) strsplit(target_stems_env, ",", fixed = TRUE)[[1]] else character()
target_sources_env <- Sys.getenv("RENDER_TARGET_SOURCES", unset = "")
target_sources <- if (nzchar(target_sources_env)) strsplit(target_sources_env, ",", fixed = TRUE)[[1]] else character()
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)
dir.create(logdir, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) {
  if (length(target_sources) && !name %in% target_sources && name != "REPAIRED_PANEL_GEOMETRY.tsv") return(invisible(NULL))
  write.table(x, file.path(srcdir, name), sep = "\t", row.names = FALSE, quote = FALSE, na = "NA")
}
save_pair <- function(p, stem, width_mm, height_mm) {
  if (length(target_stems) && !stem %in% target_stems) return(invisible(NULL))
  cairo_pdf(file.path(figdir, paste0(stem, ".pdf")), width = width_mm / 25.4,
            height = height_mm / 25.4, family = "Arial")
  print(p); dev.off()
  png(file.path(figdir, paste0(stem, ".png")), width = width_mm / 25.4,
      height = height_mm / 25.4, units = "in", res = 600, type = "cairo")
  print(p); dev.off()
}

base_theme <- theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.3),
        axis.ticks = element_line(linewidth = 0.3),
        plot.title = element_text(size = 8, face = "bold"),
        plot.subtitle = element_text(size = 7.2, lineheight = 0.95),
        plot.margin = margin(3, 4, 3, 3))

# Fig.1G: donor-level descriptive CellChat result, not a pooled-cell inference.
cc_path <- file.path(repo, "agent_work/cellchat_donor_descriptive_nobsb_20261016",
                     "formal_run_002/results/CELLCHAT_ENDPOINT_DESCRIPTIVE_SUMMARY.tsv")
cc <- read_tsv(cc_path)
stopifnot(nrow(cc) == 144L, length(unique(cc$direction_id)) == 24L,
          length(unique(cc$pathway)) == 6L, all(cc$n_total_donors == 18L),
          all(cc$claim_class == "DONOR_LEVEL_DESCRIPTIVE_CONDITIONAL_PROPENSITY"))
cc$direction <- paste(gsub("_", " ", cc$sender), "\u2192", gsub("_", " ", cc$receiver))
direction_order <- unique(cc[order(cc$endpoint_order), c("direction_id", "direction")])
cc$direction <- factor(cc$direction, levels = rev(direction_order$direction))
cc$pathway <- factor(cc$pathway, levels = unique(cc[order(cc$pathway_id), "pathway"]))
cc$cell_label <- sprintf("%d/%d", cc$n_detected, cc$n_eligible)
write_tsv(cc[, c("endpoint_id", "direction", "pathway", "n_total_donors", "n_eligible",
                 "n_ineligible", "n_detected", "n_nondetected_zero", "detection_fraction",
                 "median_probability", "claim_class")], "Fig1G_donor_cellchat_display.tsv")

p_cc <- ggplot(cc, aes(pathway, direction, fill = detection_fraction)) +
  geom_tile(color = "white", linewidth = 0.18) +
  geom_text(aes(label = cell_label), size = 2.50, family = "Arial", color = "black") +
  scale_fill_gradientn(colours = c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"),
                       limits = c(0, 1), breaks = c(0, .5, 1), name = "Detection fraction") +
  guides(fill = guide_colorbar(title.position = "top", title.hjust = .5,
                               direction = "vertical",
                               barwidth = grid::unit(2.0, "mm"),
                               barheight = grid::unit(25, "mm"))) +
  labs(title = "Donor-level CellChat summary",
       subtitle = "18 donors; cells show detected/eligible donors\nDescriptive conditional propensities",
       x = NULL, y = NULL) +
  base_theme +
  theme(plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 7.0, lineheight = 0.95),
        axis.text.x = element_text(angle = 35, hjust = 1, size = 7.0),
        axis.text.y = element_text(size = 7.0),
        legend.position = "right", legend.direction = "vertical",
        legend.title = element_text(size = 7.0, lineheight = 0.9),
        legend.text = element_text(size = 7.0),
        legend.margin = margin(0, 0, 0, 1),
        legend.box.spacing = grid::unit(0.5, "mm"))
save_pair(p_cc, "Fig1G_donor_level_CellChat", 105, 125)

# Fig.5A: registered patient-state display plus method/coverage sensitivity.
state_path <- file.path(repo, "scan_1019_remediation_20261019/inputs/IKZF1_patient_aware_AUCell.tsv")
grid_path <- file.path(module, "results/IKZF1_STATE_SENSITIVITY_GRID.tsv")
state <- read_tsv(state_path)
sens <- read_tsv(grid_path)
lev <- c("Naive/Stem-like", "Progenitor-Tex", "Intermediate-Tex", "Terminal-Tex", "Cytotoxic-Eff", "Memory")
stopifnot(nrow(state) == 108L, length(unique(state$patient)) == 18L,
          all(lev %in% state$CD8_state))
state$CD8_state <- factor(state$CD8_state, levels = lev)
needed <- c("median_all_states", "median_ncell_ge10", "median_ncell_ge20",
            "mean_ncell_ge10", "mean_ncell_ge20")
stopifnot(all(needed %in% sens$analysis))
sens_display <- sens[match(needed, sens$analysis), ]
write_tsv(state, "Fig5A_patient_state_values.tsv")
write_tsv(sens_display, "Fig5A_coverage_sensitivity.tsv")
getp <- function(id) sens$p_value[match(id, sens$analysis)]
subtitle_state <- sprintf("All-state P=%.3f; coverage-dependent\nSensitivity analyses P=%.3f\u2013%.3f",
                          getp("median_all_states"),
                          min(sens_display$p_value[-1], na.rm = TRUE),
                          max(sens_display$p_value[-1], na.rm = TRUE))
p_state <- ggplot(state, aes(CD8_state, IKZF1_AUCell_median, group = patient)) +
  geom_line(color = "grey73", linewidth = .23, alpha = .70) +
  geom_point(color = "#39769A", size = .62, alpha = .78) +
  stat_summary(aes(group = 1), fun = median, geom = "line", color = "#B33A3A", linewidth = .75) +
  stat_summary(aes(group = 1), fun = median, geom = "point", color = "#B33A3A", size = 1.45) +
  labs(title = "IKZF1 CD8-state sensitivity",
       subtitle = subtitle_state, x = NULL, y = "Patient-state median IKZF1 AUCell") +
  base_theme + theme(axis.text = element_text(size = 7.5, colour = "black"),
                     axis.title = element_text(size = 7.5),
                     plot.subtitle = element_text(size = 7.5, lineheight = .95),
                     axis.text.x = element_text(angle = 38, hjust = 1, size = 7.5))
save_pair(p_state, "Fig5A_IKZF1_state_coverage_sensitivity", 72, 66)

# SFig.3E: traceable gene-level unsigned virtual-perturbation result.
vko <- read_tsv(file.path(repo, "results/A_consistency/VKO_OUTPUT_AUDIT.tsv"))
vko <- vko[vko$gene != "IKZF1", ]
vko <- vko[order(vko$rank_distance), ][seq_len(min(12L, nrow(vko))), ]
vko$gene <- factor(vko$gene, levels = rev(vko$gene))
vko$supported <- ifelse(vko$p.adj < .05, "BH q < 0.05", "BH q \u2265 0.05")
vko$distance_scaled <- vko$distance * 1e5
stopifnot(all(vko$direction_available == FALSE), sum(vko$p.adj < .05) == 2L)
write_tsv(vko, "SFig3E_unsigned_VKO_top12.tsv")
p_vko <- ggplot(vko, aes(distance_scaled, gene, color = supported)) +
  geom_segment(aes(x = 0, xend = distance_scaled, yend = gene), color = "grey72", linewidth = .45) +
  geom_point(size = 1.9) +
  scale_color_manual(values = c("BH q < 0.05"="#B2182B", "BH q \u2265 0.05"="#6A8EAE")) +
  labs(title = "Unsigned IKZF1 perturbation",
       subtitle = "Top non-target genes; no expression direction",
       x = "Unsigned distance (x 10^-5)", y = NULL, color = NULL) +
  guides(color = guide_legend(direction = "vertical", ncol = 1,
                              byrow = TRUE,
                              keywidth = grid::unit(2.4, "mm"),
                              keyheight = grid::unit(2.4, "mm"))) +
  base_theme +
  theme(plot.title = element_text(size = 8.5, face = "bold"),
        plot.subtitle = element_text(size = 7.0, lineheight = 0.95),
        axis.title.x = element_text(size = 7.5),
        axis.text.x = element_text(size = 7.0),
        axis.text.y = element_text(size = 7.0),
        legend.position = "right",
        legend.direction = "vertical",
        legend.text = element_text(size = 7.0, lineheight = 0.9),
        legend.key.size = grid::unit(3.0, "mm"),
        legend.spacing.y = grid::unit(0.15, "mm"),
        legend.margin = margin(0, 0, 0, 1),
        legend.box.spacing = grid::unit(0.4, "mm"))
save_pair(p_vko, "SFig3E_unsigned_IKZF1_virtual_perturbation", 90, 65)

geometry <- data.frame(panel = c("Fig1G", "Fig5A", "SFig3E"),
                       width_mm = c(105, 72, 90), height_mm = c(125, 66, 65),
                       basis = rep("measured submitted composite bbox", 3))
write_tsv(geometry, "REPAIRED_PANEL_GEOMETRY.tsv")
writeLines(capture.output(sessionInfo()), file.path(logdir, "render_repaired_negative_panels_sessionInfo.txt"))
cat("Rendered Fig1G, Fig5A and SFig3E repaired panels\n")
