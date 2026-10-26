#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_fig7b_from_frozen_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_fig7b_from_frozen_source.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
  library(svglite)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
source_file <- file.path(root, "results/F_singlecell/IKZF1_CD8_FIG7B_SOURCE_DATA.tsv")
figure_dir <- file.path(root, "figures/F_singlecell")
plot_data <- read.delim(source_file, check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(plot_data) == 6L, all(plot_data$paired_donors == 18L),
          all(is.finite(plot_data$median_paired_difference_high_minus_low)),
          all(is.finite(plot_data$p_adjust_BH)))
plot_data$gene <- factor(plot_data$gene, levels = rev(plot_data$gene))
plot_data$q_label <- sprintf("q = %.3g", plot_data$p_adjust_BH)
plot_data$BH_status <- factor(plot_data$BH_status,
                              levels = c("BH q < 0.05", "BH q >= 0.05"))

fig7b <- ggplot(plot_data,
                aes(x = median_paired_difference_high_minus_low,
                    y = gene, colour = BH_status)) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "#606060") +
  geom_errorbar(aes(xmin = paired_difference_bootstrap_CI_lower,
                    xmax = paired_difference_bootstrap_CI_upper),
                orientation = "y", width = 0, linewidth = 0.65) +
  geom_point(size = 2.2) +
  geom_text(aes(x = 0.075, label = q_label), hjust = 0, colour = "#202020",
            family = "Arial", size = 2.50) +
  scale_colour_manual(values = c("BH q < 0.05" = "#2F6F9F",
                                 "BH q >= 0.05" = "#9AA0A6")) +
  coord_cartesian(xlim = c(-0.44, 0.19), clip = "off") +
  labs(x = "Median paired difference (log2 CPM)\nIKZF1-detected/high - IKZF1-undetected/zero",
       y = NULL, colour = NULL,
       subtitle = "n = 18 paired donors; 95% paired-donor bootstrap CI") +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(axis.line = element_line(linewidth = 0.35, colour = "black"),
        axis.ticks = element_line(linewidth = 0.35, colour = "black"),
        axis.text = element_text(size = 7, colour = "black"),
        axis.text.y = element_text(face = "italic"),
        axis.title.x = element_text(size = 7.2, margin = margin(t = 5)),
        plot.subtitle = element_text(size = 7.2, margin = margin(b = 6)),
        legend.position = "bottom", legend.text = element_text(size = 7.0),
        legend.key.width = grid::unit(7, "pt"),
        plot.margin = margin(6, 18, 6, 6)) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

figure_base <- file.path(figure_dir, "Fig7B_IKZF1_CD8_DONOR_PSEUDOBULK")
w <- 89 / 25.4
h <- 82 / 25.4
svglite::svglite(paste0(figure_base, ".svg"), width = w, height = h)
print(fig7b); dev.off()
grDevices::cairo_pdf(paste0(figure_base, ".pdf"), width = w, height = h,
                     family = "Arial", onefile = TRUE)
print(fig7b); dev.off()
ragg::agg_png(paste0(figure_base, ".png"), width = w, height = h,
              units = "in", res = 600, background = "white")
print(fig7b); dev.off()
