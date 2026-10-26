#!/usr/bin/env Rscript

# ==========================================================
# Script   : lymphhot_external_figure_20261020_render_fig6e_lymphhot_external.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:lymphhot_external_figure_20261020/render_fig6e_lymphhot_external.R
# ==========================================================

suppressPackageStartupMessages(library(ggplot2))
root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "lymphhot_external_figure_20261020")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "source_data"), recursive = TRUE, showWarnings = FALSE)

src <- read.delim(file.path(root, "agent_work/lymphhot_formal_close_20261020/results/LYMPHHOT_EXTERNAL_RESULTS.tsv"),
                  check.names = FALSE, stringsAsFactors = FALSE)
d <- src[src$status == "COMPUTED_PENDING_INDEPENDENT_AUDIT", , drop = FALSE]
stopifnot(nrow(d) == 2L, all(is.finite(d$hazard_ratio)), all(d$n > 0), all(d$events > 0))
d$status <- "POST_HOC_EXTERNAL_EVALUATION"
d$evaluation_status <- "POST_HOC_EXTERNAL_EVALUATION"
d$shared_family_status <- "DESCRIPTIVE_TWO_EVALUABLE_COHORT_BH_FAMILY"
d$cohort <- factor(d$cohort, levels = rev(c("GSE41613", "GSE65858")))
d$label <- sprintf("HR %.3f [%.3f-%.3f]; P=%.3f\nn=%d; %d events",
                   d$hazard_ratio, d$ci_lower, d$ci_upper, d$p_nominal,
                   d$n, d$events)

p <- ggplot(d, aes(hazard_ratio, cohort)) +
  geom_vline(xintercept = 1, linetype = 2, linewidth = 0.35, colour = "grey55") +
  geom_errorbar(aes(xmin = ci_lower, xmax = ci_upper), orientation = "y",
                width = 0.12, linewidth = 0.55, colour = "#315F7D") +
  geom_point(shape = 21, size = 2.4, stroke = 0.55, fill = "white", colour = "#315F7D") +
  geom_text(aes(x = 2.55, label = label), hjust = 1, nudge_y = -0.27,
            size = 2.50, colour = "grey20") +
  scale_x_continuous(limits = c(0.35, 2.60), breaks = c(0.5, 1.0, 1.5, 2.0, 2.5)) +
  labs(x = "LymphHot hazard ratio (95% CI)", y = NULL,
       title = "Post hoc external evaluation") +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(plot.title = element_text(size = 7.8, face = "bold", hjust = 0),
        axis.title.x = element_text(size = 7.0),
        axis.text = element_text(size = 7.0),
        plot.margin = margin(4, 4, 4, 4))

width_mm <- 60
height_mm <- 60
grDevices::cairo_pdf(file.path(out, "figures/Fig6E_LymphHot_external_evaluation.pdf"),
                     width = width_mm / 25.4, height = height_mm / 25.4,
                     family = "Arial", onefile = FALSE)
print(p)
dev.off()
ragg::agg_png(file.path(out, "figures/Fig6E_LymphHot_external_evaluation.png"),
              width = width_mm, height = height_mm, units = "mm", res = 600,
              background = "white")
print(p)
dev.off()

write.table(d, file.path(out, "source_data/Fig6E_LymphHot_external_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(c(
  "GSE159067 was not plotted because exact IKZF1 coverage failed.",
  "No proxy, fixture or simulated result was substituted.",
  "External cohorts were inspected before formal model/phenotype adoption; results are post hoc evaluations."),
  file.path(out, "source_data/Fig6E_NOT_FEASIBLE_NOTE.txt"))
capture.output(sessionInfo(), file = file.path(out, "SESSIONINFO.txt"))
