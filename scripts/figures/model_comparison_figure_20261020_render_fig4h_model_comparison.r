#!/usr/bin/env Rscript

# ==========================================================
# Script   : model_comparison_figure_20261020_render_fig4h_model_comparison.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:model_comparison_figure_20261020/render_fig4h_model_comparison.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "model_comparison_figure_20261020")
dir.create(file.path(out, "figures"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out, "source_data"), recursive = TRUE, showWarnings = FALSE)

p64 <- read.delim(file.path(root, "objective64_evidence_package_20261020/results/NESTED_CV_REPEAT_LEVEL_RESULTS.tsv"), check.names = FALSE)
ps4 <- read.delim(file.path(root, "sparse4_oof_estimand_repair_20261020/results/REPEAT_LEVEL_SCALE_VALID_C_AND_HR.tsv"), check.names = FALSE)
null <- read.delim(file.path(root, "objective64_evidence_package_20261020/results/PERMUTATION_NULL.tsv"), check.names = FALSE)
obs <- read.delim(file.path(root, "objective64_evidence_package_20261020/results/PERMUTATION_OBSERVED_STATISTIC.tsv"), check.names = FALSE)
uncertainty <- read.delim(file.path(root, "final_audit_reconciliation_20261022/production/uncertainty/SPARSE4_OOF_JACKKNIFE_INTERVAL.tsv"), check.names = FALSE)

stopifnot(nrow(p64) == 20L, nrow(ps4) == 20L, nrow(null) == 499L, nrow(obs) == 1L)
stopifnot(all(is.finite(p64$repeat_cindex)), all(is.finite(ps4$fold_stratified_harrell_c)))

repdat <- rbind(
  data.frame(procedure = "Prespecified 64-gene (lambda.min)", repeat_id = p64[["repeat"]],
             cindex = p64$repeat_cindex, stringsAsFactors = FALSE),
  data.frame(procedure = "Post hoc sparse4 (global 1-SE)", repeat_id = seq_len(nrow(ps4)),
             cindex = ps4$fold_stratified_harrell_c, stringsAsFactors = FALSE))
repdat$procedure <- factor(repdat$procedure,
  levels = c("Prespecified 64-gene (lambda.min)", "Post hoc sparse4 (global 1-SE)"))
sumdat <- aggregate(cindex ~ procedure, repdat,
                    function(x) c(mean = mean(x), sd = sd(x)))
sumdat <- data.frame(procedure = sumdat$procedure,
                     mean = sumdat$cindex[, "mean"], sd = sumdat$cindex[, "sd"])
sumdat$label <- c(
  sprintf("Mean %.3f\nSD %.3f", sumdat$mean[1], sumdat$sd[1]),
  sprintf("OOF %.3f\n95%% CI %.3f-%.3f\nSD %.3f",
          uncertainty$estimate[1], uncertainty$ci_low_normal_95[1], uncertainty$ci_high_normal_95[1], sumdat$sd[2])
)

cols <- c("Prespecified 64-gene (lambda.min)" = "#356B8C",
          "Post hoc sparse4 (global 1-SE)" = "#C47A4A")

set.seed(20261020)
p_left <- ggplot(repdat, aes(procedure, cindex, colour = procedure)) +
  geom_hline(yintercept = 0.5, linewidth = 0.35, linetype = 2, colour = "grey55") +
  geom_boxplot(width = 0.46, outlier.shape = NA, linewidth = 0.45,
               fill = "white", colour = "grey35") +
  geom_jitter(width = 0.10, height = 0, size = 1.35, alpha = 0.88) +
  geom_point(data = sumdat, aes(x = procedure, y = mean), inherit.aes = FALSE,
             shape = 23, size = 2.2, stroke = 0.45, fill = "white", colour = "black") +
  geom_text(data = sumdat,
            aes(x = procedure, y = 0.665,
                label = label),
            inherit.aes = FALSE, size = 2.60, lineheight = 0.95, colour = "black") +
  scale_colour_manual(values = cols) +
  scale_x_discrete(labels = c("64-gene\nlambda.min",
                              "sparse4\nglobal 1-SE")) +
  coord_cartesian(ylim = c(0.48, 0.69), clip = "off") +
  labs(x = NULL, y = "Fold-comparable Harrell C",
       title = "Conditional repeated CV") +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(legend.position = "none", axis.text = element_text(size = 7.5),
        axis.text.x = element_text(size = 7.5),
        axis.title.y = element_text(size = 7.5), plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(4, 4, 3, 3))

p_right <- ggplot(null, aes(mean_outer_test_cindex)) +
  geom_histogram(bins = 28, fill = "grey78", colour = "white", linewidth = 0.2) +
  geom_vline(xintercept = obs$mean_outer_test_cindex[1], colour = "#356B8C", linewidth = 0.8) +
  annotate("text", x = obs$mean_outer_test_cindex[1] - 0.004, y = Inf,
           label = sprintf("Observed %.3f\nP=0.002", obs$mean_outer_test_cindex[1]),
           hjust = 1, vjust = 1.35, size = 2.60, lineheight = 0.95, colour = "#244C66") +
  coord_cartesian(xlim = c(min(null$mean_outer_test_cindex) - 0.005,
                           obs$mean_outer_test_cindex[1] + 0.01), clip = "off") +
  labs(x = "Mean outer-test Harrell C", y = "Permutations",
       title = "64-gene null (n=499)") +
  theme_classic(base_size = 7.5, base_family = "Arial") +
  theme(axis.title = element_text(size = 7.5), axis.text = element_text(size = 7.5),
        plot.title = element_text(size = 8, face = "bold"),
        plot.margin = margin(4, 3, 3, 4))

fig <- p_left + p_right + plot_layout(widths = c(1.08, 0.92))

width_mm <- 105
height_mm <- 85
grDevices::cairo_pdf(file.path(out, "figures/Fig4H_conditional_model_comparison.pdf"),
                     width = width_mm / 25.4, height = height_mm / 25.4,
                     family = "Arial", onefile = FALSE)
print(fig)
dev.off()
ragg::agg_png(file.path(out, "figures/Fig4H_conditional_model_comparison.png"),
              width = width_mm, height = height_mm, units = "mm", res = 600,
              background = "white")
print(fig)
dev.off()

write.table(repdat, file.path(out, "source_data/Fig4H_repeat_C_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(null, file.path(out, "source_data/Fig4H_64gene_permutation_null.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(obs, file.path(out, "source_data/Fig4H_64gene_observed.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "SESSIONINFO.txt"))
