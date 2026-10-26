#!/usr/bin/env Rscript

# ==========================================================
# Script   : targeted_panel_adjudication_20261023_render_targeted_panels.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:targeted_panel_adjudication_20261023/scripts/render_targeted_panels.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(patchwork)
  library(grid)
})

root <- normalizePath(getwd(), mustWork = TRUE)
task <- file.path(root, "targeted_panel_adjudication_20261023")
out <- file.path(task, "figures")
srcout <- file.path(task, "source_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(srcout, recursive = TRUE, showWarnings = FALSE)

pt_mm <- function(x) x / 2.845276
save_pair <- function(plot, id, width_mm, height_mm) {
  grDevices::cairo_pdf(file.path(out, paste0(id, ".pdf")),
                       width = width_mm / 25.4, height = height_mm / 25.4,
                       family = "Arial", onefile = FALSE)
  print(plot)
  dev.off()
  ragg::agg_png(file.path(out, paste0(id, ".png")),
                width = width_mm, height = height_mm, units = "mm", res = 600,
                background = "white")
  print(plot)
  dev.off()
}

theme_pub <- function(base_size = 7.5) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5,
                                margin = margin(b = 2)),
      axis.title = element_text(size = 7.5),
      axis.text = element_text(size = 7),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 7),
      plot.margin = margin(3, 3, 3, 3)
    )
}

# Fig.5C: the analysis is unchanged; this rerender makes the already reported
# log-rank result visible in the standalone panel.
d <- read.delim(file.path(root, "legacy_fig5_fig6_recompute_20261020/source_data/PATIENT_LEVEL_MASTER_476.tsv"),
                check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(d) == 476L, sum(d$OS_status) == 220L,
          all(c("OS_time_months", "OS_status", "IKZF1", "CD79A", "quadrant") %in% names(d)))
d$ikzf1_group <- factor(ifelse(d$IKZF1 >= median(d$IKZF1), "High", "Low"),
                         levels = c("High", "Low"))
stopifnot(identical(as.integer(table(d$ikzf1_group)), c(238L, 238L)))
fit <- survfit(Surv(OS_time_months, OS_status) ~ ikzf1_group, data = d)
lr <- survdiff(Surv(OS_time_months, OS_status) ~ ikzf1_group, data = d)
p_lr <- unname(pchisq(lr$chisq, df = 1, lower.tail = FALSE))
stopifnot(abs(p_lr - 0.00466) < 0.00002)
ss <- summary(fit)
km <- data.frame(time = ss$time, survival = ss$surv,
                 group = sub("^ikzf1_group=", "", ss$strata),
                 n_censor = ss$n.censor, stringsAsFactors = FALSE)
km$group <- factor(km$group, levels = c("High", "Low"))
cens <- km[km$n_censor > 0, , drop = FALSE]
times <- c(0, 50, 100, 150, 200)
rs <- summary(fit, times = times, extend = TRUE)
risk <- data.frame(time = rs$time, n_risk = rs$n.risk,
                   group = sub("^ikzf1_group=", "", rs$strata),
                   stringsAsFactors = FALSE)
risk$group <- factor(risk$group, levels = c("Low", "High"))
cols <- c(High = "#D73027", Low = "#4575B4")
p_main <- ggplot(km, aes(time, survival, colour = group)) +
  geom_step(linewidth = 0.55) +
  geom_point(data = cens, shape = 3, size = 1.15, stroke = 0.45) +
  geom_hline(yintercept = 0.5, linetype = 2, linewidth = 0.35) +
  scale_colour_manual(values = cols, name = NULL) +
  scale_x_continuous(breaks = times, expand = c(0.01, 0)) +
  scale_y_continuous(breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  coord_cartesian(xlim = c(-18, 210), ylim = c(0, 1.02)) +
  labs(title = "IKZF1 Expression and\nOverall Survival",
       subtitle = sprintf("Log-rank P = %.5f; n = 476", p_lr),
       x = NULL, y = "Overall survival") +
  theme_pub() +
  theme(plot.subtitle = element_text(size = 7, hjust = 0.5, colour = "grey20",
                                     margin = margin(b = 1.5)),
        legend.position = c(0.82, 0.80),
        legend.background = element_blank(), legend.key = element_blank())
p_risk <- ggplot(risk, aes(time, group, label = n_risk, colour = group)) +
  geom_text(size = pt_mm(7), family = "Arial", show.legend = FALSE) +
  scale_colour_manual(values = cols) +
  scale_x_continuous(limits = c(-18, 210), breaks = times, expand = c(0.01, 0)) +
  labs(title = "Number at risk", x = "Time (months)", y = NULL) +
  theme_pub() +
  theme(plot.title = element_text(size = 7.5, face = "bold", hjust = 0.5),
        axis.text.y = element_text(colour = "grey20", margin = margin(r = 8)),
        panel.grid = element_blank())
p5c <- p_main / p_risk + plot_layout(heights = c(3.25, 1.15))
save_pair(p5c, "Fig5C", 60, 60)
write.table(data.frame(n = nrow(d), events = sum(d$OS_status),
                       n_high = sum(d$ikzf1_group == "High"),
                       n_low = sum(d$ikzf1_group == "Low"),
                       logrank_p = p_lr),
            file.path(srcout, "Fig5C_statistics.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)
write.table(risk, file.path(srcout, "Fig5C_risk_table.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

# Fig.5D: retain the four-group definition and medians but remove the duplicate
# Spearman annotation already reported in Fig.5B.
qd <- d[, c("patient_id", "sample_id", "CD79A", "IKZF1", "quadrant")]
qd$quadrant <- factor(qd$quadrant,
                      levels = c("LymphHot", "B-skewed", "T-skewed", "Cold"))
stopifnot(identical(as.integer(table(qd$quadrant)), c(171L, 67L, 67L, 171L)))
med_cd <- median(qd$CD79A)
med_ik <- median(qd$IKZF1)
qcols <- c("LymphHot" = "#D73027", "B-skewed" = "#FC8D59",
           "T-skewed" = "#91BFDB", "Cold" = "#4575B4")
p5d <- ggplot(qd, aes(CD79A, IKZF1, colour = quadrant)) +
  geom_hline(yintercept = med_ik, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = med_cd, linetype = 2, colour = "grey50", linewidth = 0.4) +
  geom_point(size = 0.75, alpha = 0.5) +
  scale_colour_manual(values = qcols, name = NULL,
                      labels = c("LymphHot\n(CD79A+ / IKZF1+)",
                                 "B-skewed\n(CD79A+ / IKZF1-)",
                                 "T-skewed\n(CD79A- / IKZF1+)",
                                 "Cold\n(CD79A- / IKZF1-)")) +
  guides(colour = guide_legend(ncol = 2, override.aes = list(size = 2, alpha = 1))) +
  labs(title = "CD79A-by-IKZF1\nQuadrant Classification (TCGA-HNSC)",
       x = "CD79A expression (VST)", y = "IKZF1 expression (VST)") +
  theme_pub() +
  theme(legend.position = "bottom", legend.key.size = unit(3, "mm"))
save_pair(p5d, "Fig5D", 80, 80)
write.table(data.frame(cd79a_median = med_cd, ikzf1_median = med_ik,
                       LymphHot = 171L, B_skewed = 67L,
                       T_skewed = 67L, Cold = 171L),
            file.path(srcout, "Fig5D_cutoffs_and_counts.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(task, "SESSIONINFO.txt"))
