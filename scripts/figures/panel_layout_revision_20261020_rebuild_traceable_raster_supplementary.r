#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rebuild_traceable_raster_supplementary.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rebuild_traceable_raster_supplementary.R
# ==========================================================

# Rebuild three previously raster-only Supplementary panels from traceable
# reback inputs. No statistical definition, patient roster, cell assignment or
# marker set is changed; only vector rendering and typography are repaired.

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(survminer)
  library(ragg)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
out <- file.path(root, "panel_layout_revision_20261020/figures/rebuilt_raster_supplementary")
source_out <- file.path(root, "panel_layout_revision_20261020/source_data")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(source_out, recursive = TRUE, showWarnings = FALSE)

save_pair <- function(plot, stem, width_mm, height_mm) {
  ragg::agg_png(paste0(stem, ".png"), width = width_mm / 25.4, height = height_mm / 25.4,
                units = "in", res = 600, background = "white")
  print(plot); dev.off()
  grDevices::cairo_pdf(paste0(stem, ".pdf"), width = width_mm / 25.4, height = height_mm / 25.4,
                       family = "Arial", onefile = TRUE)
  print(plot); dev.off()
}

theme_pub <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      axis.text = element_text(size = 7.5, colour = "grey10"),
      axis.title = element_text(size = 8),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7.5),
      panel.grid = element_blank()
    )
}

# -------------------------------------------------------------------------
# SFig2E: original TLS21 definition reconstructed from the frozen legacy VST
# and 476-patient survival table. This exactly reproduces the legacy 238/238
# split and risk-table counts at 0/50/100/150/200 months.
# -------------------------------------------------------------------------
master <- read.csv(file.path(root, "inputs/legacy_results/master_survival_IES_table.csv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
vst <- readRDS(file.path(root, "inputs/wgcna/legacy/vst_top25mad.rds"))
tls12 <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
bcell_tls <- c("MS4A1", "CD79B", "IGHM", "IGHG1", "IGHG3", "IGKC", "JCHAIN", "MZB1", "CD38", "CXCL13")
tls21 <- unique(c(tls12, bcell_tls))
stopifnot(length(tls21) == 21L, all(tls21 %in% rownames(vst)), nrow(master) == 476L)
sample_index <- match(master$sample_id, substr(colnames(vst), 1, 15))
stopifnot(!anyNA(sample_index))
z_by_gene <- apply(vst[tls21, sample_index, drop = FALSE], 1, function(value) (value - mean(value)) / sd(value))
master$TLS21_score <- rowMeans(z_by_gene)
master$TLS21_group <- factor(
  ifelse(master$TLS21_score >= median(master$TLS21_score), "TLS-High", "TLS-Low"),
  levels = c("TLS-High", "TLS-Low")
)
stopifnot(identical(as.integer(table(master$TLS21_group)), c(238L, 238L)))

fit <- survfit(Surv(OS_time_months, OS_status) ~ TLS21_group, data = master)
curve <- survminer::surv_summary(fit, data = master)
curve$strata <- sub("TLS21_group=", "", curve$strata)
curve$strata <- factor(curve$strata, levels = c("TLS-High", "TLS-Low"))
curve <- rbind(
  data.frame(time = 0, n.risk = 238, n.event = 0, n.censor = 0, surv = 1,
             std.err = 0, upper = 1, lower = 1, strata = factor(c("TLS-High"), levels = levels(curve$strata)),
             TLS21_group = "TLS-High"),
  data.frame(time = 0, n.risk = 238, n.event = 0, n.censor = 0, surv = 1,
             std.err = 0, upper = 1, lower = 1, strata = factor(c("TLS-Low"), levels = levels(curve$strata)),
             TLS21_group = "TLS-Low"),
  curve
)
curve$strata <- factor(as.character(curve$strata), levels = c("TLS-High", "TLS-Low"))

risk_times <- c(0, 50, 100, 150, 200)
risk_summary <- summary(fit, times = risk_times, extend = TRUE)
risk <- data.frame(
  group = factor(sub("TLS21_group=", "", risk_summary$strata), levels = c("TLS-High", "TLS-Low")),
  time = risk_summary$time,
  n_risk = risk_summary$n.risk
)
expected_risk <- c(238, 40, 8, 5, 1, 238, 35, 7, 2, 0)
stopifnot(identical(as.integer(risk$n_risk), as.integer(expected_risk)))

logrank <- survdiff(Surv(OS_time_months, OS_status) ~ TLS21_group, data = master)
logrank_p <- pchisq(logrank$chisq, df = 1, lower.tail = FALSE)
cox <- coxph(Surv(OS_time_months, OS_status) ~ relevel(TLS21_group, ref = "TLS-Low"), data = master)
cox_ci <- summary(cox)$conf.int[1, ]
write.table(
  data.frame(
    patient_id = master$patient_id, OS_time_months = master$OS_time_months,
    OS_status = master$OS_status, TLS21_score = master$TLS21_score,
    TLS21_group = master$TLS21_group
  ),
  file.path(source_out, "SFIG2E_TLS21_SURVIVAL_SOURCE.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  data.frame(n = nrow(master), events = sum(master$OS_status), logrank_p = logrank_p,
             HR_TLS_high_vs_low = cox_ci["exp(coef)"], CI_lower = cox_ci["lower .95"],
             CI_upper = cox_ci["upper .95"]),
  file.path(source_out, "SFIG2E_TLS21_SURVIVAL_SUMMARY.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

palette_surv <- c("TLS-High" = "#E64B35", "TLS-Low" = "#4DBBD5")
p_curve <- ggplot(curve, aes(time, surv, colour = strata)) +
  geom_step(linewidth = 0.7) +
  geom_point(data = subset(curve, n.censor > 0), shape = 3, size = 1.2, stroke = 0.45) +
  scale_colour_manual(values = palette_surv, labels = c("TLS-High", "TLS-Low")) +
  scale_x_continuous(limits = c(0, 215), breaks = risk_times, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0, 0)) +
  labs(title = "TLS21 score and overall survival", subtitle = sprintf("n = 476; log-rank P = %.3g", logrank_p),
       x = "Time (months)", y = "Overall survival probability", colour = NULL) +
  theme_pub(8) +
  theme(legend.position = c(0.78, 0.82), legend.background = element_rect(fill = "white", colour = NA),
        plot.margin = margin(3, 4, 1, 2))

risk$y <- ifelse(risk$group == "TLS-High", 2, 1)
risk$text_time <- ifelse(risk$time == 0, 1, risk$time)
risk$text_hjust <- ifelse(risk$time == 0, 0, 0.5)
p_risk <- ggplot(risk, aes(text_time, y, label = n_risk, hjust = text_hjust)) +
  geom_text(size = 2.8, family = "Arial", show.legend = FALSE) +
  scale_x_continuous(limits = c(0, 215), breaks = risk_times, expand = c(0, 0)) +
  scale_y_continuous(limits = c(0.5, 2.5), breaks = c(1, 2), labels = c("TLS-Low", "TLS-High")) +
  labs(title = "Number at risk", x = NULL, y = NULL) +
  coord_cartesian(xlim = c(0, 215), clip = "off") +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(axis.text.x = element_text(size = 7.5), axis.text.y = element_text(size = 7.5, face = "bold"),
        axis.ticks.y = element_blank(), axis.line.y = element_blank(),
        plot.title = element_text(size = 8, face = "bold", hjust = 0),
        plot.margin = margin(0, 5, 5, 5))

p_sfig2e <- p_curve / p_risk + plot_layout(heights = c(0.80, 0.20), guides = "collect")
save_pair(p_sfig2e, file.path(out, "SFig2E_TLS21_survival_vector"), 82, 92)

# -------------------------------------------------------------------------
# SFig3A/B: use the frozen CD8 subcluster object (20,522 cells), its saved UMAP
# embedding, saved subtype assignments and normalized RNA layer.
# -------------------------------------------------------------------------
cd8 <- readRDS(file.path(root, "inputs/singlecell/subcluster/cd8_subclusters.rds"))
state_levels <- c("Naive/Stem-like", "Progenitor-Tex", "Intermediate-Tex", "Terminal-Tex", "Cytotoxic-Eff", "Memory")
state_colours <- c(
  "Naive/Stem-like" = "#2C7FB8", "Progenitor-Tex" = "#7FCDBB",
  "Intermediate-Tex" = "#F2C84B", "Terminal-Tex" = "#D95F4E",
  "Cytotoxic-Eff" = "#E59F38", "Memory" = "#78C679"
)
stopifnot(ncol(cd8) == 20522L, setequal(as.character(unique(cd8$cd8_subtype)), state_levels))
embedding <- Embeddings(cd8, "umap")
umap <- data.frame(UMAP_1 = embedding[, 1], UMAP_2 = embedding[, 2],
                   subtype = factor(as.character(cd8$cd8_subtype), levels = state_levels))
write.table(umap, file.path(source_out, "SFIG3A_CD8_UMAP_SOURCE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_sfig3a <- ggplot(umap, aes(UMAP_1, UMAP_2, colour = subtype)) +
  geom_point(size = 0.20, alpha = 0.78, stroke = 0) +
  scale_colour_manual(values = state_colours, drop = FALSE) +
  coord_equal() +
  labs(title = "CD8+ T cell subtypes", subtitle = "Six states; n = 20,522 cells",
       x = "UMAP 1", y = "UMAP 2", colour = "Subtype") +
  theme_pub(8) +
  theme(legend.position = "right", legend.key.height = grid::unit(3.6, "mm"),
        plot.margin = margin(3, 4, 3, 3))
save_pair(p_sfig3a, file.path(out, "SFig3A_CD8_state_UMAP_vector"), 82, 78)

markers <- c("TCF7", "LEF1", "CCR7", "TOX", "PDCD1", "TIGIT", "HAVCR2", "LAYN", "ENTPD1", "GZMB", "PRF1", "CX3CR1", "GZMK", "IL7R")
stopifnot(all(markers %in% rownames(cd8)))
expression <- as.matrix(GetAssayData(cd8, assay = "RNA", layer = "data")[markers, , drop = FALSE])
marker_long <- do.call(rbind, lapply(markers, function(gene) {
  data.frame(gene = gene, subtype = factor(as.character(cd8$cd8_subtype), levels = state_levels),
             expression = as.numeric(expression[gene, ]), stringsAsFactors = FALSE)
}))
marker_long$gene <- factor(marker_long$gene, levels = rev(markers))
write.table(marker_long, file.path(source_out, "SFIG3B_CD8_MARKER_SOURCE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p_sfig3b <- ggplot(marker_long, aes(subtype, expression, fill = subtype)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25, alpha = 0.88) +
  facet_grid(gene ~ ., scales = "free_y", switch = "y") +
  scale_fill_manual(values = state_colours, guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "CD8+ T cell marker expression across subtypes", x = NULL, y = "Log-normalized expression") +
  theme_classic(base_family = "Arial", base_size = 8) +
  theme(
    axis.text.x = element_text(size = 7.5, angle = 35, hjust = 1, colour = "grey10"),
    axis.text.y = element_text(size = 7.5, colour = "grey10"),
    axis.title.y = element_text(size = 8),
    strip.text.y.left = element_text(size = 7.5, angle = 0, face = "italic", hjust = 1),
    strip.placement = "outside", strip.background = element_blank(),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
    panel.spacing.y = grid::unit(0.08, "lines"),
    plot.margin = margin(3, 5, 4, 4)
  )
save_pair(p_sfig3b, file.path(out, "SFig3B_CD8_marker_violin_vector"), 110, 125)

cat("Traceable vector rebuild complete for SFig2E, SFig3A and SFig3B\n")
