#!/usr/bin/env Rscript

# ==========================================================
# Script   : legacy_fig5_fig6_recompute_20261020_recompute_and_render_legacy_fig5_fig6.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:legacy_fig5_fig6_recompute_20261020/scripts/recompute_and_render_legacy_fig5_fig6.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(ragg)
})

ROOT <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
MOD <- file.path(ROOT, "legacy_fig5_fig6_recompute_20261020")
SRC <- file.path(MOD, "source_data")
FIG <- file.path(MOD, "figures")
REP <- file.path(MOD, "reports")
LOG <- file.path(MOD, "logs")
for (d in c(SRC, FIG, REP, LOG)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

stopif <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(SRC, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
z_rows <- function(m) {
  z <- t(scale(t(m)))
  stopif(!any(!is.finite(z)), "Non-finite row z-score")
  z
}
fmt_p <- function(p) ifelse(p < 0.001, "P < 0.001", sprintf("P = %.3f", p))

pal <- c("Cold" = "#4575B4", "B-skewed" = "#FC8D59", "T-skewed" = "#74ADD1", "LymphHot" = "#D73027")
lev <- names(pal)

theme_pub <- function() {
  theme_classic(base_size = 7.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(size = 7, colour = "black"),
      axis.title = element_text(size = 7.5),
      legend.title = element_text(size = 7),
      legend.text = element_text(size = 7),
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 7, face = "bold"),
      plot.margin = margin(2.5, 3, 2.5, 3, "mm")
    )
}

save_pair <- function(p, stem, width_mm, height_mm) {
  pdf_file <- file.path(FIG, paste0(stem, ".pdf"))
  png_file <- file.path(FIG, paste0(stem, ".png"))
  grDevices::cairo_pdf(pdf_file, width = width_mm / 25.4, height = height_mm / 25.4,
                       family = "Arial", onefile = FALSE)
  print(p)
  dev.off()
  ragg::agg_png(png_file, width = width_mm / 25.4, height = height_mm / 25.4,
                units = "in", res = 600, background = "white")
  print(p)
  dev.off()
}

# Guarded inputs and exact patient mapping.
surv <- read.csv(file.path(ROOT, "inputs/tcga/clinical/survival_data.csv"), check.names = FALSE, stringsAsFactors = FALSE)
vst <- readRDS(file.path(ROOT, "inputs/tcga/expression/vst_matrix.rds"))
tme <- read.csv(file.path(ROOT, "inputs/tcga/features/tme_scores.csv"), check.names = FALSE, stringsAsFactors = FALSE)
score <- read.delim(file.path(ROOT, "results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/SCORE_TABLE.tsv"), check.names = FALSE, stringsAsFactors = FALSE)

stopif(nrow(surv) == 476L && sum(surv$OS_status) == 220L, "Survival roster drift")
stopif(!anyNA(surv$sample_id) && !anyDuplicated(surv$sample_id), "Survival sample IDs invalid")
stopif(identical(dim(vst), c(29451L, 520L)), "VST dimensions drift")
stopif(!anyDuplicated(colnames(vst)), "VST sample IDs duplicated")
vi <- match(surv$sample_id, colnames(vst))
stopif(!anyNA(vi) && !anyDuplicated(vi), "Survival-to-VST mapping failed")
ti <- match(substr(surv$sample_id, 1, 15), substr(tme$sample_id, 1, 15))
stopif(!anyNA(ti) && !anyDuplicated(ti), "Survival-to-immune-score mapping failed")
zi <- match(surv$sample_id, score$sample_id)
stopif(!anyNA(zi) && !anyDuplicated(zi), "Survival-to-sparse4 mapping failed")
stopif(all(is.finite(score$linear_predictor[zi])), "Sparse4 score contains missing/non-finite values")

required_genes <- c("CD79A", "IKZF1", "PDCD1", "CD274", "TIGIT", "HAVCR2", "LAG3", "CTLA4",
                    "GZMB", "PRF1", "GNLY", "IFNG", "NKG7",
                    "CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21",
                    "CXCL9", "CXCL10", "CXCL11", "CXCL13")
stopif(all(required_genes %in% rownames(vst)), "Required gene coverage failed")
v <- vst[, vi, drop = FALSE]

cd <- as.numeric(v["CD79A", ])
ik <- as.numeric(v["IKZF1", ])
cd_cut <- median(cd)
ik_cut <- median(ik)
grp <- ifelse(cd >= cd_cut & ik >= ik_cut, "LymphHot",
              ifelse(cd >= cd_cut, "B-skewed", ifelse(ik >= ik_cut, "T-skewed", "Cold")))
grp <- factor(grp, levels = lev)
stopif(identical(as.integer(table(grp)), c(171L, 67L, 67L, 171L)), "Quadrant counts drift")

checkpoint_genes <- c("PDCD1", "CD274", "TIGIT", "HAVCR2", "LAG3", "CTLA4")
cytotoxic_genes <- c("GZMB", "PRF1", "GNLY", "IFNG", "NKG7")
tls12_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
checkpoint_z <- z_rows(v[checkpoint_genes, , drop = FALSE])
cytotoxic_z <- z_rows(v[cytotoxic_genes, , drop = FALSE])
tls12_z <- z_rows(v[tls12_genes, , drop = FALSE])

d <- data.frame(
  patient_id = substr(surv$sample_id, 1, 12), sample_id = surv$sample_id,
  OS_time_months = as.numeric(surv$OS_time), OS_status = as.integer(surv$OS_status),
  CD79A = cd, IKZF1 = ik, quadrant = grp,
  sparse4 = as.numeric(score$linear_predictor[zi]),
  checkpoint_score = colMeans(checkpoint_z), cytotoxic_score = colMeans(cytotoxic_z),
  project_immune_ssGSEA = as.numeric(tme$ImmuneScore[ti]), TLS12 = colMeans(tls12_z),
  stringsAsFactors = FALSE
)
stopif(all(complete.cases(d)), "Master table contains missing values")
write_tsv(d, "PATIENT_LEVEL_MASTER_476.tsv")
write_tsv(data.frame(marker = c("CD79A", "IKZF1"), cutoff = c(cd_cut, ik_cut),
                     rule = "high if value >= full-cohort median", n = nrow(d)), "QUADRANT_CUTOFFS.tsv")
write_tsv(data.frame(signature = c("Checkpoint", "Cytotoxic", "TLS12"),
                     genes = c(paste(checkpoint_genes, collapse = ";"), paste(cytotoxic_genes, collapse = ";"), paste(tls12_genes, collapse = ";")),
                     calculation = "gene-wise z score across the 476 patients, followed by equal-weight mean",
                     stringsAsFactors = FALSE), "SIGNATURE_DICTIONARY.tsv")

# Fig.5B: patient-level bulk association.
ct <- suppressWarnings(cor.test(d$CD79A, d$IKZF1, method = "spearman", exact = FALSE))
fig5b_stats <- data.frame(n = nrow(d), spearman_rho = unname(ct$estimate), p_value = ct$p.value,
                          unit = "patient", alternative_explanation = "shared immune-cell abundance")
write_tsv(fig5b_stats, "Fig5B_statistics.tsv")
write_tsv(d[, c("patient_id", "sample_id", "CD79A", "IKZF1", "quadrant")], "Fig5B_source_data.tsv")
p5b <- ggplot(d, aes(CD79A, IKZF1)) +
  geom_point(size = 0.8, alpha = 0.5, colour = "#3B6EA8") +
  geom_smooth(method = "lm", formula = y ~ x, se = TRUE, linewidth = 0.55, colour = "#B34E4E", fill = "#E7B4B4") +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.35,
           label = sprintf("Spearman rho = %.3f\nn = 476, %s", unname(ct$estimate), fmt_p(ct$p.value)),
           size = 2.5, family = "Arial") +
  labs(title = "Bulk CD79A and IKZF1 expression", x = "CD79A (VST)", y = "IKZF1 (VST)") + theme_pub()
save_pair(p5b, "Fig5B_CD79A_IKZF1_bulk_association", 60, 60)

# Fig.5E: four-group descriptive Kaplan-Meier curve.
sf <- survfit(Surv(OS_time_months, OS_status) ~ quadrant, data = d, conf.type = "log")
sd <- summary(sf)
km <- data.frame(time = sd$time, n_risk = sd$n.risk, n_event = sd$n.event,
                 survival = sd$surv, lower = sd$lower, upper = sd$upper,
                 quadrant = sub("^quadrant=", "", as.character(sd$strata)), stringsAsFactors = FALSE)
km$quadrant <- factor(km$quadrant, levels = lev)
lr <- survdiff(Surv(OS_time_months, OS_status) ~ quadrant, data = d)
lr_p <- pchisq(lr$chisq, df = length(lr$n) - 1L, lower.tail = FALSE)
write_tsv(km, "Fig5E_KM_curve_source_data.tsv")
write_tsv(data.frame(n = nrow(d), events = sum(d$OS_status), chisq = unname(lr$chisq), df = length(lr$n) - 1L,
                     p_value = lr_p, test = "four-group log-rank", unit = "patient"), "Fig5E_statistics.tsv")
p5e <- ggplot(km, aes(time, survival, colour = quadrant)) +
  geom_step(linewidth = 0.65) +
  scale_colour_manual(values = pal, name = NULL) +
  scale_x_continuous(breaks = c(0, 100, 200)) +
  coord_cartesian(xlim = c(0, min(200, max(d$OS_time_months))), ylim = c(0, 1)) +
  annotate("text", x = 6, y = 0.10, hjust = 0, label = paste("Log-rank", fmt_p(lr_p)), size = 2.5, family = "Arial") +
  labs(title = "Four-group overall survival", x = "Time (months)", y = "Overall survival probability") +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE)) +
  theme_pub() + theme(legend.position = "bottom", legend.key.width = unit(5, "mm"), legend.text = element_text(size = 7.0))
save_pair(p5e, "Fig5E_CD79A_IKZF1_four_group_KM", 70, 65)

# Fig.5F: current frozen sparse4 distribution, not the retired score.
kw5f <- kruskal.test(sparse4 ~ quadrant, data = d)
pair_targets <- c("Cold", "B-skewed", "T-skewed")
pair5f <- do.call(rbind, lapply(pair_targets, function(g) {
  x <- d$sparse4[d$quadrant == "LymphHot"]
  y <- d$sparse4[d$quadrant == g]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.frame(comparison = paste("LymphHot vs", g), n_LymphHot = length(x), n_other = length(y),
             median_LymphHot = median(x), median_other = median(y), W = unname(wt$statistic), p_nominal = wt$p.value)
}))
pair5f$p_holm <- p.adjust(pair5f$p_nominal, method = "holm")
write_tsv(pair5f, "Fig5F_pairwise_Holm.tsv")
write_tsv(data.frame(n = nrow(d), statistic = unname(kw5f$statistic), df = unname(kw5f$parameter), p_value = kw5f$p.value,
                     test = "Kruskal-Wallis", multiplicity = "three LymphHot-versus-other contrasts use Holm correction"), "Fig5F_statistics.tsv")
write_tsv(d[, c("patient_id", "sample_id", "quadrant", "sparse4")], "Fig5F_source_data.tsv")
p5f <- ggplot(d, aes(quadrant, sparse4, fill = quadrant)) +
  geom_violin(trim = FALSE, colour = "grey35", linewidth = 0.3, alpha = 0.72) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  scale_fill_manual(values = pal, guide = "none") +
  annotate("text", x = 2.5, y = max(d$sparse4) + 0.03 * diff(range(d$sparse4)),
           label = paste("Kruskal-Wallis", fmt_p(kw5f$p.value)), size = 2.5, family = "Arial") +
  labs(title = "Sparse4 by expression group", x = NULL, y = "Frozen sparse4 linear predictor") +
  theme_pub() + theme(axis.text.x = element_text(angle = 32, hjust = 1))
save_pair(p5f, "Fig5F_sparse4_by_CD79A_IKZF1_groups", 60, 60)

# Fig.5G: within-CD79A-high descriptive immune context.
sub <- d[d$quadrant %in% c("B-skewed", "LymphHot"), ]
sub$quadrant <- factor(sub$quadrant, levels = c("B-skewed", "LymphHot"))
features <- c("checkpoint_score", "cytotoxic_score", "project_immune_ssGSEA", "TLS12")
feature_labels <- c("Checkpoint", "Cytotoxic", "Project immune", "TLS12")
long5g <- do.call(rbind, lapply(seq_along(features), function(i) data.frame(
  patient_id = sub$patient_id, quadrant = sub$quadrant, feature = feature_labels[i], value = sub[[features[i]]]
)))
long5g$feature <- factor(long5g$feature, levels = feature_labels)
st5g <- do.call(rbind, lapply(seq_along(features), function(i) {
  x <- sub[[features[i]]][sub$quadrant == "LymphHot"]
  y <- sub[[features[i]]][sub$quadrant == "B-skewed"]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.frame(feature = feature_labels[i], n_LymphHot = length(x), n_B_skewed = length(y),
             median_LymphHot = median(x), median_B_skewed = median(y), W = unname(wt$statistic), p_nominal = wt$p.value)
}))
st5g$p_holm <- p.adjust(st5g$p_nominal, method = "holm")
write_tsv(long5g, "Fig5G_source_data.tsv")
write_tsv(st5g, "Fig5G_statistics_Holm.tsv")
ann5g <- data.frame(feature = factor(st5g$feature, levels = feature_labels), quadrant = "B-skewed",
                    y = vapply(feature_labels, function(f) max(long5g$value[long5g$feature == f]), numeric(1)),
                    label = vapply(st5g$p_holm, fmt_p, character(1)))
p5g <- ggplot(long5g, aes(quadrant, value, fill = quadrant)) +
  geom_violin(trim = FALSE, alpha = 0.7, colour = "grey35", linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  geom_text(data = ann5g, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE,
            vjust = -0.35, size = 2.5, family = "Arial") +
  facet_wrap(~feature, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = pal[c("B-skewed", "LymphHot")], guide = "none") +
  labs(title = "Immune context within CD79A-high tumours", x = NULL, y = "Signature score") +
  theme_pub() + theme(axis.text.x = element_text(angle = 25, hjust = 1))
save_pair(p5g, "Fig5G_immune_features_CD79A_high", 180, 60)

# Fig.6A: median OS and model-based 95% CI from the same KM object as Fig.5E.
mt <- as.data.frame(summary(sf)$table)
mt$quadrant <- sub("^quadrant=", "", rownames(mt))
med6a <- data.frame(quadrant = factor(mt$quadrant, levels = lev), n = as.integer(mt$records), events = as.integer(mt$events),
                    median_months = mt$median, ci_lower = mt$`0.95LCL`, ci_upper = mt$`0.95UCL`, row.names = NULL)
write_tsv(med6a, "Fig6A_median_OS.tsv")
p6a_finite <- med6a[is.finite(med6a$ci_upper), ]
p6a_open <- med6a[!is.finite(med6a$ci_upper), ]
p6a <- ggplot(med6a, aes(quadrant, median_months, colour = quadrant)) +
  geom_errorbar(data = p6a_finite, aes(ymin = ci_lower, ymax = ci_upper), width = 0.18, linewidth = 0.6) +
  geom_segment(data = p6a_open, aes(x = quadrant, xend = quadrant, y = ci_lower, yend = 114),
               arrow = arrow(length = unit(1.4, "mm"), type = "open"), linewidth = 0.6) +
  geom_point(size = 2.5) +
  geom_text(aes(label = sprintf("%.1f", median_months)), vjust = -0.8, size = 2.5, family = "Arial", colour = "black") +
  scale_colour_manual(values = pal, guide = "none") +
  scale_x_discrete(limits = lev) +
  annotate("text", x = 2.5, y = 118, label = "Arrow: upper 95% CI not reached", size = 2.50, family = "Arial") +
  coord_cartesian(ylim = c(15, 120)) +
  labs(title = "Median overall survival", x = NULL, y = "Median overall survival (months)") +
  theme_pub() + theme(axis.text.x = element_text(angle = 32, hjust = 1))
save_pair(p6a, "Fig6A_median_OS_by_group", 60, 60)

# Fig.6B: checkpoint-expression heatmap; inferential family is six gene-wise KW tests.
gene_z <- checkpoint_z
mean6b <- do.call(rbind, lapply(checkpoint_genes, function(g) do.call(rbind, lapply(lev, function(q) {
  idx <- d$quadrant == q
  data.frame(gene = g, quadrant = q, n = sum(idx), mean_z = mean(gene_z[g, idx]), median_z = median(gene_z[g, idx]))
}))))
test6b <- do.call(rbind, lapply(checkpoint_genes, function(g) {
  kt <- kruskal.test(as.numeric(gene_z[g, ]) ~ d$quadrant)
  data.frame(gene = g, statistic = unname(kt$statistic), df = unname(kt$parameter), p_nominal = kt$p.value)
}))
test6b$p_BH <- p.adjust(test6b$p_nominal, method = "BH")
mean6b$gene <- factor(mean6b$gene, levels = rev(checkpoint_genes))
mean6b$quadrant <- factor(mean6b$quadrant, levels = lev)
write_tsv(mean6b, "Fig6B_group_mean_checkpoint_z.tsv")
write_tsv(test6b, "Fig6B_gene_KW_BH.tsv")
p6b <- ggplot(mean6b, aes(quadrant, gene, fill = mean_z)) +
  geom_tile(colour = "white", linewidth = 0.5) +
  geom_text(aes(label = sprintf("%.2f", mean_z)), size = 2.5, family = "Arial") +
  scale_fill_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0, name = "Mean z") +
  labs(title = "Checkpoint expression", x = NULL, y = NULL) +
  theme_pub() + theme(axis.text.x = element_text(angle = 32, hjust = 1), legend.position = "right")
save_pair(p6b, "Fig6B_checkpoint_expression_heatmap", 60, 60)

# Fig.6C: project-defined immune ssGSEA score across groups.
kw6c <- kruskal.test(project_immune_ssGSEA ~ quadrant, data = d)
pw <- combn(lev, 2, simplify = FALSE)
pair6c <- do.call(rbind, lapply(pw, function(z) {
  x <- d$project_immune_ssGSEA[d$quadrant == z[1]]
  y <- d$project_immune_ssGSEA[d$quadrant == z[2]]
  wt <- wilcox.test(x, y, exact = FALSE)
  data.frame(comparison = paste(z, collapse = " vs "), n_1 = length(x), n_2 = length(y), W = unname(wt$statistic), p_nominal = wt$p.value)
}))
pair6c$p_holm <- p.adjust(pair6c$p_nominal, method = "holm")
write_tsv(d[, c("patient_id", "sample_id", "quadrant", "project_immune_ssGSEA")], "Fig6C_source_data.tsv")
write_tsv(data.frame(n = nrow(d), statistic = unname(kw6c$statistic), df = unname(kw6c$parameter), p_value = kw6c$p.value,
                     test = "Kruskal-Wallis", score_identity = "project-defined immune ssGSEA; not official ESTIMATE"), "Fig6C_statistics.tsv")
write_tsv(pair6c, "Fig6C_pairwise_Holm.tsv")
p6c <- ggplot(d, aes(quadrant, project_immune_ssGSEA, fill = quadrant)) +
  geom_violin(trim = FALSE, alpha = 0.72, colour = "grey35", linewidth = 0.3) +
  geom_boxplot(width = 0.16, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  scale_fill_manual(values = pal, guide = "none") +
  annotate("text", x = 2.5, y = max(d$project_immune_ssGSEA), vjust = -0.35,
           label = paste("Kruskal-Wallis", fmt_p(kw6c$p.value)), size = 2.5, family = "Arial") +
  labs(title = "Project immune score", x = NULL, y = "Project immune ssGSEA score") +
  theme_pub() + theme(axis.text.x = element_text(angle = 32, hjust = 1))
save_pair(p6c, "Fig6C_project_immune_score_by_group", 60, 60)

# Geometry and explicit scientific disposition.
geometry <- data.frame(
  panel = c("Fig5B", "Fig5E", "Fig5F", "Fig5G", "Fig6A", "Fig6B", "Fig6C"),
  width_mm = c(60, 60, 60, 180, 60, 60, 60), height_mm = rep(60, 7),
  geometry_basis = c(rep("recoverable Fig.5 production theme", 4), rep("recoverable Fig.6 sibling production contract", 3)),
  stringsAsFactors = FALSE
)
write_tsv(geometry, "FIGURE_GEOMETRY.tsv")

disp <- data.frame(
  panel = geometry$panel,
  disposition = "RECOMPUTED_PENDING_DIFFERENT_ROLE_AUDIT",
  scientific_task = c("patient-level bulk CD79A-IKZF1 association", "four-group descriptive OS KM",
                      "current sparse4 distribution by independent marker-defined groups",
                      "immune context within CD79A-high tumours", "four-group median OS with 95% CI",
                      "six checkpoint genes across groups", "project immune ssGSEA across groups"),
  independent_unit = "patient",
  claim_ceiling = c("same-cohort association; shared immune abundance remains an alternative",
                    "same-cohort descriptive grouping; not subtype validation",
                    "association only; groups do not derive the sparse4 model",
                    "same-cohort descriptive context; not independent of immune abundance",
                    "descriptive medians from the Fig.5E KM object",
                    "descriptive bulk expression context; not treatment prediction",
                    "descriptive immune abundance; not an independent phenotype or purity measure"),
  stringsAsFactors = FALSE
)
write.table(disp, file.path(MOD, "PANEL_DISPOSITION.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(LOG, "SESSIONINFO.txt"))
cat("Completed seven-panel recomputation.\n")
