#!/usr/bin/env Rscript

# ==========================================================
# Script   : modified_render_all_replacement_figures.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:replacement_figures_20261018/scripts/modified/render_all_replacement_figures.R
# ==========================================================

# Post-freeze display integration adapted from the latest dated scripts copied
# from HNSCC/final. This script does not fit or tune a model and does not run
# bootstrap analyses. Every displayed estimate is checked against a frozen table.

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ragg)
  library(jsonlite)
  library(ggrepel)
})

args0 <- commandArgs(trailingOnly = FALSE)
script_path <- normalizePath(sub("^--file=", "", grep("^--file=", args0, value = TRUE)[1]))
ROOT <- normalizePath(file.path(dirname(script_path), "../../.."))
BASE <- file.path(ROOT, "replacement_figures_20261018")
FIG <- file.path(BASE, "figures")
SRC <- file.path(BASE, "source_data")
LOG <- file.path(BASE, "logs")
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG, recursive = TRUE, showWarnings = FALSE)

read_tsv <- function(path) read.delim(file.path(ROOT, path), check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) {
  # Plot labels may contain line breaks, but exported TSV fields must remain
  # one physical line per record for robust downstream parsing.
  y <- x
  for (j in seq_along(y)) {
    if (is.character(y[[j]]) || is.factor(y[[j]])) {
      y[[j]] <- gsub("[\r\n]+", " / ", as.character(y[[j]]))
    }
  }
  write.table(y, file.path(SRC, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
stopif <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
close_to <- function(x, y, tol = 1e-10) isTRUE(all.equal(as.numeric(x), as.numeric(y), tolerance = tol, check.attributes = FALSE))

theme_pub <- function() {
  theme_classic(base_family = "Arial", base_size = 8) +
    theme(
      plot.title = element_text(size = 8.5, face = "bold", hjust = 0.5, lineheight = 0.95),
      plot.subtitle = element_text(size = 7.5, hjust = 0.5, lineheight = 0.95),
      plot.caption = element_text(size = 7.5, colour = "grey25", hjust = 0, lineheight = 0.95),
      axis.title = element_text(size = 8), axis.text = element_text(size = 7.5, colour = "black"),
      legend.title = element_text(size = 7.5), legend.text = element_text(size = 7.5),
      plot.margin = margin(4, 5, 4, 5)
    )
}

save_plot <- function(p, name, w_mm, h_mm) {
  pdf_path <- file.path(FIG, paste0(name, ".pdf"))
  png_path <- file.path(FIG, paste0(name, ".png"))
  grDevices::cairo_pdf(pdf_path, width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
  print(p); grDevices::dev.off()
  ragg::agg_png(png_path, width = w_mm / 25.4, height = h_mm / 25.4,
                units = "in", res = 600, background = "white")
  print(p); grDevices::dev.off()
}

save_surv <- function(g, name, w_mm, h_mm) {
  pdf_path <- file.path(FIG, paste0(name, ".pdf"))
  png_path <- file.path(FIG, paste0(name, ".png"))
  grDevices::cairo_pdf(pdf_path, width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
  print(g, newpage = FALSE); grDevices::dev.off()
  ragg::agg_png(png_path, width = w_mm / 25.4, height = h_mm / 25.4,
                units = "in", res = 600, background = "white")
  print(g, newpage = FALSE); grDevices::dev.off()
}

# Authorities and display-level patient data.
lock <- jsonlite::fromJSON(file.path(ROOT, "results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/MODEL_LOCK.json"))
coef <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/COEFFICIENTS.tsv")
score <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/SCORE_TABLE.tsv")
clinical <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/CLINICAL_MASTER_476.tsv")
stopif(nrow(score) == 476 && nrow(clinical) == 476, "Frozen TCGA roster must contain 476 patients")
stopif(identical(score$sample_id, clinical$sample_id), "Score/clinical order mismatch")
stopif(max(abs(score$linear_predictor - clinical$score_raw_linear_predictor)) < 1e-12, "Score authority mismatch")

vst <- readRDS(file.path(ROOT, "inputs/tcga/expression/vst_matrix.rds"))
mi <- match(score$sample_id, colnames(vst))
stopif(!anyNA(mi) && !anyDuplicated(mi), "VST-to-score sample mapping failed")
vst476 <- vst[, mi, drop = FALSE]
zscore <- function(x) as.numeric((x - mean(x)) / sd(x))

# Fig. 2E: real-data description of all eligible Yellow-module candidates.
# Selection frequency is descriptive output from the same nested-CV procedure;
# it was not used as an additional screening step for the frozen full-data fit.
prov <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/CANDIDATE_PROVENANCE.tsv")
selection_frequency <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_SELECTION_FREQUENCY.tsv")
nested_coef <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_COEFFICIENTS.tsv")
nested_tasks <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_TASK_METADATA.tsv")
hub <- read.csv(file.path(ROOT, "inputs/wgcna/legacy/hub_genes_all.csv"), check.names = FALSE, stringsAsFactors = FALSE)
stopif(nrow(prov) == 539 && sum(prov$eligible_for_model) == 492 && nrow(coef) == 4, "Candidate flow counts drift")
eligible <- prov$gene_symbol[prov$eligible_for_model]
stopif(nrow(selection_frequency) == 492 && !anyDuplicated(selection_frequency$gene_symbol), "Nested selection-frequency roster drift")
stopif(nrow(nested_tasks) == 100 && length(unique(nested_tasks$repeat_id)) == 20 &&
         length(unique(nested_tasks$outer_fold)) == 5, "Nested task roster drift")
stopif(all(is.finite(nested_tasks$heldout_harrell_c)), "Nested task failure/nonfinite held-out C")
stopif(sum(nested_tasks$selected_n) == nrow(nested_coef), "Nested nonzero-coefficient accounting drift")
stopif(setequal(selection_frequency$gene_symbol, eligible), "Nested candidate roster differs from frozen eligible roster")
yellow_hub <- hub[hub$module == "yellow", c("gene", "kME")]
fig2e <- merge(selection_frequency, yellow_hub, by.x = "gene_symbol", by.y = "gene", all.x = TRUE, sort = FALSE)
fig2e <- fig2e[match(selection_frequency$gene_symbol, fig2e$gene_symbol), ]
stopif(!anyNA(fig2e$kME), "Fig2E candidate-to-kME mapping failed")
fig2e$frozen_gene <- fig2e$gene_symbol %in% coef$gene_symbol
write_tsv(fig2e, "Fig2E_candidate_kME_selection_frequency.tsv")
write_tsv(prov[!prov$eligible_for_model, , drop = FALSE], "Fig2E_47_excluded_candidates.tsv")
write_tsv(coef, "Fig2E_model_genes_and_coefficients.tsv")
p2e <- ggplot(fig2e, aes(kME, selection_frequency)) +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = fig2e[!fig2e$frozen_gene, ], size = 0.9, alpha = 0.45, colour = "grey45") +
  geom_point(data = fig2e[fig2e$frozen_gene, ], size = 2.1, shape = 21,
             fill = "#B2182B", colour = "white", stroke = 0.35) +
  ggrepel::geom_text_repel(data = fig2e[fig2e$frozen_gene, ], aes(label = gene_symbol),
                           family = "Arial", size = 2.65, fontface = "bold", colour = "#8B1A1A",
                           box.padding = 0.22, point.padding = 0.18, min.segment.length = 0,
                           segment.size = 0.25, seed = 20261019, max.overlaps = Inf) +
  scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, .2), expand = expansion(mult = c(0, .04))) +
  labs(title = "Candidate-gene selection stability",
       x = "Yellow-module membership (kME)", y = "Outer-model selection frequency") + theme_pub()
save_plot(p2e, "Fig2E_model_development_flow_sparse4", 70, 50)

# Fig. 2F: frozen full-fit inner selection surface; no re-fitting.
surface <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_FULL_SELECTION_SURFACE.tsv")
stopif(all(is.finite(surface$lambda)) && all(surface$lambda > 0), "Selection surface contains non-positive lambda")
surface$alpha <- factor(surface$alpha)
selected <- surface[surface$selected %in% c(TRUE, "TRUE"), , drop = FALSE]
stopif(nrow(selected) == 1 && close_to(selected$lambda, lock$selection_rule$selected_lambda, 1e-12), "Selected lambda drift")
write_tsv(surface, "Fig2F_frozen_selection_surface.tsv")
selection_summary <- data.frame(
  alpha = lock$selection_rule$selected_alpha,
  lambda = lock$selection_rule$selected_lambda,
  nonzero_genes = nrow(coef),
  inner_folds = lock$selection_rule$inner_folds,
  rule_name = lock$selection_rule$rule_name,
  band_definition = "cvm +/- cvsd",
  stringsAsFactors = FALSE
)
write_tsv(selection_summary, "Fig2F_selection_rule_summary.tsv")
p2f <- ggplot(surface, aes(log(lambda), cvm, colour = alpha, group = alpha)) +
  geom_ribbon(aes(ymin = cvlo, ymax = cvup, fill = alpha), alpha = 0.08, colour = NA) +
  geom_line(linewidth = 0.55) +
  geom_vline(xintercept = log(selected$lambda), linetype = "dotted", linewidth = 0.4, colour = "grey35") +
  geom_point(data = selected, shape = 21, size = 2.2, stroke = 0.6, fill = "white") +
  scale_colour_manual(values = c("0.5" = "#2166AC", "0.75" = "#7B3294", "1" = "#B2182B")) +
  scale_fill_manual(values = c("0.5" = "#2166AC", "0.75" = "#7B3294", "1" = "#B2182B")) +
  labs(title = "Elastic-net model selection",
       subtitle = sprintf("alpha=%.2f; lambda=%.4f; %d genes",
                          lock$selection_rule$selected_alpha, lock$selection_rule$selected_lambda, nrow(coef)),
       x = "log(lambda)", y = "Inner Harrell C", colour = "alpha", fill = "alpha") +
  guides(colour = guide_legend(title = NULL, nrow = 1), fill = "none") +
  theme_pub() + theme(legend.position = "top", legend.direction = "horizontal",
                      legend.key.width = grid::unit(4, "mm"), legend.spacing.x = grid::unit(1, "mm"),
                      legend.background = element_blank(), legend.margin = margin(0, 0, 0, 0),
                      legend.box.spacing = grid::unit(0, "pt"))
save_plot(p2f, "Fig2F_frozen_elastic_net_selection", 60, 50)

# Fig. 2G: CD79A versus current frozen score.
cd79a_tab <- read_tsv("agent_work/sparse4_bio_bulk_20261016/results/BIO_BULK_SPARSE4_V3/CD79A_SCORE_ASSOCIATION.tsv")
scatter_cd79a <- data.frame(patient_id = score$patient_id, CD79A_VST = as.numeric(vst476["CD79A", ]), score = score$linear_predictor)
rho_cd79a <- cor(scatter_cd79a$CD79A_VST, scatter_cd79a$score, method = "spearman")
stopif(close_to(rho_cd79a, cd79a_tab$spearman_rho), "CD79A-score rho mismatch")
write_tsv(scatter_cd79a, "Fig2G_CD79A_sparse4_display_data.tsv")
p2g <- ggplot(scatter_cd79a, aes(CD79A_VST, score)) +
  geom_point(size = 0.8, alpha = 0.48, colour = "#2166AC") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20", fill = "grey70", linewidth = 0.5, alpha = 0.25) +
  labs(title = "CD79A and B-cell score", x = "CD79A expression (VST)", y = "B-cell expression score") +
  theme_pub() + theme(plot.margin = margin(8, 5, 4, 5))
save_plot(p2g, "Fig2G_CD79A_vs_sparse4_score", 60, 50)

# Shared descriptive median-split Kaplan-Meier display.
km_effect <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/MEDIAN_KM_DESCRIPTIVE_EFFECT.tsv")
km_group <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/MEDIAN_KM_GROUP_SUMMARY.tsv")
km <- clinical[, c("patient_id", "OS_time_months", "OS_status", "score_median_group")]
km$score_median_group <- factor(km$score_median_group, levels = c("Low", "High"))
stopif(identical(as.integer(table(km$score_median_group)), c(238L, 238L)), "KM groups drift")
fit_km <- survfit(Surv(OS_time_months, OS_status) ~ score_median_group, data = km)
write_tsv(km, "Fig2H_Bcell_score_KM_display_data.tsv")
make_km <- function(title) {
  g <- ggsurvplot(
    fit_km, data = km, risk.table = TRUE, risk.table.height = 0.30, conf.int = FALSE,
    palette = c("#2166AC", "#B2182B"), censor.shape = 124, censor.size = 2.2,
    break.time.by = 24, xlim = c(0, 120),
    legend.title = NULL, legend.labs = c("Low", "High"),
    title = title, xlab = "Overall survival (months)", ylab = "Survival probability",
    pval = sprintf("Log-rank p = %.2e\nMedian cutoff = %.4f", km_effect$logrank_p, km_effect$cutoff_raw_global_median),
    pval.coord = c(30, 0.94),
    pval.size = 2.65, fontsize = 2.65, risk.table.fontsize = 2.65,
    ggtheme = theme_pub() + theme(legend.position = c(0.72, 0.82)),
    tables.theme = theme_cleantable(base_family = "Arial", base_size = 7.5)
  )
  # The risk-table scale is drawn bottom-to-top as High, Low.
  g$plot <- g$plot + theme(legend.title = element_blank())
  g$table <- g$table +
    scale_y_discrete(labels = c("H", "L")) +
    scale_x_continuous(limits = c(0, 120), breaks = seq(0, 120, 24),
                       expand = expansion(mult = c(0.12, 0.04))) +
    labs(x = "Time (months)") +
    theme(axis.text.x = element_text(size = 7.0, colour = "black"),
          axis.ticks.x = element_line(linewidth = 0.3, colour = "black"),
          axis.ticks.length.x = grid::unit(1.0, "mm"),
          axis.line.x = element_line(linewidth = 0.3, colour = "black"),
          axis.title.x = element_text(size = 7.5, colour = "black", margin = margin(t = 2)))
  g
}
save_surv(make_km("B-cell score and overall survival"), "Fig2H_sparse4_KM", 60, 80)

# Fig. 4C: replace unsupported subgroup-specific HRs with the frozen formal
# score-by-treatment interaction estimate.
tx <- read_tsv("agent_work/analytical_gap_wave1_repair_v2_execution_20261016/formal_run_001/dg/G_TREATMENT_INTERACTION_V2.tsv")
stopif(all(c(tx$interaction_hazard_ratio, tx$wald_ci_lower, tx$wald_ci_upper) > 0), "Treatment HR scale is non-positive")
write_tsv(tx, "Fig4C_treatment_interaction.tsv")
p4c <- ggplot(tx, aes(y = "Surgery + adjuvant\nvs non-surgical", x = interaction_hazard_ratio,
                      xmin = wald_ci_lower, xmax = wald_ci_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.16, linewidth = 0.55, colour = "#2166AC") +
  geom_point(size = 2.5, shape = 21, fill = "#2166AC", colour = "white", stroke = 0.45) +
  annotate("text", x = 0.95, y = 1.22, hjust = 0.5,
           label = sprintf("HR %.2f [%.2f, %.2f]\nLR p = %.3f\nn = %d; events = %d",
                           tx$interaction_hazard_ratio, tx$wald_ci_lower, tx$wald_ci_upper,
                           tx$p_lr, tx$n, tx$events), family = "Arial", size = 2.65) +
  scale_x_log10(limits = c(0.55, 1.55), breaks = c(0.6, 0.8, 1, 1.2, 1.5)) +
  labs(title = "Score-by-treatment interaction", x = "Interaction HR (log scale)", y = NULL,
       caption = "Per 1-SD score; adjusted for treatment main effects\nObservational interaction; no treatment-effect claim") + theme_pub()
save_plot(p4c, "Fig4C_sparse4_treatment_interaction", 90, 60)

# Fig. 4D: common-roster M1-M4 score estimates.
m14 <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/M1_M4_COMMON_ROSTER_COEFFICIENTS.tsv")
m14 <- m14[m14$term == "score_per_global_sd", ]
m14_label_map <- c(M1 = "M1: score", M2 = "M2: + immune", M3 = "M3: + age and stage", M4 = "M4: + HPV, site,\nsurgery, stroma")
m14$model_label <- unname(m14_label_map[m14$model_id])
m14$model_label <- factor(m14$model_label, levels = rev(m14$model_label))
ph <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/PH_DIAGNOSTICS_COMMON_ROSTER.tsv")
m4_ph <- ph[ph$model_id == "M4" & ph$term == "score_per_global_sd", , drop = FALSE]
m4_global_ph <- ph[ph$model_id == "M4" & ph$term == "GLOBAL", , drop = FALSE]
stopif(nrow(m14) == 4 && all(m14$n == 320) && all(m14$events == 106), "M1-M4 common roster drift")
stopif(nrow(m4_ph) == 1 && nrow(m4_global_ph) == 1, "M4 PH diagnostic missing")
write_tsv(m14, "Fig4D_common_roster_M1_M4_score.tsv")
p4d <- ggplot(m14, aes(hazard_ratio, model_label, xmin = ci_lower, xmax = ci_upper)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.18, linewidth = 0.55, colour = "#2166AC") +
  geom_point(shape = 21, size = 2.5, fill = "#2166AC", colour = "white", stroke = 0.45) +
  geom_text(aes(x = hazard_ratio,
                label = sprintf("%.2f [%.2f, %.2f]", hazard_ratio, ci_lower, ci_upper)),
            nudge_y = 0.22, hjust = 0.5, size = 2.65, family = "Arial") +
  scale_x_log10(limits = c(0.9, 2.45), breaks = c(1, 1.25, 1.5, 2, 2.4)) +
  labs(title = "Sequential adjusted score models",
       x = "Hazard ratio per score SD", y = NULL) + theme_pub()
save_plot(p4d, "Fig4D_sparse4_M1_M4_common_roster", 90, 60)

# Fig. 4E: a transparent changed-estimand replacement using analytic paired-C
# uncertainty for clinical-base versus clinical-plus-score models.
inc <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/CLINICAL_INCREMENTAL_PAIRED_C_NONBOOTSTRAP.tsv")
inc$se <- sqrt(inc$variance_delta); inc$ci_lower <- inc$delta_c - 1.96 * inc$se; inc$ci_upper <- inc$delta_c + 1.96 * inc$se
inc$label <- sprintf("%s\nHolm p = %.3f",
                     c("Clinical base + score", "Extended base + score"),
                     inc$p_holm_C_INCREMENTAL_C_2)
write_tsv(inc, "Fig4E_incremental_paired_C_nonbootstrap.tsv")
p4e <- ggplot(inc, aes(delta_c, factor(label, levels = rev(label)), xmin = ci_lower, xmax = ci_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.18, linewidth = 0.55, colour = "#2166AC") +
  geom_point(shape = 21, size = 2.6, fill = "#2166AC", colour = "white", stroke = 0.45) +
  geom_text(aes(x = delta_c,
                label = sprintf("Delta C %.3f [%.3f, %.3f]", delta_c, ci_lower, ci_upper)),
            nudge_y = 0.22, hjust = 0.5, family = "Arial", size = 2.65) +
  scale_x_continuous(limits = c(-0.015, 0.078), breaks = c(0, 0.02, 0.04, 0.06)) +
  labs(title = "Incremental score discrimination",
       x = "Paired delta Harrell C (n = 320; events = 106)", y = NULL) + theme_pub()
save_plot(p4e, "Fig4E_sparse4_incremental_delta_C", 90, 60)

# Fig. 4F: current prespecified TLS12 transcriptional score versus frozen score.
tls_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
stopif(all(tls_genes %in% rownames(vst476)), "TLS12 coverage failed")
tls_score <- colMeans(t(apply(vst476[tls_genes, , drop = FALSE], 1, zscore)))
tls_frozen <- read_tsv("agent_work/sparse4_bio_bulk_20261016/results/BIO_BULK_SPARSE4_V3/TLS12_SCORE_ASSOCIATION.tsv")
tls_display <- data.frame(patient_id = score$patient_id, score = score$linear_predictor, TLS12 = tls_score)
rho_tls <- cor(tls_display$score, tls_display$TLS12, method = "spearman")
stopif(close_to(rho_tls, tls_frozen$spearman_rho), "TLS12-score rho mismatch")
write_tsv(tls_display, "Fig4F_TLS12_sparse4_display_data.tsv")
p4f <- ggplot(tls_display, aes(score, TLS12)) +
  geom_point(size = 0.8, alpha = 0.48, colour = "#2166AC") +
  geom_smooth(method = "lm", se = TRUE, colour = "grey20", fill = "grey70", linewidth = 0.5, alpha = 0.25) +
  annotate("label", x = Inf, y = Inf, hjust = 1.05, vjust = 1.15,
           label = sprintf("Spearman rho = %.3f\np = %.2e; n = 476", rho_tls, tls_frozen$p_value),
           family = "Arial", size = 2.65, linewidth = 0.2) +
  labs(title = "B-cell score and TLS12",
       x = "B-cell expression score", y = "TLS12 chemokine score") + theme_pub()
save_plot(p4f, "Fig4F_TLS12_vs_sparse4_score", 60, 60)

# Fig. 4B: exploratory responder/non-responder display using the exact same
# frozen GSE159067 two-gene partial score and response labels as Fig. 6E.
external_ledger <- read_tsv("agent_work/yellow_sparse_enet_external_20261016/EXTERNAL_EVALUATION/EXTERNAL_SAMPLE_ENDPOINT_SCORE_LEDGER.tsv")
g159_meta <- read.csv(file.path(ROOT, "inputs/external/legacy_derivatives/gse159067_ies.csv"), check.names = FALSE, stringsAsFactors = FALSE)
g159_score <- external_ledger[external_ledger$cohort == "GSE159067", c("sample_id", "score")]
stopif(nrow(g159_score) == 102 && identical(g159_score$sample_id, g159_meta$sample_id), "GSE159067 score/response mapping drift")
response_raw <- toupper(trimws(g159_meta$response))
stopif(all(response_raw %in% c("CR", "PR", "SD", "PD")), "Unexpected GSE159067 response category")
g159_box <- data.frame(
  sample_id = g159_score$sample_id,
  partial_score = g159_score$score,
  response = response_raw,
  responder = response_raw %in% c("CR", "PR"),
  stringsAsFactors = FALSE
)
g159_box$response_group <- factor(ifelse(g159_box$responder, "Responder", "Non-responder"),
                                  levels = c("Responder", "Non-responder"))
g159_wilcox <- wilcox.test(partial_score ~ response_group, data = g159_box, exact = FALSE)
stopif(sum(g159_box$responder) == 11 && sum(!g159_box$responder) == 91, "GSE159067 response counts drift")
write_tsv(g159_box, "Fig4B_GSE159067_partial_score_by_response.tsv")
write_tsv(data.frame(test = "two-sided Wilcoxon rank-sum with continuity correction",
                     n = nrow(g159_box), responders = sum(g159_box$responder),
                     nonresponders = sum(!g159_box$responder), statistic_W = unname(g159_wilcox$statistic),
                     p_value = g159_wilcox$p.value,
                     claim_boundary = "exploratory two-gene partial score; not complete-model validation"),
          "Fig4B_GSE159067_partial_score_response_test.tsv")
fig4b_labels <- c("Responder" = "Resp.\n(n=11)",
                  "Non-responder" = "Non-resp.\n(n=91)")
p4b <- ggplot(g159_box, aes(response_group, partial_score, fill = response_group, colour = response_group)) +
  geom_boxplot(width = 0.48, alpha = 0.30, outlier.shape = NA, linewidth = 0.4) +
  geom_jitter(width = 0.11, size = 0.85, alpha = 0.60) +
  scale_x_discrete(labels = fig4b_labels) +
  scale_fill_manual(values = c("Responder" = "#2166AC", "Non-responder" = "#B2182B")) +
  scale_colour_manual(values = c("Responder" = "#2166AC", "Non-responder" = "#B2182B")) +
  labs(title = "ICI response and\npartial score", x = NULL, y = "Two-gene partial score",
       subtitle = sprintf("Wilcoxon p = %.3f", g159_wilcox$p.value)) +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(lineheight = 0.9))
save_plot(p4b, "Fig4B_GSE159067_partial_score_by_response", 50, 60)

# Fig. 5F: current score across retained CD79A-by-IKZF1 quadrants.
cd <- as.numeric(vst476["CD79A", ]); ik <- as.numeric(vst476["IKZF1", ])
quad <- ifelse(cd >= median(cd) & ik >= median(ik), "LymphHot",
        ifelse(cd >= median(cd) & ik < median(ik), "B-skewed",
        ifelse(cd < median(cd) & ik >= median(ik), "T-skewed", "Cold")))
qdat <- data.frame(patient_id = score$patient_id, quadrant = factor(quad, levels = c("Cold", "B-skewed", "T-skewed", "LymphHot")), score = score$linear_predictor)
qsum <- read_tsv("agent_work/sparse4_bio_bulk_20261016/results/BIO_BULK_SPARSE4_V3/LYMPHHOT_QUADRANT_SCORE_SUMMARY.tsv")
obs <- aggregate(score ~ quadrant, qdat, function(x) c(n = length(x), median = median(x), q1 = quantile(x, .25), q3 = quantile(x, .75)))
stopif(all(sort(as.integer(table(qdat$quadrant))) == sort(qsum$n)), "Quadrant counts drift")
stopif(max(abs(sort(obs$score[, "median"]) - sort(qsum$median_score))) < 1e-12, "Quadrant medians drift")
write_tsv(qdat, "Fig5F_sparse4_by_quadrant_display_data.tsv")
q_label_name <- c("Cold" = "Cold", "B-skewed" = "B-\nskewed", "T-skewed" = "T-\nskewed", "LymphHot" = "Lymph-\nHot")
q_labels <- setNames(sprintf("%s\n(n=%d)", q_label_name[qsum$quadrant], qsum$n), qsum$quadrant)
p5f <- ggplot(qdat, aes(quadrant, score, fill = quadrant)) +
  geom_violin(trim = TRUE, scale = "width", alpha = 0.72, linewidth = 0.25) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.35) +
  scale_fill_manual(values = c("Cold" = "#BDBDBD", "B-skewed" = "#80B1D3", "T-skewed" = "#FDB462", "LymphHot" = "#8DD3C7")) +
  scale_x_discrete(labels = q_labels) +
  labs(title = "Score across marker states", x = NULL, y = "B-cell expression score",
       subtitle = sprintf("Kruskal-Wallis p = %.2e", qsum$omnibus_kruskal_p[1]),
       caption = "Median CD79A/IKZF1 splits\nOmnibus test only") +
  theme_pub() + theme(legend.position = "none", axis.text.x = element_text(angle = 0, hjust = 0.5, lineheight = 0.9))
save_plot(p5f, "Fig5F_sparse4_by_CD79A_IKZF1_quadrant", 60, 60)

# Fig. 6D: current adjusted marker-comparison delta-C estimates.
marker_all <- read_tsv("agent_work/analytical_gap_wave1_20261016/formal_run_006/results/C_MARKER_PAIRED_DELTA_C.tsv")
marker_all$displayed_in_panel <- marker_all$model_1 == "SCORE" & marker_all$model_2 %in% c("CD79A", "MS4A1")
stopif(nrow(marker_all) == 3 && sum(marker_all$displayed_in_panel) == 2, "Fig6D three-comparison Holm family drift")
write_tsv(marker_all, "Fig6D_score_marker_paired_delta_C.tsv")
marker <- marker_all[marker_all$displayed_in_panel, ]
marker$label <- sprintf("Extended + score vs\nextended + %s\nHolm p = %.3f", marker$model_2, marker$p_adjusted_Holm)
p6d <- ggplot(marker, aes(delta_c, factor(label, levels = rev(label)), xmin = ci_lower, xmax = ci_upper)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.18, linewidth = 0.55, colour = "#2166AC") +
  geom_point(shape = 21, size = 2.6, fill = "#2166AC", colour = "white", stroke = 0.45) +
  geom_label(aes(x = 0.049,
                 label = sprintf("Delta C %.3f\nCI [%.3f, %.3f]\nC %.3f vs %.3f", delta_c, ci_lower, ci_upper, c_model_1, c_model_2)),
             nudge_y = 0.36, hjust = 1, family = "Arial", size = 2.65, lineheight = 0.9,
             linewidth = 0.2, label.padding = grid::unit(0.55, "mm"), fill = "white") +
  scale_x_continuous(limits = c(-0.008, 0.052), breaks = c(0, 0.02, 0.04)) +
  labs(title = "Score vs marker models",
       x = "Paired delta Harrell C", y = NULL) + theme_pub()
save_plot(p6d, "Fig6D_sparse4_vs_markers_delta_C", 60, 60)

# Fig. 6E: current exploratory partial-score response AUC. The historical
# responder/non-responder distribution (Fig. 4B) is not regenerated.
resp <- read_tsv("agent_work/yellow_sparse_enet_external_20261016/EXTERNAL_EVALUATION/EXTERNAL_RESPONSE_RESULTS.tsv")
coverage <- read_tsv("agent_work/yellow_sparse_enet_external_20261016/EXTERNAL_GENE_COVERAGE.tsv")
g159 <- coverage[coverage$cohort == "GSE159067" & coverage$eligible_feature %in% c(TRUE, "TRUE"), , drop = FALSE]
stopif(identical(g159$gene_symbol, c("CAMK2N1", "MS4A1")), "GSE159067 partial-score gene identity drift")
write_tsv(resp, "Fig6E_partial_score_response_AUC.tsv")
write_tsv(g159[, c("cohort", "gene_symbol", "selected_feature", "gamma_external_z_weight")], "Fig6E_partial_score_formula.tsv")
p6e <- ggplot(resp, aes(x = "GSE159067\nCAMK2N1 + MS4A1", y = AUC, ymin = CI95_lower, ymax = CI95_upper)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey55") +
  geom_errorbar(width = 0.12, linewidth = 0.55, colour = "#B2182B") +
  geom_point(shape = 21, size = 3, fill = "#B2182B", colour = "white", stroke = 0.45) +
  scale_y_continuous(limits = c(0.42, 1.00), breaks = seq(0.5, 1.0, 0.1)) +
  labs(title = "ICI-response discrimination", x = NULL, y = "AUC",
       subtitle = sprintf("AUC %.3f [%.3f, %.3f]\nn=%d; responders=%d",
                          resp$AUC, resp$CI95_lower, resp$CI95_upper, resp$n, resp$responders),
       caption = "Exploratory 2/4-gene partial score\nDeLong CI; within-cohort z scores\nTCGA-derived weights") + theme_pub()
save_plot(p6e, "Fig6E_GSE159067_partial_score_AUC", 60, 60)

# Supplementary Fig. 2E: current cohort-wise survival associations.
ext <- read_tsv("agent_work/yellow_sparse_enet_external_20261016/EXTERNAL_EVALUATION/EXTERNAL_SURVIVAL_RESULTS.tsv")
qfull <- read_tsv("agent_work/yellow_sparse_external_independent_audit_20261016/INDEPENDENT_FULL_OS_MULTIPLICITY.tsv")
ext$BH_q <- qfull$BH_q_across_two_FULL_OS_tests[match(ext$cohort, qfull$cohort)]
ext$label <- sprintf("%s\n%s %d/%d; n=%d, e=%d", ext$cohort, ext$coverage_class, ext$score_gene_n, nrow(coef), ext$n, ext$events)
external_order <- c("GSE41613", "GSE65858", "GSE159067")
ext$label <- factor(ext$label, levels = rev(ext$label[match(external_order, ext$cohort)]))
write_tsv(ext, "SFig2E_external_survival_forest.tsv")
pse <- ggplot(ext, aes(HR_per_score_SD, label, xmin = CI95_lower, xmax = CI95_upper, colour = coverage_class)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey55") +
  geom_errorbarh(height = 0.18, linewidth = 0.6) + geom_point(size = 2.6) +
  geom_text(aes(x = 2.28,
                label = ifelse(coverage_class == "FULL",
                               sprintf("HR %.2f [%.2f, %.2f]\np=%.3f; BH q=%.3f", HR_per_score_SD, CI95_lower, CI95_upper, Wald_p, BH_q),
                               sprintf("HR %.2f [%.2f, %.2f]\np=%.3f; exploratory", HR_per_score_SD, CI95_lower, CI95_upper, Wald_p))),
            nudge_y = 0.23, hjust = 1, family = "Arial", size = 2.65, lineheight = 0.9, colour = "black") +
  scale_colour_manual(values = c("FULL" = "#2166AC", "PARTIAL" = "#B2182B")) +
  scale_x_log10(limits = c(0.82, 2.35), breaks = c(1, 1.2, 1.5, 2)) +
  labs(title = "External survival associations", x = "HR per score SD (log scale)", y = NULL, colour = "Coverage",
       caption = NULL) + theme_pub() + theme(legend.position = "none")
save_plot(pse, "SFig2E_sparse4_external_survival_forest", 100, 80)

# Supplementary Fig. 2F: repeated nested-CV distribution replaces the retired
# fixed-coefficient random-split stability analysis. SD is descriptive, not CI.
ncv <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_REPEAT_CONCORDANCE.tsv")
summary_ncv <- jsonlite::fromJSON(file.path(ROOT, "results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_SUMMARY.json"))
write_tsv(ncv, "SFig2F_nested_repeat_concordance.tsv")
write_tsv(data.frame(
  repeat_mean_C = summary_ncv$repeat_mean_harrell_c,
  repeat_SD_descriptive = summary_ncv$repeat_sd_harrell_c,
  patient_mean_ensemble_C_distinct_estimand = summary_ncv$patient_mean_heldout_harrell_c,
  outer_folds = summary_ncv$outer_folds,
  outer_repeats = summary_ncv$outer_repeats,
  selected_n_min = summary_ncv$selected_n_min,
  selected_n_median = summary_ncv$selected_n_median,
  selected_n_max = summary_ncv$selected_n_max,
  claim_boundary = summary_ncv$claim_boundary,
  stringsAsFactors = FALSE
), "SFig2F_validation_estimands.tsv")
psf <- ggplot(ncv, aes(x = "20 repeats", y = heldout_harrell_c)) +
  geom_violin(width = 0.55, fill = "#DCEAF5", colour = "#2166AC", linewidth = 0.45) +
  geom_boxplot(width = 0.18, outlier.shape = NA, fill = "white", linewidth = 0.4) +
  geom_jitter(width = 0.07, size = 1.4, alpha = 0.75, colour = "#2166AC") +
  geom_hline(yintercept = summary_ncv$repeat_mean_harrell_c, colour = "#B2182B", linewidth = 0.55) +
  annotate("label", x = 1.26, y = max(ncv$heldout_harrell_c), hjust = 1,
           label = sprintf("Repeat mean C = %.3f\nSD = %.3f (descriptive)",
                           summary_ncv$repeat_mean_harrell_c, summary_ncv$repeat_sd_harrell_c),
           family = "Arial", size = 2.65, linewidth = 0.2) +
  scale_y_continuous(limits = c(0.50, 0.64), breaks = seq(0.50, 0.64, 0.02)) +
  labs(title = "Repeated nested cross-validation",
       x = NULL, y = "Held-out Harrell C per repeat",
       caption = "5-fold x 20 repeats; conditional on the fixed Yellow module and VST candidate set") + theme_pub()
save_plot(psf, "SFig2F_sparse4_repeated_nested_CV", 130, 85)

# Supplementary Fig. 3F: descriptive same-cohort joint median stratification
# using the current frozen score and directly measured IKZF1 expression.
ikzf1_median <- median(as.numeric(vst476["IKZF1", ]))
s3 <- data.frame(
  patient_id = clinical$patient_id,
  OS_time_months = clinical$OS_time_months,
  OS_status = clinical$OS_status,
  score_group = ifelse(clinical$score_raw_linear_predictor >= median(clinical$score_raw_linear_predictor), "Score-high", "Score-low"),
  IKZF1_group = ifelse(as.numeric(vst476["IKZF1", ]) >= ikzf1_median, "IKZF1-high", "IKZF1-low"),
  stringsAsFactors = FALSE
)
s3$joint_group <- paste(s3$score_group, s3$IKZF1_group, sep = " / ")
s3_levels <- c("Score-low / IKZF1-low", "Score-low / IKZF1-high",
               "Score-high / IKZF1-low", "Score-high / IKZF1-high")
s3$joint_group <- factor(s3$joint_group, levels = s3_levels)
stopif(nrow(s3) == 476 && !anyNA(s3$joint_group), "SFig3F joint-group mapping failed")
s3_fit <- survfit(Surv(OS_time_months, OS_status) ~ joint_group, data = s3)
s3_lr <- survdiff(Surv(OS_time_months, OS_status) ~ joint_group, data = s3)
s3_p <- pchisq(s3_lr$chisq, df = length(s3_lr$n) - 1L, lower.tail = FALSE)
s3_sum <- data.frame(
  joint_group = names(s3_lr$n), n = as.integer(s3_lr$n),
  events = as.integer(s3_lr$obs),
  median_OS_months = unname(summary(s3_fit)$table[, "median"]),
  global_logrank_chisq = s3_lr$chisq, global_logrank_df = length(s3_lr$n) - 1L,
  global_logrank_p = s3_p, score_cutoff = median(clinical$score_raw_linear_predictor),
  IKZF1_cutoff = ikzf1_median, stringsAsFactors = FALSE
)
write_tsv(s3, "SFig3F_sparse4_IKZF1_joint_KM_display_data.tsv")
write_tsv(s3_sum, "SFig3F_sparse4_IKZF1_joint_KM_summary.tsv")
s3_counts <- setNames(as.integer(table(s3$joint_group)), levels(s3$joint_group))
s3_legend <- c(sprintf("Low / low (n=%d)", s3_counts[[1]]),
               sprintf("Low / high (n=%d)", s3_counts[[2]]),
               sprintf("High / low (n=%d)", s3_counts[[3]]),
               sprintf("High / high (n=%d)", s3_counts[[4]]))
g_s3 <- ggsurvplot(
  s3_fit, data = s3, risk.table = FALSE, conf.int = FALSE,
  palette = c("#BDBDBD", "#2166AC", "#E08214", "#B2182B"),
  censor.shape = 124, censor.size = 2.0, break.time.by = 24, xlim = c(0, 120),
  legend.title = NULL, legend.labs = s3_legend,
  title = "Score x IKZF1 survival",
  xlab = "Overall survival (months)", ylab = "Survival probability",
  pval = sprintf("Global log-rank p = %.3g", s3_p), pval.coord = c(4, 0.12),
  pval.size = 2.65,
  ggtheme = theme_pub() + theme(legend.position = "top", legend.title = element_blank(),
                                legend.key.width = grid::unit(4, "mm"),
                                legend.margin = margin(0, 0, 0, 0), legend.box.spacing = grid::unit(0, "pt"))
)
g_s3$plot <- g_s3$plot +
  guides(colour = guide_legend(ncol = 2, byrow = TRUE), linetype = guide_legend(ncol = 2, byrow = TRUE))
save_surv(g_s3, "SFig3F_sparse4_IKZF1_joint_KM", 80, 80)

stats <- data.frame(
  panel = c("Fig2E", "Fig2F", "Fig2G", "Fig2H", "Fig4B", "Fig4C", "Fig4D", "Fig4E", "Fig4F", "Fig5F", "Fig6D", "Fig6E", "SFig2E", "SFig2F", "SFig3F"),
  status = "RENDERED_FROM_FROZEN_EVIDENCE",
  model_id = lock$model_id,
  stringsAsFactors = FALSE
)
write_tsv(stats, "RENDERED_PANEL_STATUS.tsv")
writeLines(capture.output(sessionInfo()), file.path(LOG, "sessionInfo.txt"))
cat(sprintf("Rendered %d active replacement panels as PDF+PNG; Fig4A duplicate and Fig4H unsupported panel retired; Fig7D user-owned\n", nrow(stats)))
