#!/usr/bin/env Rscript

# ==========================================================
# Script   : clinical_comparison_panel_repairs_20261019_recompute_and_render.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:clinical_comparison_panel_repairs_20261019/scripts/recompute_and_render.R
# ==========================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(ragg)
})

ROOT <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
MOD <- file.path(ROOT, "clinical_comparison_panel_repairs_20261019")
SRC <- file.path(MOD, "source_data")
FIG <- file.path(MOD, "figures")
REP <- file.path(MOD, "reports")
dir.create(SRC, recursive = TRUE, showWarnings = FALSE)
dir.create(FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(REP, recursive = TRUE, showWarnings = FALSE)

stopif <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
read_tsv <- function(path) read.delim(file.path(ROOT, path), check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, name) write.table(x, file.path(SRC, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
znum <- function(x) as.numeric((x - mean(x)) / sd(x))

clinical <- read_tsv("agent_work/sparse4_tcga_clinical_20261016/CLINICAL_MASTER_476.tsv")
score <- read_tsv("results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/SCORE_TABLE.tsv")
vst <- readRDS(file.path(ROOT, "inputs/tcga/expression/vst_matrix.rds"))
wgcna <- readRDS(file.path(ROOT, "inputs/wgcna/legacy/wgcna_results.rds"))
traits <- readRDS(file.path(ROOT, "inputs/wgcna/legacy/trait_matrix.rds"))
vst25 <- readRDS(file.path(ROOT, "inputs/wgcna/legacy/vst_top25mad.rds"))

stopif(nrow(clinical) == 476L && sum(clinical$OS_status) == 220L, "Frozen clinical roster drift")
stopif(identical(clinical$sample_id, score$sample_id), "Clinical/score order mismatch")
stopif(max(abs(clinical$score_raw_linear_predictor - score$linear_predictor)) < 1e-12, "Sparse4 authority mismatch")

vi <- match(clinical$sample_id, colnames(vst))
wi <- match(clinical$sample_id, wgcna$sample_order)
ti <- match(clinical$sample_id, traits$sample_id)
stopif(!anyNA(vi) && !anyDuplicated(vi), "VST mapping failed")
stopif(!anyNA(ti) && !anyDuplicated(ti), "Trait mapping failed")
stopif(sum(is.na(wi)) == 11L, "Historical WGCNA overlap changed")

v476 <- vst[, vi, drop = FALSE]
cd79a <- as.numeric(v476["CD79A", ])
ikzf1 <- as.numeric(v476["IKZF1", ])
ms4a1 <- as.numeric(v476["MS4A1", ])
cd_hi <- cd79a >= median(cd79a)
ik_hi <- ikzf1 >= median(ikzf1)
lymphhot <- cd_hi & ik_hi
stopif(identical(as.integer(table(cd_hi)), c(238L, 238L)), "CD79A median split drift")
stopif(sum(lymphhot) == 171L, "LymphHot count drift")

tls_genes <- c("CCL2", "CCL3", "CCL4", "CCL5", "CCL8", "CCL18", "CCL19", "CCL21", "CXCL9", "CXCL10", "CXCL11", "CXCL13")
stopif(all(tls_genes %in% rownames(v476)), "TLS12 gene coverage failed")
tls12 <- colMeans(t(apply(v476[tls_genes, , drop = FALSE], 1, znum)))

# Recover the historical Yellow eigengene on its original 507 training samples,
# then project the 11 current-roster samples excluded as WGCNA outliers using
# fixed training centering/scaling and PC1 loadings. This is outcome-blind and
# reproduces every stored training-sample ME to numerical precision.
yellow_genes <- names(wgcna$module_colors)[wgcna$module_colors == "yellow"]
stopif(length(yellow_genes) == 539L && all(yellow_genes %in% rownames(vst25)), "Yellow membership/coverage drift")
yellow_train <- t(vst25[yellow_genes, wgcna$sample_order, drop = FALSE])
yellow_pc <- prcomp(yellow_train, center = TRUE, scale. = TRUE, rank. = 1)
yellow_calibration <- lm(wgcna$MEs$MEyellow ~ yellow_pc$x[, 1])
stopif(abs(cor(yellow_pc$x[, 1], wgcna$MEs$MEyellow)) > 0.999999999999,
       "PC1 does not reproduce historical Yellow ME")
yellow_current_pc <- predict(yellow_pc, newdata = t(vst25[yellow_genes, clinical$sample_id, drop = FALSE]))[, 1]
yellow_projected <- unname(coef(yellow_calibration)[1] + coef(yellow_calibration)[2] * yellow_current_pc)
yellow_saved_current <- as.numeric(wgcna$MEs$MEyellow[wi])
yellow_overlap_diff <- abs(yellow_projected[!is.na(wi)] - yellow_saved_current[!is.na(wi)])
stopif(max(yellow_overlap_diff) < 1e-12, "Fixed-loading Yellow projection does not reproduce stored ME")

master <- data.frame(
  patient_id = clinical$patient_id,
  sample_id = clinical$sample_id,
  OS_time_months = clinical$OS_time_months,
  OS_status = clinical$OS_status,
  sparse4 = clinical$score_raw_linear_predictor,
  MS4A1 = ms4a1,
  CD79A = cd79a,
  IKZF1 = ikzf1,
  CD79A_high = cd_hi,
  LymphHot = lymphhot,
  ImmuneScore = as.numeric(traits$Immune_Score[ti]),
  TLS12 = tls12,
  Yellow_ME = yellow_projected,
  stringsAsFactors = FALSE
)

# Risk-direction contract was frozen before outcome analysis. Positive values
# below always denote greater risk, so concordance uses reverse=TRUE.
master$risk_sparse4 <- master$sparse4
master$risk_MS4A1 <- -master$MS4A1
master$risk_CD79A <- -master$CD79A
master$risk_Yellow_ME <- -master$Yellow_ME
master$risk_LymphHot <- -as.numeric(master$LymphHot)
master$risk_ImmuneScore <- -master$ImmuneScore
master$risk_TLS12 <- -master$TLS12
master$risk_CD79A_high <- -as.numeric(master$CD79A_high)
master$risk_IKZF1 <- -master$IKZF1

stopif(cor(master$Yellow_ME, master$CD79A, use = "complete.obs") > 0.8,
       "Yellow eigengene no longer tracks the B-cell axis")

concordance_family <- function(d, predictors, labels, family_id) {
  stopif(all(complete.cases(d[, c("OS_time_months", "OS_status", predictors)])),
         paste(family_id, "has missing predictor/outcome values"))
  f <- as.formula(paste("Surv(OS_time_months, OS_status) ~", paste(predictors, collapse = " + ")))
  cc <- concordance(f, data = d, reverse = TRUE, timewt = "n")
  est <- as.numeric(coef(cc))
  vv <- vcov(cc)
  stopif(length(est) == length(predictors) && all(dim(vv) == length(predictors)),
         paste(family_id, "concordance covariance dimension failure"))
  se <- sqrt(diag(vv))
  absolute <- data.frame(
    family_id = family_id, predictor_id = predictors, label = labels,
    n = nrow(d), events = sum(d$OS_status), harrell_c = est,
    se_ij = se, ci_lower = pmax(0, est - qnorm(.975) * se),
    ci_upper = pmin(1, est + qnorm(.975) * se),
    variance_method = "survival::concordance infinitesimal jackknife; no bootstrap",
    stringsAsFactors = FALSE
  )
  list(absolute = absolute, covariance = vv)
}

paired_against <- function(abs, vv, reference_id, comparison_ids, family_name) {
  ri <- match(reference_id, abs$predictor_id)
  do.call(rbind, lapply(comparison_ids, function(cid) {
    ci <- match(cid, abs$predictor_id)
    contrast <- rep(0, nrow(abs)); contrast[ri] <- 1; contrast[ci] <- -1
    delta <- sum(contrast * abs$harrell_c)
    se <- sqrt(drop(t(contrast) %*% vv %*% contrast))
    p <- 2 * pnorm(-abs(delta / se))
    data.frame(
      multiplicity_family = family_name, reference_id = reference_id, comparator_id = cid,
      n = abs$n[1], events = abs$events[1], c_reference = abs$harrell_c[ri],
      c_comparator = abs$harrell_c[ci], delta_c = delta, se_ij = se,
      ci_lower = delta - qnorm(.975) * se, ci_upper = delta + qnorm(.975) * se,
      p_nominal = p, stringsAsFactors = FALSE
    )
  }))
}

# -----------------------------------------------------------------------------
# Fig.4E: preserve the submitted seven-predictor task on a true common roster.
# Eleven current patients lacked stored MEs because they were excluded as WGCNA
# outliers; the outcome-blind fixed-loading projection restores the full roster.
# -----------------------------------------------------------------------------
f4_vars <- c("risk_sparse4", "risk_MS4A1", "risk_CD79A", "risk_Yellow_ME",
             "risk_LymphHot", "risk_ImmuneScore", "risk_TLS12")
f4_labels <- c("B-cell score", "MS4A1", "CD79A", "Yellow ME", "LymphHot", "ImmuneScore", "TLS12")
f4_ok <- complete.cases(master[, c("OS_time_months", "OS_status", f4_vars)])
f4 <- master[f4_ok, , drop = FALSE]
stopif(nrow(f4) == 476L && sum(f4$OS_status) == 220L, "Fig4E common roster drift")
f4_fit <- concordance_family(f4, f4_vars, f4_labels, "FIG4E_SEVEN_PREDICTOR_COMMON_ROSTER")
f4_abs <- f4_fit$absolute
f4_pair <- paired_against(f4_abs, f4_fit$covariance, "risk_sparse4", setdiff(f4_vars, "risk_sparse4"), "FIG4E_SPARSE4_VS_SIX_HOLM6")
f4_pair$p_holm6 <- p.adjust(f4_pair$p_nominal, method = "holm")
f4_pair$call <- ifelse(f4_pair$p_holm6 < 0.05, "DIFFERENT_AFTER_HOLM6", "NO_DIFFERENCE_AFTER_HOLM6")
write_tsv(f4_abs, "Fig4E_seven_predictor_apparent_C.tsv")
write_tsv(f4_pair, "Fig4E_sparse4_paired_contrasts.tsv")
write_tsv(master[, c("patient_id", "sample_id", "OS_time_months", "OS_status", f4_vars)], "Fig4E_patient_predictor_ledger.tsv")
write_tsv(data.frame(
  patient_id = master$patient_id, sample_id = master$sample_id,
  stored_historical_ME = !is.na(wi), projected_with_fixed_training_loadings = is.na(wi),
  Yellow_ME = master$Yellow_ME,
  overlap_absolute_difference = ifelse(is.na(wi), NA_real_, yellow_overlap_diff[match(seq_len(nrow(master)), which(!is.na(wi)))]),
  stringsAsFactors = FALSE
), "Fig4E_Yellow_ME_projection_audit.tsv")

# -----------------------------------------------------------------------------
# Fig.6D: preserve the submitted four-object discrimination task on n=476/e=220.
# -----------------------------------------------------------------------------
f6_vars <- c("risk_sparse4", "risk_LymphHot", "risk_CD79A_high", "risk_CD79A")
f6_labels <- c("B-cell score", "LymphHot", "CD79A-high", "CD79A")
f6_fit <- concordance_family(master, f6_vars, f6_labels, "FIG6D_ORIGINAL_FOUR_OBJECT_COMMON_ROSTER")
f6_abs <- f6_fit$absolute
f6_pair <- paired_against(f6_abs, f6_fit$covariance, "risk_LymphHot", setdiff(f6_vars, "risk_LymphHot"), "FIG6D_LYMPHHOT_VS_THREE_HOLM3")
f6_pair$p_holm3 <- p.adjust(f6_pair$p_nominal, method = "holm")
f6_pair$call <- ifelse(f6_pair$p_holm3 < 0.05, "DIFFERENT_AFTER_HOLM3", "NO_DIFFERENCE_AFTER_HOLM3")
write_tsv(f6_abs, "Fig6D_four_predictor_apparent_C.tsv")
write_tsv(f6_pair, "Fig6D_LymphHot_paired_contrasts.tsv")
write_tsv(master[, c("patient_id", "sample_id", "OS_time_months", "OS_status", "sparse4", "CD79A", "IKZF1", "CD79A_high", "LymphHot")],
          "Fig6D_patient_predictor_ledger.tsv")

# Supporting incremental analysis: same frozen M4 complete-case roster and
# established extended clinical base. It does not replace the primary panel.
aug <- cbind(clinical, master[, c("CD79A", "IKZF1", "CD79A_high", "LymphHot")])
aug$z_CD79A <- znum(aug$CD79A)
aug$z_IKZF1 <- znum(aug$IKZF1)
need_cols <- c("OS_time_months", "OS_status", "score_per_global_sd", "age_years", "stage_factor",
               "HPV_p16_3level", "CDE3108203_subsite", "surgery_factor",
               "project_immune_ssGSEA", "project_stromal_ssGSEA", "z_CD79A", "z_IKZF1",
               "CD79A_high", "LymphHot")
common <- aug[complete.cases(aug[, need_cols]), , drop = FALSE]
stopif(nrow(common) == 320L && sum(common$OS_status) == 106L, "Adjusted common roster drift")
base_formula <- Surv(OS_time_months, OS_status) ~ age_years + stage_factor + HPV_p16_3level +
  CDE3108203_subsite + surgery_factor + project_immune_ssGSEA + project_stromal_ssGSEA
aug_terms <- c(sparse4 = "score_per_global_sd", LymphHot = "LymphHot", CD79A_high = "CD79A_high",
               CD79A_continuous = "z_CD79A", IKZF1_continuous = "z_IKZF1")
base_fit <- coxph(base_formula, data = common, x = TRUE, y = TRUE, model = TRUE, ties = "efron")
aug_fits <- lapply(aug_terms, function(term) {
  coxph(update(base_formula, paste(". ~ . +", term)), data = common, x = TRUE, y = TRUE, model = TRUE, ties = "efron")
})
fits_all <- c(list(clinical_base = base_fit), aug_fits)
lp <- as.data.frame(lapply(fits_all, predict, type = "lp"), check.names = FALSE)
lp$OS_time_months <- common$OS_time_months; lp$OS_status <- common$OS_status
lp_vars <- names(fits_all)
lp_fit <- concordance_family(lp, lp_vars, lp_vars, "FIG6D_ADJUSTED_INCREMENTAL_COMMON_ROSTER")
lp_abs <- lp_fit$absolute
lp_delta_raw <- paired_against(lp_abs, lp_fit$covariance, "clinical_base", names(aug_terms), "FIG6D_FIVE_AUGMENTED_VS_BASE_HOLM5")
lp_delta <- data.frame(
  multiplicity_family = lp_delta_raw$multiplicity_family,
  predictor_id = lp_delta_raw$comparator_id,
  n = lp_delta_raw$n, events = lp_delta_raw$events,
  c_augmented = lp_delta_raw$c_comparator, c_base = lp_delta_raw$c_reference,
  delta_c_augmented_minus_base = -lp_delta_raw$delta_c,
  se_ij = lp_delta_raw$se_ij,
  ci_lower = -lp_delta_raw$ci_upper, ci_upper = -lp_delta_raw$ci_lower,
  p_nominal = lp_delta_raw$p_nominal,
  stringsAsFactors = FALSE
)
lp_delta$p_holm5 <- p.adjust(lp_delta$p_nominal, method = "holm")
lrt <- do.call(rbind, lapply(names(aug_fits), function(id) {
  a <- anova(base_fit, aug_fits[[id]], test = "LRT")
  data.frame(predictor_id = id, n = nrow(common), events = sum(common$OS_status),
             c_augmented = lp_abs$harrell_c[match(id, lp_abs$predictor_id)],
             c_base = lp_abs$harrell_c[match("clinical_base", lp_abs$predictor_id)],
             delta_lr_chisq = a[2, "Chisq"], delta_df = a[2, "Df"], p_lr = a[2, "Pr(>|Chi|)"],
             stringsAsFactors = FALSE)
}))
lrt$p_lr_holm5 <- p.adjust(lrt$p_lr, method = "holm")
inc <- merge(lp_delta, lrt, by = "predictor_id", sort = FALSE, suffixes = c("_pairedC", "_LRT"))
inc <- inc[match(names(aug_terms), inc$predictor_id), ]
inc$incremental_call <- ifelse(inc$p_holm5 < 0.05 & inc$p_lr_holm5 < 0.05,
                               "SUPPORTED_BY_BOTH_C_AND_LR_AFTER_HOLM5", "NOT_SUPPORTED_BY_BOTH_AFTER_HOLM5")
write_tsv(lp_abs, "Fig6D_adjusted_model_apparent_C.tsv")
write_tsv(inc, "Fig6D_adjusted_incremental_support.tsv")
write_tsv(common[, c("patient_id", "sample_id", need_cols)], "Fig6D_adjusted_common_roster.tsv")

ph <- do.call(rbind, lapply(names(aug_fits), function(id) {
  z <- cox.zph(aug_fits[[id]], transform = "km", terms = TRUE, singledf = TRUE)
  data.frame(model_id = id, term = rownames(z$table), chisq = z$table[, "chisq"],
             df = z$table[, "df"], p = z$table[, "p"], stringsAsFactors = FALSE)
}))
write_tsv(ph, "Fig6D_adjusted_PH_diagnostics.tsv")

# -----------------------------------------------------------------------------
# Rendering: original panel sizes, native visible text >=7.5 pt.
# -----------------------------------------------------------------------------
theme_pub <- function() {
  theme_classic(base_family = "Arial", base_size = 7.5) +
    theme(axis.text = element_text(size = 7.5, colour = "black"),
          axis.title = element_text(size = 7.5, colour = "black"),
          plot.subtitle = element_text(size = 7.5, hjust = 0.5, lineheight = 0.95),
          plot.caption = element_text(size = 7.5, hjust = 0, lineheight = 0.92, colour = "grey20"),
          plot.margin = margin(2.5, 2.5, 2.5, 2.5))
}
save_panel <- function(p, stem, w_mm, h_mm) {
  cairo_pdf(file.path(FIG, paste0(stem, ".pdf")), width = w_mm / 25.4, height = h_mm / 25.4, family = "Arial")
  print(p); dev.off()
  agg_png(file.path(FIG, paste0(stem, ".png")), width = round(w_mm / 25.4 * 600),
          height = round(h_mm / 25.4 * 600), units = "px", res = 600, bg = "white")
  print(p); dev.off()
}

f4_abs$label <- factor(f4_abs$label, levels = rev(f4_labels))
x4 <- range(c(f4_abs$ci_lower, f4_abs$ci_upper))
p4 <- ggplot(f4_abs, aes(harrell_c, label, xmin = ci_lower, xmax = ci_upper,
                         colour = predictor_id == "risk_sparse4")) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_errorbar(orientation = "y", width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#C0614E"), guide = "none") +
  scale_x_continuous(limits = c(min(0.48, x4[1] - 0.01), x4[2] + 0.01)) +
  labs(x = "Apparent Harrell C (95% IJ CI)", y = NULL) +
  theme_pub()
save_panel(p4, "Fig4E_seven_predictor_apparent_C", 90, 60)

f6_abs$label <- factor(f6_abs$label, levels = rev(f6_labels))
p6 <- ggplot(f6_abs, aes(harrell_c, label, xmin = ci_lower, xmax = ci_upper,
                         colour = predictor_id == "risk_LymphHot")) +
  geom_vline(xintercept = 0.5, linetype = "dashed", colour = "grey55", linewidth = 0.3) +
  geom_errorbar(orientation = "y", width = 0.16, linewidth = 0.55) +
  geom_point(size = 2.2) +
  scale_colour_manual(values = c(`FALSE` = "#4E79A7", `TRUE` = "#1FA187"), guide = "none") +
  scale_x_continuous(limits = c(0.47, 0.695), breaks = c(0.50, 0.55, 0.60, 0.65)) +
  labs(x = "Apparent C (95% IJ CI)", y = NULL) +
  theme_pub() + theme(axis.text.y = element_text(size = 7.5))
save_panel(p6, "Fig6D_LymphHot_simple_predictor_C", 60, 60)

qa <- data.frame(
  check_id = c("frozen_roster", "frozen_score", "vst_mapping", "trait_mapping", "yellow_overlap",
               "lymphhot_definition", "yellow_projection", "fig4_common_roster", "fig6_common_roster", "fig4_family",
               "fig6_family", "adjusted_roster", "no_bootstrap", "native_text_floor"),
  status = "PASS",
  observed = c("476/220", "exact", "476/476", "476/476", "465/476",
               "171 LymphHot", "max overlap error <1e-12", "476/220", "476/220", "6 Holm contrasts",
               "3 Holm contrasts", "320/106", "TRUE", "7.5 pt"),
  expected = c("476/220", "exact", "476/476", "476/476", "465/476",
               "171 LymphHot", "max overlap error <1e-12", "476/220", "476/220", "6 Holm contrasts",
               "3 Holm contrasts", "320/106", "TRUE", ">=7.5 pt"),
  stringsAsFactors = FALSE
)
write.table(qa, file.path(REP, "NUMERICAL_QA.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(REP, "SESSIONINFO.txt"))

cat("Clinical comparison panels complete\n")
print(f4_abs)
print(f4_pair)
print(f6_abs)
print(f6_pair)
print(inc)
