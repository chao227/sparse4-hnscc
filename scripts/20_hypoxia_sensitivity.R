#!/usr/bin/env Rscript

# ==========================================================
# Script   : 20_hypoxia_sensitivity.R
# Purpose  : Compute registered hypoxia-score sensitivity analyses
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: hypoxia sensitivity; Supplementary Table 13
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/10_sensitivities/run_hypoxia_sensitivity.R
# ==========================================================

suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

args_all <- commandArgs(trailingOnly = FALSE)
script_arg <- grep("^--file=", args_all, value = TRUE)
if (length(script_arg) != 1L) stop("Invoke this script with Rscript", call. = FALSE)
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
out_dir <- dirname(script_path)
repo <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)

need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
wtsv <- function(x, name) write.table(x, file.path(out_dir, name), sep = "\t",
                                      quote = FALSE, row.names = FALSE, na = "NA")

input_contract <- data.frame(
  role = c("analysis_contract", "full_vst", "clinical_master_476",
           "authoritative_m4_coefficients", "common_roster_crosscheck",
           "sparse4_score_crosscheck", "hallmark_json", "buffa_json",
           "prior_hypoxia_coverage_contract"),
  relative_path = c(
    "final_audit_reconciliation_20261022/ANALYSIS_CONTRACT.md",
    "inputs/tcga/expression/vst_matrix.rds",
    "agent_work/sparse4_tcga_clinical_20261016/CLINICAL_MASTER_476.tsv",
    "agent_work/sparse4_tcga_clinical_20261016/M1_M4_COMMON_ROSTER_COEFFICIENTS.tsv",
    "clinical_comparison_panel_repairs_20261019/source_data/Fig6D_adjusted_common_roster.tsv",
    "scan_1019_remediation_20261019/output/tables/STable3_TCGA_sparse4_scores.tsv",
    "protocol/gene_sets/v2026.1.Hs/HALLMARK_HYPOXIA.v2026.1.Hs.json",
    "protocol/gene_sets/v2026.1.Hs/BUFFA_HYPOXIA_METAGENE.v2026.1.Hs.json",
    "results/G_confounders/HYPOXIA_DATA_CONTRACT.tsv"
  ),
  expected_sha256 = c(
    "7e8c35c43de32be25055d202e0282cd61fc0f8398034df95ec0ba71b4cf8ae72",
    "0fdd4ea8f1697fe7714e968649b80b001a5dc89c306c47b06a7ae469542d5c84",
    "c3b30fe718beb71af4def9b84afd86d0c52b63d11ba064d7cb9bd6ae412163f4",
    "32c0427df5c35f7a76ecad6815a0b60418137e79f8cc9232c7c146f635d4e0e1",
    "b7df83ea0984345eeea1886c8dd0cc580171b22a4a12e9bbef14096233b2176b",
    "5bd81fb69dbfbb425427ed740391235d475730de0222102937e2bb075e114786",
    "ff32259a58f8ffc8266bbc783f6f40bfb18b9dc3719df00fce0e8f8cfc79e29d",
    "095e0d753b764555826b3571a69aa2797ef6800c2ab28f1363214bd2d04b8783",
    "316fcc78bcb6a5df443cbb466cf2511de16d14837909ef97eb2580690e97a067"
  ),
  stringsAsFactors = FALSE
)

# Correct the separately registered Buffa hash; the explicit assignment keeps the
# input order above readable and is checked before any result is calculated.
input_contract$expected_sha256[input_contract$role == "buffa_json"] <-
  "095e0d753b764555826b3571a69aa2797ef6800c2ab28f1363214bd2d04b8783"

input_contract$absolute_path <- vapply(input_contract$relative_path, function(p) {
  normalizePath(file.path(repo, p), mustWork = TRUE)
}, character(1))
input_contract$observed_sha256 <- vapply(input_contract$absolute_path, sha256, character(1))
input_contract$bytes <- file.info(input_contract$absolute_path)$size
input_contract$status <- ifelse(input_contract$expected_sha256 == input_contract$observed_sha256,
                                "PASS", "FAIL")
wtsv(input_contract, "INPUT_RECEIPT.tsv")
need(all(input_contract$status == "PASS"), "One or more frozen input hashes differ")

path_for <- function(role) input_contract$absolute_path[match(role, input_contract$role)]

vst <- readRDS(path_for("full_vst"))
need(is.matrix(vst) && identical(dim(vst), c(29451L, 520L)),
     "Full VST shape differs from 29451 genes x 520 samples")
need(!anyDuplicated(rownames(vst)) && !anyDuplicated(colnames(vst)),
     "Full VST identifiers are duplicated")
need(all(is.finite(vst)), "Full VST contains non-finite values")

master <- read.delim(path_for("clinical_master_476"), check.names = FALSE)
required_master <- c("patient_id", "sample_id", "OS_time_months", "OS_status",
                     "score_raw_linear_predictor", "score_per_global_sd", "age_years",
                     "stage_factor", "HPV_p16_3level", "CDE3108203_subsite",
                     "surgery_factor", "project_immune_ssGSEA", "project_stromal_ssGSEA")
need(all(required_master %in% names(master)), "Clinical master schema is incomplete")
need(nrow(master) == 476L && sum(master$OS_status) == 220L,
     "Clinical master roster/event count differs from 476/220")
need(!anyDuplicated(master$patient_id) && !anyDuplicated(master$sample_id),
     "Clinical master identifiers are duplicated")
need(identical(substr(master$sample_id, 1L, 12L), master$patient_id),
     "Clinical master sample-to-patient mapping differs")
master$stage_factor <- factor(master$stage_factor, levels = c("I_II", "III_IV"))
master$HPV_p16_3level <- factor(master$HPV_p16_3level, levels = c("Negative", "Positive", "Unknown"))
master$CDE3108203_subsite <- factor(master$CDE3108203_subsite,
                                    levels = c("OralCavity", "Oropharynx", "Larynx_Hypopharynx"))
master$surgery_factor <- factor(master$surgery_factor, levels = c("No", "Yes"))

expr_idx <- match(master$sample_id, colnames(vst))
need(!anyNA(expr_idx), "At least one locked survival sample is absent from full VST")
need(!anyDuplicated(expr_idx), "A full-VST sample was reused in the survival mapping")
vst476 <- vst[, expr_idx, drop = FALSE]
need(identical(colnames(vst476), master$sample_id), "VST order does not match clinical master")

score_ref <- read.delim(path_for("sparse4_score_crosscheck"), check.names = FALSE)
score_idx <- match(master$sample_id, score_ref$sample_id)
need(!anyNA(score_idx) && !anyDuplicated(score_ref$sample_id), "Sparse4 score crosscheck mapping failed")
need(max(abs(master$score_raw_linear_predictor - score_ref$linear_predictor[score_idx])) <= 1e-14,
     "Sparse4 raw score differs from the registered score table")
need(max(abs(master$score_per_global_sd - score_ref$score_per_global_sd[score_idx])) <= 1e-14,
     "Sparse4 standardized score differs from the registered score table")

read_signature <- function(path, object_name) {
  z <- jsonlite::fromJSON(path, simplifyVector = TRUE)
  need(object_name %in% names(z), paste("Signature object absent:", object_name))
  genes <- as.character(z[[object_name]]$geneSymbols)
  need(length(genes) > 0L && !anyNA(genes) && !anyDuplicated(genes),
       paste("Signature membership invalid:", object_name))
  genes
}
hallmark_genes <- read_signature(path_for("hallmark_json"), "HALLMARK_HYPOXIA")
buffa_genes <- read_signature(path_for("buffa_json"), "BUFFA_HYPOXIA_METAGENE")
need(length(hallmark_genes) == 200L && length(buffa_genes) == 50L,
     "Registered signature sizes differ from 200 and 50")

score_signature <- function(genes, label) {
  present <- genes[genes %in% rownames(vst476)]
  missing <- setdiff(genes, present)
  need(length(present) / length(genes) >= 0.70,
       paste("Coverage gate failed for", label))
  x <- vst476[present, , drop = FALSE]
  gene_sd <- apply(x, 1L, sd)
  need(all(is.finite(gene_sd) & gene_sd > 0), paste("Non-variable gene in", label))
  z <- sweep(sweep(x, 1L, rowMeans(x), "-"), 1L, gene_sd, "/")
  s <- colMeans(z)
  need(all(is.finite(s)), paste("Non-finite patient score in", label))
  gene_stats <- data.frame(signature = label, gene_symbol = genes,
                           covered = genes %in% present,
                           mean_476 = NA_real_, sample_sd_476 = NA_real_,
                           stringsAsFactors = FALSE)
  mm <- match(present, gene_stats$gene_symbol)
  gene_stats$mean_476[mm] <- rowMeans(x)
  gene_stats$sample_sd_476[mm] <- gene_sd
  list(score = s, present = present, missing = missing, gene_stats = gene_stats)
}

hallmark <- score_signature(hallmark_genes, "HALLMARK_HYPOXIA")
buffa <- score_signature(buffa_genes, "BUFFA_HYPOXIA_METAGENE")
coverage <- rbind(
  data.frame(signature = "HALLMARK_HYPOXIA", release = "MSigDB v2026.1.Hs",
             systematic_name = "M5891", registered_n = length(hallmark_genes),
             covered_n = length(hallmark$present), coverage_fraction = length(hallmark$present) / length(hallmark_genes),
             missing_genes = paste(hallmark$missing, collapse = ";"), coverage_gate = ">=0.70", status = "PASS"),
  data.frame(signature = "BUFFA_HYPOXIA_METAGENE", release = "MSigDB v2026.1.Hs",
             systematic_name = "M34030", registered_n = length(buffa_genes),
             covered_n = length(buffa$present), coverage_fraction = length(buffa$present) / length(buffa_genes),
             missing_genes = paste(buffa$missing, collapse = ";"), coverage_gate = ">=0.70", status = "PASS")
)
wtsv(coverage, "SIGNATURE_COVERAGE.tsv")
wtsv(rbind(hallmark$gene_stats, buffa$gene_stats), "HYPOXIA_GENE_STANDARDIZATION.tsv")

scores476 <- data.frame(
  patient_id = master$patient_id,
  sample_id = master$sample_id,
  hallmark_hypoxia = unname(hallmark$score[master$sample_id]),
  buffa_hypoxia = unname(buffa$score[master$sample_id]),
  stringsAsFactors = FALSE
)
need(all(is.finite(scores476$hallmark_hypoxia)) && all(is.finite(scores476$buffa_hypoxia)),
     "Patient hypoxia scores are incomplete")
wtsv(scores476, "PATIENT_HYPOXIA_SCORES_476.tsv")
score_summary <- rbind(
  data.frame(signature = "HALLMARK_HYPOXIA", population = "locked_476",
             n = nrow(scores476), mean = mean(scores476$hallmark_hypoxia),
             sample_sd = sd(scores476$hallmark_hypoxia), minimum = min(scores476$hallmark_hypoxia),
             maximum = max(scores476$hallmark_hypoxia)),
  data.frame(signature = "BUFFA_HYPOXIA_METAGENE", population = "locked_476",
             n = nrow(scores476), mean = mean(scores476$buffa_hypoxia),
             sample_sd = sd(scores476$buffa_hypoxia), minimum = min(scores476$buffa_hypoxia),
             maximum = max(scores476$buffa_hypoxia))
)
wtsv(score_summary, "SIGNATURE_SCORE_SUMMARY.tsv")

predictor_genes <- c("CAMK2N1", "MASP1", "MIR9-3HG", "MS4A1")
need(all(predictor_genes %in% rownames(vst476)), "One or more sparse4 genes are absent from full VST")
predictors <- cbind(sparse4_score = master$score_per_global_sd,
                    t(vst476[predictor_genes, , drop = FALSE]))
hypoxia_scores <- cbind(HALLMARK_HYPOXIA = scores476$hallmark_hypoxia,
                        BUFFA_HYPOXIA_METAGENE = scores476$buffa_hypoxia)
corr <- do.call(rbind, lapply(colnames(predictors), function(pred) {
  do.call(rbind, lapply(colnames(hypoxia_scores), function(sig) {
    z <- suppressWarnings(cor.test(predictors[, pred], hypoxia_scores[, sig],
                                   method = "spearman", exact = FALSE))
    data.frame(predictor = pred, hypoxia_signature = sig, n = nrow(master),
               spearman_rho = unname(z$estimate), p_nominal = z$p.value,
               stringsAsFactors = FALSE)
  }))
}))
need(nrow(corr) == 10L, "The frozen correlation family does not contain 10 tests")
corr$p_bh_10 <- p.adjust(corr$p_nominal, method = "BH", n = 10L)
corr$multiplicity_family <- "HYPOXIA_SPEARMAN_10"
corr$claim_boundary <- "post hoc same-cohort descriptive association; not specificity or causality"
wtsv(corr, "HYPOXIA_SPEARMAN_BH10.tsv")

m4_fields <- c("OS_time_months", "OS_status", "score_per_global_sd",
               "project_immune_ssGSEA", "age_years", "stage_factor",
               "HPV_p16_3level", "CDE3108203_subsite", "surgery_factor",
               "project_stromal_ssGSEA")
keep <- complete.cases(master[, m4_fields, drop = FALSE])
common <- droplevels(master[keep, , drop = FALSE])
need(nrow(common) == 320L && sum(common$OS_status) == 106L,
     "Reconstructed M4 roster differs from 320 patients/106 events")
common_ref <- read.delim(path_for("common_roster_crosscheck"), check.names = FALSE)
need(nrow(common_ref) == 320L && setequal(common$patient_id, common_ref$patient_id),
     "M4 patient membership differs from the registered common-roster table")
need(max(abs(common$OS_time_months[match(common_ref$patient_id, common$patient_id)] - common_ref$OS_time_months)) <= 1e-12,
     "M4 endpoint time differs from common-roster reference")
need(identical(as.integer(common$OS_status[match(common_ref$patient_id, common$patient_id)]),
               as.integer(common_ref$OS_status)), "M4 event status differs from common-roster reference")

hs_idx <- match(common$patient_id, scores476$patient_id)
need(!anyNA(hs_idx), "Hypoxia score-to-M4 mapping failed")
common$hallmark_hypoxia_raw <- scores476$hallmark_hypoxia[hs_idx]
common$buffa_hypoxia_raw <- scores476$buffa_hypoxia[hs_idx]
common$hallmark_hypoxia_per_m4_sd <- as.numeric(scale(common$hallmark_hypoxia_raw))
common$buffa_hypoxia_per_m4_sd <- as.numeric(scale(common$buffa_hypoxia_raw))
need(abs(mean(common$hallmark_hypoxia_per_m4_sd)) < 1e-12 && abs(sd(common$hallmark_hypoxia_per_m4_sd) - 1) < 1e-12,
     "Hallmark M4-roster standardization failed")
need(abs(mean(common$buffa_hypoxia_per_m4_sd)) < 1e-12 && abs(sd(common$buffa_hypoxia_per_m4_sd) - 1) < 1e-12,
     "Buffa M4-roster standardization failed")
wtsv(common, "M4_HYPOXIA_ROSTER_320.tsv")

f_base <- Surv(OS_time_months, OS_status) ~ score_per_global_sd +
  project_immune_ssGSEA + age_years + stage_factor + HPV_p16_3level +
  CDE3108203_subsite + surgery_factor + project_stromal_ssGSEA
fit_cox <- function(formula, data) {
  z <- coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  need(is.null(z$fail) && all(is.finite(coef(z))), "Cox fit failed or has non-finite coefficients")
  z
}
base_fit <- fit_cox(f_base, common)
fit_h <- fit_cox(update(f_base, . ~ . + hallmark_hypoxia_per_m4_sd), common)
fit_b <- fit_cox(update(f_base, . ~ . + buffa_hypoxia_per_m4_sd), common)

reference_coef <- read.delim(path_for("authoritative_m4_coefficients"), check.names = FALSE)
reference_coef <- reference_coef[reference_coef$model_id == "M4", , drop = FALSE]
need(length(coef(base_fit)) == nrow(reference_coef), "M4 coefficient count differs from authority")
base_cmp <- data.frame(
  term = names(coef(base_fit)),
  reproduced_coefficient = unname(coef(base_fit)),
  authoritative_coefficient = reference_coef$coefficient[match(names(coef(base_fit)), reference_coef$term)]
)
need(!anyNA(base_cmp$authoritative_coefficient), "M4 coefficient term mapping failed")
base_cmp$absolute_difference <- abs(base_cmp$reproduced_coefficient - base_cmp$authoritative_coefficient)
base_cmp$status <- ifelse(base_cmp$absolute_difference <= 1e-12, "PASS", "FAIL")
wtsv(base_cmp, "M4_BASELINE_REPRODUCTION.tsv")
need(all(base_cmp$status == "PASS"), "Baseline M4 reproduction differs from authority")

coef_rows <- function(fit, model_id, signature_term) {
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients)
  ci <- as.data.frame(sm$conf.int)
  data.frame(
    model_id = model_id, signature_term = signature_term,
    term = rownames(cf), n = fit$n, events = fit$nevent,
    coefficient = cf$coef, se = cf$`se(coef)`,
    hazard_ratio = ci$`exp(coef)`, ci_lower_95 = ci$`lower .95`,
    ci_upper_95 = ci$`upper .95`, p_nominal = cf$`Pr(>|z|)`,
    stringsAsFactors = FALSE, row.names = NULL
  )
}
coefs <- rbind(coef_rows(fit_h, "M4_PLUS_HALLMARK", "hallmark_hypoxia_per_m4_sd"),
               coef_rows(fit_b, "M4_PLUS_BUFFA", "buffa_hypoxia_per_m4_sd"))
coefs$claim_boundary <- "post hoc same-cohort confounding sensitivity; Wald inference is descriptive"
wtsv(coefs, "M4_HYPOXIA_COEFFICIENTS.tsv")

lrt_row <- function(fit, id, term) {
  aa <- anova(base_fit, fit, test = "LRT")
  data.frame(model_id = id, added_term = term, n = fit$n, events = fit$nevent,
             base_loglik = as.numeric(logLik(base_fit)), augmented_loglik = as.numeric(logLik(fit)),
             lr_chisq = aa[2L, "Chisq"], added_df = aa[2L, "Df"], p_nominal = aa[2L, "Pr(>|Chi|)"],
             delta_AIC_base_minus_augmented = AIC(base_fit) - AIC(fit),
             multiplicity = "none; two registered post hoc sensitivities reported separately",
             claim_boundary = "descriptive nested-model comparison; not confirmatory", stringsAsFactors = FALSE)
}
lr <- rbind(lrt_row(fit_h, "M4_PLUS_HALLMARK", "hallmark_hypoxia_per_m4_sd"),
            lrt_row(fit_b, "M4_PLUS_BUFFA", "buffa_hypoxia_per_m4_sd"))
wtsv(lr, "M4_HYPOXIA_LR.tsv")

ph_rows <- function(fit, id) {
  z <- cox.zph(fit, transform = "km", terms = TRUE, singledf = TRUE)
  tab <- as.data.frame(z$table)
  tab$term <- rownames(tab)
  rownames(tab) <- NULL
  term_idx <- tab$term != "GLOBAL"
  tab$p_bh_within_model_terms <- NA_real_
  tab$p_bh_within_model_terms[term_idx] <- p.adjust(tab$p[term_idx], method = "BH")
  data.frame(model_id = id, term = tab$term, chisq = tab$chisq, df = tab$df,
             p_nominal = tab$p, p_bh_within_model_terms = tab$p_bh_within_model_terms,
             transform = "km", stringsAsFactors = FALSE)
}
ph <- rbind(ph_rows(fit_h, "M4_PLUS_HALLMARK"), ph_rows(fit_b, "M4_PLUS_BUFFA"))
wtsv(ph, "M4_HYPOXIA_PH_DIAGNOSTICS.tsv")

qa <- data.frame(
  check_id = c("input_hashes", "vst_shape", "vst_unique_ids", "locked_roster",
               "sample_mapping", "score_identity", "hallmark_identity_coverage",
               "buffa_identity_coverage", "correlation_family", "m4_roster",
               "m4_event_count", "m4_membership", "baseline_m4_reproduction",
               "hallmark_and_buffa_separate", "no_bootstrap"),
  observed = c("9/9 PASS", paste(dim(vst), collapse = "x"), "TRUE", "476 patients;220 events",
               "476/476 exact sample IDs", "raw and standardized sparse4 exact",
               sprintf("M5891;%d/%d", length(hallmark$present), length(hallmark_genes)),
               sprintf("M34030;%d/%d", length(buffa$present), length(buffa_genes)),
               sprintf("%d tests", nrow(corr)), sprintf("%d patients", nrow(common)),
               sprintf("%d events", sum(common$OS_status)), "exact set equality",
               sprintf("max abs diff %.3g", max(base_cmp$absolute_difference)),
               "two distinct augmented fits", "TRUE"),
  acceptance_rule = c("all frozen hashes exact", "29451x520", "no duplicate gene/sample IDs",
                      "476 patients;220 events", "all survival aliquots map exactly once",
                      "difference <=1e-14", ">=70% of registered 200 genes",
                      ">=70% of registered 50 genes", "exactly 10; BH across all 10",
                      "320 patients", "106 events", "matches registered common roster",
                      "all coefficient differences <=1e-12", "never jointly entered", "no bootstrap call"),
  status = "PASS", stringsAsFactors = FALSE
)
wtsv(qa, "QA_CHECKS.tsv")

key_coef <- coefs[coefs$term %in% c("score_per_global_sd", "hallmark_hypoxia_per_m4_sd",
                                    "buffa_hypoxia_per_m4_sd"), , drop = FALSE]
hallmark_score <- key_coef[key_coef$model_id == "M4_PLUS_HALLMARK" &
                             key_coef$term == "score_per_global_sd", , drop = FALSE]
hallmark_hyp <- key_coef[key_coef$model_id == "M4_PLUS_HALLMARK" &
                           key_coef$term == "hallmark_hypoxia_per_m4_sd", , drop = FALSE]
buffa_score <- key_coef[key_coef$model_id == "M4_PLUS_BUFFA" &
                          key_coef$term == "score_per_global_sd", , drop = FALSE]
buffa_hyp <- key_coef[key_coef$model_id == "M4_PLUS_BUFFA" &
                        key_coef$term == "buffa_hypoxia_per_m4_sd", , drop = FALSE]

report <- c(
  "# Producer report: post hoc hypoxia sensitivity",
  "",
  "Status: `COMPUTED_PENDING_INDEPENDENT_AUDIT`.",
  "",
  "This production run implements frozen contract C only. It uses the locked 476-patient TCGA survival roster and full VST matrix, the registered MSigDB v2026.1.Hs HALLMARK_HYPOXIA (M5891) and BUFFA_HYPOXIA_METAGENE (M34030) memberships, and the exact 320-patient/106-event M4 common roster.",
  "",
  sprintf("Coverage was %d/%d (%.1f%%) for Hallmark and %d/%d (%.1f%%) for Buffa.",
          length(hallmark$present), length(hallmark_genes), 100 * length(hallmark$present) / length(hallmark_genes),
          length(buffa$present), length(buffa_genes), 100 * length(buffa$present) / length(buffa_genes)),
  "Each covered gene was standardized across the locked 476 patients and signatures were equal-weight means. For Cox sensitivity, each final signature score was standardized once within the exact M4 roster, so its HR is per one M4-roster SD.",
  "",
  "The descriptive correlation family contains exactly 10 tests (sparse4 plus its four genes against two signatures), with BH adjustment across all 10. Hallmark and Buffa were fitted in separate M4-augmented Cox models. The output includes all augmented-model coefficients, nested LR comparisons, and complete term/global proportional-hazards diagnostics.",
  "",
  sprintf("Sparse4 correlated with Hallmark at rho=%.3f (BH q=%.3g) and with Buffa at rho=%.3f (BH q=%.3g). These are same-cohort descriptive associations.",
          corr$spearman_rho[corr$predictor == "sparse4_score" & corr$hypoxia_signature == "HALLMARK_HYPOXIA"],
          corr$p_bh_10[corr$predictor == "sparse4_score" & corr$hypoxia_signature == "HALLMARK_HYPOXIA"],
          corr$spearman_rho[corr$predictor == "sparse4_score" & corr$hypoxia_signature == "BUFFA_HYPOXIA_METAGENE"],
          corr$p_bh_10[corr$predictor == "sparse4_score" & corr$hypoxia_signature == "BUFFA_HYPOXIA_METAGENE"]),
  sprintf("After separate Hallmark adjustment, sparse4 HR=%.3f (95%% CI %.3f-%.3f); Hallmark HR=%.3f (%.3f-%.3f), LR P=%.3g.",
          hallmark_score$hazard_ratio, hallmark_score$ci_lower_95, hallmark_score$ci_upper_95,
          hallmark_hyp$hazard_ratio, hallmark_hyp$ci_lower_95, hallmark_hyp$ci_upper_95,
          lr$p_nominal[lr$model_id == "M4_PLUS_HALLMARK"]),
  sprintf("After separate Buffa adjustment, sparse4 HR=%.3f (95%% CI %.3f-%.3f); Buffa HR=%.3f (%.3f-%.3f), LR P=%.3g.",
          buffa_score$hazard_ratio, buffa_score$ci_lower_95, buffa_score$ci_upper_95,
          buffa_hyp$hazard_ratio, buffa_hyp$ci_lower_95, buffa_hyp$ci_upper_95,
          lr$p_nominal[lr$model_id == "M4_PLUS_BUFFA"]),
  sprintf("Global proportional-hazards P values were %.3f (Hallmark-augmented) and %.3f (Buffa-augmented).",
          ph$p_nominal[ph$model_id == "M4_PLUS_HALLMARK" & ph$term == "GLOBAL"],
          ph$p_nominal[ph$model_id == "M4_PLUS_BUFFA" & ph$term == "GLOBAL"]),
  "",
  "No bootstrap, refitting of sparse4, gene-set switching, joint Hallmark/Buffa model, imputation, or causal/specificity claim was performed. Results are post hoc same-cohort confounding sensitivities and require a different-role independent audit before use."
)
writeLines(report, file.path(out_dir, "PRODUCER_REPORT.md"))

self_review <- c(
  "# Producer Level-3 self-review",
  "",
  "- Data contract: all nine frozen inputs were hash-checked before analysis; full VST shape and unique identifiers passed.",
  "- Mapping: all 476 survival aliquots mapped exactly once; the 320/106 M4 roster matched the separate registered common-roster table.",
  "- Estimand and unit: correlations use 476 patients; Cox sensitivities use the exact 320 patients, with patient as the independent unit.",
  "- Statistics: exactly 10 registered Spearman tests use one BH family; Cox LR/Wald results are descriptive post hoc sensitivities; complete term/global PH diagnostics are exported.",
  "- Reproducibility: registered memberships, covered/missing genes, per-gene 476-patient means/SDs, patient scores, code, input hashes and session information are retained.",
  "- Claim ceiling: association with hypoxia signatures and persistence of the sparse4 coefficient after adjustment do not prove B-cell specificity, absence of confounding, causality or external validity.",
  "- Independent validation: not performed by this producer; status remains COMPUTED_PENDING_INDEPENDENT_AUDIT."
)
writeLines(self_review, file.path(out_dir, "SELF_REVIEW.md"))

sink(file.path(out_dir, "SESSIONINFO.txt"))
cat("Run time: ", format(Sys.time(), tz = "Asia/Shanghai", usetz = TRUE), "\n", sep = "")
cat("R executable: ", normalizePath(file.path(R.home("bin"), "R"), mustWork = TRUE), "\n", sep = "")
sessionInfo()
sink()

output_files <- setdiff(list.files(out_dir, full.names = TRUE),
                        c(file.path(out_dir, "SHA256SUMS.txt"), file.path(out_dir, "run.log")))
output_files <- output_files[file.info(output_files)$isdir %in% FALSE]
output_hashes <- data.frame(
  sha256 = vapply(output_files, sha256, character(1)),
  file = basename(output_files), stringsAsFactors = FALSE
)
write.table(output_hashes, file.path(out_dir, "SHA256SUMS.txt"), sep = "  ",
            quote = FALSE, row.names = FALSE, col.names = FALSE)

cat("COMPUTED_PENDING_INDEPENDENT_AUDIT\n")
