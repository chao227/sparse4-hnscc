#!/usr/bin/env Rscript

# ==========================================================
# Script   : 13_survival_models_m1_m4.R
# Purpose  : Fit common-roster M1-M4 Cox models and diagnostics
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: survival models; Fig. 4E; Main Table 2
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/07_clinical/run_v3_1_p1_06_authoritative_cox.R
# ==========================================================

suppressPackageStartupMessages({
  library(survival)
  library(jsonlite)
})

options(stringsAsFactors = FALSE)
all_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", all_args, value = TRUE)
if (length(file_arg) != 1L) stop("Runner must be invoked by Rscript from its regular file", call. = FALSE)
runner_path <- normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE)
out <- dirname(runner_path)
repo <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)

need <- function(ok, msg, ...) if (!isTRUE(ok)) stop(sprintf(msg, ...), call. = FALSE)
sha <- function(path) {
  need(file.exists(path) && !isTRUE(file.info(path)$isdir), "Missing file for hash: %s", path)
  digest::digest(file = path, algo = "sha256", serialize = FALSE)
}
write_tsv <- function(x, name) write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write_json <- function(x, name) jsonlite::write_json(x, file.path(out, name), auto_unbox = TRUE, pretty = TRUE, digits = NA, na = "null")
runtime_contract <- read.delim(file.path(out, "RUNTIME_INPUT_CONTRACT.tsv"), check.names = FALSE)
need(identical(names(runtime_contract), c("input_role", "relative_path", "expected_sha256")), "Runtime input contract schema changed")
need(!anyDuplicated(runtime_contract$input_role) && !anyDuplicated(runtime_contract$relative_path), "Runtime input contract contains duplicate roles or paths")
contract_path <- function(role) {
  i <- match(role, runtime_contract$input_role)
  need(!is.na(i), "Runtime input role is absent: %s", role)
  p <- normalizePath(file.path(repo, runtime_contract$relative_path[i]), mustWork = TRUE)
  need(startsWith(paste0(p, .Platform$file.sep), paste0(repo, .Platform$file.sep)), "Runtime input escapes the repository: %s", role)
  need(identical(sha(p), runtime_contract$expected_sha256[i]), "Runtime input hash differs for role: %s", role)
  p
}
fit_cox <- function(formula, data) {
  z <- survival::coxph(formula, data = data, ties = "efron", x = TRUE, y = TRUE, model = TRUE)
  need(is.null(z$fail) && all(is.finite(coef(z))), "Cox fit failed or contains aliased/nonfinite coefficients: %s", paste(deparse(formula), collapse = " "))
  z
}
cox_rows <- function(fit, model_id, population) {
  sm <- summary(fit)
  cf <- as.data.frame(sm$coefficients); ci <- as.data.frame(sm$conf.int)
  need(identical(rownames(cf), rownames(ci)), "Coefficient/CI row mismatch")
  data.frame(model_id = model_id, population = population, term = rownames(cf), n = fit$n, events = fit$nevent,
             coefficient = cf$coef, se = cf$`se(coef)`, hazard_ratio = ci$`exp(coef)`,
             ci_lower = ci$`lower .95`, ci_upper = ci$`upper .95`, p_value = cf$`Pr(>|z|)`, row.names = NULL)
}
cindex <- function(time, status, marker) {
  z <- survival::concordance(Surv(time, status) ~ marker, reverse = TRUE)
  cnt <- z$count
  if (is.matrix(cnt)) cnt <- colSums(cnt)
  data.frame(harrell_c = as.numeric(z$concordance), comparable_pairs = sum(cnt[c("concordant", "discordant", "tied.x")], na.rm = TRUE))
}
model_metrics <- function(fit, data, model_id, population) {
  lp <- as.numeric(predict(fit, newdata = data, type = "lp"))
  cc <- cindex(data$OS_time_months, data$OS_status, lp)
  data.frame(model_id = model_id, population = population, n = nrow(data), events = sum(data$OS_status),
             model_df = sum(!is.na(coef(fit))), log_likelihood = as.numeric(logLik(fit)), AIC = AIC(fit),
             harrell_c = cc$harrell_c, comparable_pairs = cc$comparable_pairs)
}
stage34 <- function(x) {
  z <- trimws(as.character(x))
  ifelse(grepl("^Stage I([AB])?$|^Stage II([AB])?$", z, ignore.case = TRUE), 0L,
         ifelse(grepl("^Stage III|^Stage IV", z, ignore.case = TRUE), 1L, NA_integer_))
}

# -----------------------------------------------------------------------------
# Authority and score contract. No model gene is named in this script and no
# alternate embedded gene/coefficient vector is permitted.
# -----------------------------------------------------------------------------
model_lock_path <- contract_path("model_lock")
model_root <- dirname(model_lock_path)
downstream_contract_path <- contract_path("downstream_contract")
model_lock <- fromJSON(model_lock_path, simplifyVector = FALSE)
score_path <- normalizePath(file.path(repo, model_lock$authoritative_artifacts$score_table$path), mustWork = TRUE)
coefficient_path <- normalizePath(file.path(repo, model_lock$authoritative_artifacts$coefficients$path), mustWork = TRUE)
required_authority <- c(model_lock_path, score_path, coefficient_path, downstream_contract_path)
need(all(file.exists(required_authority)), "One or more sparse4 authority files are missing")
need(identical(model_lock$status, "VALIDATED_MODEL_LOCK"), "Sparse4 MODEL_LOCK is not validated")
need(length(model_lock$genes_and_coefficients) >= 1L, "Sparse4 MODEL_LOCK has no selected coefficients")
need(identical(sha(score_path), model_lock$authoritative_artifacts$score_table$sha256), "SCORE_TABLE hash differs from MODEL_LOCK")
need(identical(sha(coefficient_path), model_lock$authoritative_artifacts$coefficients$sha256), "COEFFICIENTS hash differs from MODEL_LOCK")

coef_table <- read.delim(coefficient_path, check.names = FALSE)
lock_gene <- vapply(model_lock$genes_and_coefficients, function(z) z$gene_symbol, character(1))
lock_beta <- vapply(model_lock$genes_and_coefficients, function(z) as.numeric(z$coefficient), numeric(1))
need(nrow(coef_table) == length(lock_gene) && !anyDuplicated(coef_table$gene_symbol), "COEFFICIENTS row count or uniqueness differs from MODEL_LOCK")
need(identical(as.character(coef_table$gene_symbol), lock_gene), "COEFFICIENTS gene order differs from MODEL_LOCK")
need(max(abs(as.numeric(coef_table$coefficient) - lock_beta)) <= 1e-15, "COEFFICIENTS values differ from MODEL_LOCK")

score <- read.delim(score_path, check.names = FALSE)
need(identical(names(score), c("sample_id", "patient_id", "OS_time", "OS_status", "linear_predictor")), "Score schema changed")
need(nrow(score) == as.integer(model_lock$roster$patients) && !anyDuplicated(score$sample_id) && !anyDuplicated(score$patient_id) && sum(score$OS_status) == as.integer(model_lock$roster$events), "Score roster/event contract failed")
need(all(is.finite(score$OS_time) & score$OS_time > 0) && all(score$OS_status %in% c(0L, 1L)) && all(is.finite(score$linear_predictor)), "Score values are invalid")
need(all(grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", score$patient_id)), "Score patient IDs are not 12-character TCGA barcodes")
need(identical(substr(score$sample_id, 1L, 12L), score$patient_id), "Score sample-to-patient projection failed")

# -----------------------------------------------------------------------------
# Rebuild the clinical master from guarded clinical inputs.
# -----------------------------------------------------------------------------
survival_path <- contract_path("survival")
se_path <- contract_path("clinical_model_source")
cde_path <- contract_path("anatomic_source")
subsite_dictionary_path <- contract_path("anatomic_dictionary")
hpv_path <- contract_path("hpv")
tx_path <- contract_path("treatment")
tme_path <- contract_path("tme")
clinical_inputs <- c(survival_path, se_path, cde_path, subsite_dictionary_path, hpv_path, tx_path, tme_path)
need(all(file.exists(clinical_inputs)), "One or more clinical inputs are missing")

surv <- read.csv(survival_path, check.names = FALSE)
need(nrow(surv) == as.integer(model_lock$roster$patients) && !anyDuplicated(surv$sample_id) && sum(surv$OS_status) == as.integer(model_lock$roster$events), "Survival authority contract failed")
surv$patient_id <- substr(surv$sample_id, 1L, 12L)
surv$sample15 <- substr(surv$sample_id, 1L, 15L)
idx <- match(score$patient_id, surv$patient_id)
need(!anyNA(idx) && !anyDuplicated(surv$patient_id), "Score-survival patient mapping failed")
surv <- surv[idx, , drop = FALSE]
need(identical(as.character(surv$sample_id), as.character(score$sample_id)), "Score-survival aliquot identity differs")
need(max(abs(surv$OS_time - score$OS_time)) <= 1e-12 && identical(as.integer(surv$OS_status), as.integer(score$OS_status)), "Score endpoint differs from survival authority")

se <- readRDS(se_path)
need(inherits(se, "SummarizedExperiment"), "Clinical model source is not SummarizedExperiment")
cd <- as.data.frame(SummarizedExperiment::colData(se), stringsAsFactors = FALSE)
cd$sample_id <- rownames(cd)
cidx <- match(surv$sample_id, cd$sample_id)
need(!anyNA(cidx), "Clinical model-source join failed")
age <- suppressWarnings(as.numeric(cd$age_at_index[cidx]))
stage <- stage34(cd$ajcc_pathologic_stage[cidx])

cde <- read.csv(cde_path, check.names = FALSE)
need(all(c("CDE_ID.2003301", "CDE_ID.3108203") %in% names(cde)), "CDE subsite fields missing")
cde <- cde[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", cde$CDE_ID.2003301), c("CDE_ID.2003301", "CDE_ID.3108203")]
need(!anyDuplicated(cde$CDE_ID.2003301), "CDE patient IDs duplicated")
sub_idx <- match(score$patient_id, cde$CDE_ID.2003301)
need(!anyNA(sub_idx), "CDE join incomplete")
raw_subsite <- cde$CDE_ID.3108203[sub_idx]
subsite_dictionary <- read.delim(subsite_dictionary_path, check.names = FALSE)
need(all(c("raw_CDE_ID_3108203_value", "frozen_subsite", "model_factor_order") %in% names(subsite_dictionary)), "Frozen subsite dictionary schema changed")
need(!anyDuplicated(subsite_dictionary$raw_CDE_ID_3108203_value), "Frozen subsite dictionary contains duplicate raw values")
subsite_map <- setNames(as.character(subsite_dictionary$frozen_subsite), as.character(subsite_dictionary$raw_CDE_ID_3108203_value))
factor_orders <- unique(as.character(subsite_dictionary$model_factor_order))
need(length(factor_orders) == 1L, "Frozen subsite dictionary has inconsistent factor order")
subsite_levels <- strsplit(factor_orders, ";", fixed = TRUE)[[1L]]
need(all(raw_subsite %in% names(subsite_map)), "Unmapped or missing subsite")
subsite <- factor(unname(subsite_map[raw_subsite]), levels = subsite_levels)

hpv0 <- read.csv(hpv_path, check.names = FALSE)
hpv0 <- hpv0[grepl("^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$", hpv0$bcr_patient_barcode), ]
need(!anyDuplicated(hpv0$bcr_patient_barcode), "HPV patient IDs duplicated")
hidx <- match(score$patient_id, hpv0$bcr_patient_barcode)
hraw <- hpv0$hpv_status_by_p16_testing[hidx]
hpv <- rep("Unknown", nrow(score)); hpv[!is.na(hraw) & tolower(hraw) == "positive"] <- "Positive"; hpv[!is.na(hraw) & tolower(hraw) == "negative"] <- "Negative"
hpv <- factor(hpv, levels = c("Negative", "Positive", "Unknown"))

tx <- read.delim(tx_path, check.names = FALSE)
need(nrow(tx) == as.integer(model_lock$roster$patients) && !anyDuplicated(tx$patient_id), "Treatment authority roster differs from MODEL_LOCK")
tidx <- match(score$patient_id, tx$patient_id); need(!anyNA(tidx), "Treatment join failed")
surgery <- suppressWarnings(as.numeric(tx$surgery_binary[tidx])); tx_status <- tx$treatment_record_status[tidx]
need(all(is.na(surgery[tx_status == "record_absent"])), "Absent treatment records were encoded as surgery=no")
surgery_factor <- factor(ifelse(surgery == 0, "No", ifelse(surgery == 1, "Yes", NA_character_)), levels = c("No", "Yes"))

tme <- read.csv(tme_path, check.names = FALSE)
need(identical(names(tme), c("sample_id", "ImmuneScore", "StromalScore")), "TME schema changed")
need(!anyDuplicated(substr(tme$sample_id, 1L, 15L)), "TME 15-character IDs duplicated")
midx <- match(surv$sample15, substr(tme$sample_id, 1L, 15L)); need(!anyNA(midx), "TME join failed")

score_mean <- mean(score$linear_predictor); score_sd <- sd(score$linear_predictor); score_median <- median(score$linear_predictor)
need(is.finite(score_sd) && score_sd > 0, "Score SD invalid")
master <- data.frame(
  patient_id = score$patient_id, sample_id = surv$sample_id,
  OS_time_months = as.numeric(score$OS_time), OS_status = as.integer(score$OS_status),
  score_raw_linear_predictor = score$linear_predictor,
  score_per_global_sd = (score$linear_predictor - score_mean) / score_sd,
  score_median_group = factor(ifelse(score$linear_predictor >= score_median, "High", "Low"), levels = c("Low", "High")),
  age_years = age,
  stage_factor = factor(ifelse(stage == 0L, "I_II", ifelse(stage == 1L, "III_IV", NA_character_)), levels = c("I_II", "III_IV")),
  HPV_p16_3level = hpv, CDE3108203_subsite = subsite,
  surgery_factor = surgery_factor, treatment_record_status = tx_status,
  project_immune_ssGSEA = as.numeric(tme$ImmuneScore[midx]), project_stromal_ssGSEA = as.numeric(tme$StromalScore[midx])
)
need(nrow(master) == as.integer(model_lock$roster$patients) && !anyDuplicated(master$patient_id) && sum(master$OS_status) == as.integer(model_lock$roster$events), "Clinical master contract failed")
need(abs(mean(master$score_per_global_sd)) <= 1e-12 && abs(sd(master$score_per_global_sd) - 1) <= 1e-12, "Global score standardization failed")
write_tsv(master, "CLINICAL_MASTER_476.tsv")

contract_checks <- data.frame(
  check_id = c("model_lock_status", "model_lock_id", "model_lock_hash", "score_lock_hash", "coefficient_lock_hash", "coefficient_identity", "patient_roster", "event_count", "sample_patient_mapping", "endpoint_mapping", "score_finite", "global_score_z", "retired_score_authority_not_read"),
  status = "PASS",
  observed = c(model_lock$status, model_lock$model_id, sha(model_lock_path), sha(score_path), sha(coefficient_path), sprintf("n=%d;max_abs_beta_diff=%.3g", nrow(coef_table), max(abs(as.numeric(coef_table$coefficient) - lock_beta))), nrow(master), sum(master$OS_status), "exact aliquot and 12-character patient projection", "exact time/status", sum(is.finite(master$score_raw_linear_predictor)), sprintf("mean=%.17g;sd=%.17g", mean(master$score_per_global_sd), sd(master$score_per_global_sd)), "TRUE"),
  acceptance_rule = c("VALIDATED_MODEL_LOCK", model_lock$model_id, sha(model_lock_path), model_lock$authoritative_artifacts$score_table$sha256, model_lock$authoritative_artifacts$coefficients$sha256, sprintf("%d unique full-precision coefficient rows identical to MODEL_LOCK", length(lock_gene)), sprintf("%d unique patients", as.integer(model_lock$roster$patients)), sprintf("%d events", as.integer(model_lock$roster$events)), "exact", "exact score-vs-survival time/status", as.character(as.integer(model_lock$roster$patients)), "mean=0 and sample SD=1 within 1e-12", "only runtime-bound score authority opened")
)
write_tsv(contract_checks, "DATA_CONTRACT_CHECKS.tsv")
score_authority <- data.frame(model_lock_sha256 = sha(model_lock_path), model_id = model_lock$model_id,
                              biological_chain = model_lock$biological_chain,
                              candidate_n = model_lock$candidate_universe$exact_locked_vst_intersection,
                              selected_n = length(model_lock$genes_and_coefficients), score_sha256 = sha(score_path),
                              n = nrow(master), events = sum(master$OS_status), score_mean = score_mean, score_sample_sd = score_sd,
                              score_median = score_median, score_direction = model_lock$score_contract$direction,
                              claim_ceiling = "same-cohort post-lock clinical association and apparent discrimination only; not external validation, calibration, causality, clinical utility, or prospective performance")
write_tsv(score_authority, "SCORE_AUTHORITY.tsv")

# -----------------------------------------------------------------------------
# Common complete-case primary and maximum-available sensitivity M1-M4.
# -----------------------------------------------------------------------------
forms <- list(
  M1 = Surv(OS_time_months, OS_status) ~ score_per_global_sd,
  M2 = Surv(OS_time_months, OS_status) ~ score_per_global_sd + project_immune_ssGSEA,
  M3 = Surv(OS_time_months, OS_status) ~ score_per_global_sd + project_immune_ssGSEA + age_years + stage_factor,
  M4 = Surv(OS_time_months, OS_status) ~ score_per_global_sd + project_immune_ssGSEA + age_years + stage_factor + HPV_p16_3level + CDE3108203_subsite + surgery_factor + project_stromal_ssGSEA
)
required_fields <- list(
  M1 = c("OS_time_months", "OS_status", "score_per_global_sd"),
  M2 = c("OS_time_months", "OS_status", "score_per_global_sd", "project_immune_ssGSEA"),
  M3 = c("OS_time_months", "OS_status", "score_per_global_sd", "project_immune_ssGSEA", "age_years", "stage_factor"),
  M4 = c("OS_time_months", "OS_status", "score_per_global_sd", "project_immune_ssGSEA", "age_years", "stage_factor", "HPV_p16_3level", "CDE3108203_subsite", "surgery_factor", "project_stromal_ssGSEA")
)
native_data <- lapply(required_fields, function(v) droplevels(master[complete.cases(master[, v, drop = FALSE]), , drop = FALSE]))
common <- native_data$M4
need(nrow(common) > 0L && sum(common$OS_status) > 0L, "M4 common roster empty")
native_fits <- Map(fit_cox, forms, native_data)
common_fits <- lapply(forms, fit_cox, data = common)
native_coef <- do.call(rbind, lapply(names(forms), function(id) cox_rows(native_fits[[id]], id, "maximum-available sensitivity")))
common_coef <- do.call(rbind, lapply(names(forms), function(id) cox_rows(common_fits[[id]], id, "common M4 complete-case primary")))
native_metrics <- do.call(rbind, lapply(names(forms), function(id) model_metrics(native_fits[[id]], native_data[[id]], id, "maximum-available sensitivity")))
common_metrics <- do.call(rbind, lapply(names(forms), function(id) model_metrics(common_fits[[id]], common, id, "common M4 complete-case primary")))
write_tsv(native_coef, "M1_M4_MAXIMUM_AVAILABLE_COEFFICIENTS.tsv")
write_tsv(native_metrics, "M1_M4_MAXIMUM_AVAILABLE_METRICS.tsv")
write_tsv(common_coef, "M1_M4_COMMON_ROSTER_COEFFICIENTS.tsv")
write_tsv(common_metrics, "M1_M4_COMMON_ROSTER_METRICS.tsv")

flow <- do.call(rbind, lapply(names(required_fields), function(id) {
  k <- complete.cases(master[, required_fields[[id]], drop = FALSE])
  data.frame(model_id = id, source_n = nrow(master), retained_n = sum(k), excluded_n = sum(!k),
             retained_events = sum(master$OS_status[k]), excluded_events = sum(master$OS_status[!k]),
             required_fields = paste(required_fields[[id]], collapse = ";"))
}))
write_tsv(flow, "COMPLETE_CASE_FLOW.tsv")

score_rows_common <- common_coef[common_coef$term == "score_per_global_sd", ]
score_rows_native <- native_coef[native_coef$term == "score_per_global_sd", ]
primary_score <- rbind(transform(score_rows_common, estimand = "common-roster primary"), transform(score_rows_native, estimand = "maximum-available sensitivity"))
is_primary_family <- primary_score$estimand == "common-roster primary" & primary_score$model_id %in% c("M1", "M4")
primary_score$multiplicity_family <- ifelse(is_primary_family, "REVISION_PRIMARY_3", ifelse(primary_score$estimand == "maximum-available sensitivity", "NONE_SENSITIVITY", "M1_M4_MODEL_SEQUENCE_DESCRIPTIVE"))
primary_score$adjustment_status <- ifelse(is_primary_family, "PENDING_THIRD_PRESPECIFIED_COMPONENT; no two-only adjustment", ifelse(primary_score$estimand == "maximum-available sensitivity", "SENSITIVITY_NOT_A_PRIMARY_FAMILY_MEMBER", "DESCRIPTIVE_NO_CLAIM"))
write_tsv(primary_score, "CONTINUOUS_SCORE_PRIMARY_AND_SENSITIVITY.tsv")

# Median split is descriptive only and does not replace the continuous estimand.
km_fit <- survfit(Surv(OS_time_months, OS_status) ~ score_median_group, data = master)
km_tab <- summary(km_fit)$table
km_groups <- data.frame(group = sub("score_median_group=", "", rownames(km_tab)), n = km_tab[, "records"], events = km_tab[, "events"],
                        median_survival_months = km_tab[, "median"], ci_lower_95 = km_tab[, "0.95LCL"], ci_upper_95 = km_tab[, "0.95UCL"], row.names = NULL)
km_cox <- fit_cox(Surv(OS_time_months, OS_status) ~ score_median_group, master)
km_sm <- summary(km_cox)
km_effect <- data.frame(comparison = "High versus Low", n = nrow(master), events = sum(master$OS_status), cutoff_raw_global_median = score_median,
                        hazard_ratio = km_sm$conf.int[1L, "exp(coef)"], ci_lower = km_sm$conf.int[1L, "lower .95"], ci_upper = km_sm$conf.int[1L, "upper .95"],
                        nominal_cox_p = km_sm$coefficients[1L, "Pr(>|z|)"],
                        logrank_p = pchisq(survdiff(Surv(OS_time_months, OS_status) ~ score_median_group, master)$chisq, 1, lower.tail = FALSE),
                        role = "DESCRIPTIVE_ONLY; no optimized cutpoint; continuous score is primary")
write_tsv(km_groups, "MEDIAN_KM_GROUP_SUMMARY.tsv"); write_tsv(km_effect, "MEDIAN_KM_DESCRIPTIVE_EFFECT.tsv")

# -----------------------------------------------------------------------------
# Missingness and M4 inclusion comparison (existing C_MISSINGNESS family).
# -----------------------------------------------------------------------------
keep_m4 <- complete.cases(master[, required_fields$M4, drop = FALSE])
master$m4_included <- factor(ifelse(keep_m4, "Included", "Excluded"), levels = c("Included", "Excluded"))
inc <- master[keep_m4, ]; exc <- master[!keep_m4, ]
missingness <- do.call(rbind, lapply(required_fields$M4, function(v) data.frame(variable = v, n_missing = sum(is.na(master[[v]])), fraction_missing = mean(is.na(master[[v]])))))
write_tsv(missingness, "MISSINGNESS_BY_M4_FIELD.tsv")
tests <- list()
add_cont <- function(v) {
  a0 <- inc[[v]]; b0 <- exc[[v]]; a <- a0[is.finite(a0)]; b <- b0[is.finite(b0)]
  need(length(a) >= 2 && length(b) >= 2, "Too few observed values for inclusion comparison: %s", v)
  psd <- sqrt(((length(a)-1)*var(a)+(length(b)-1)*var(b))/(length(a)+length(b)-2))
  data.frame(variable=v, test="two-sided Wilcoxon", included_total_n=length(a0), included_observed_n=length(a), included_missing_n=sum(!is.finite(a0)),
             excluded_total_n=length(b0), excluded_observed_n=length(b), excluded_missing_n=sum(!is.finite(b0)),
             included_summary=sprintf("median=%.6g;IQR=%.6g",median(a),IQR(a)), excluded_summary=sprintf("median=%.6g;IQR=%.6g",median(b),IQR(b)),
             effect=(mean(b)-mean(a))/psd, effect_definition="standardized mean difference: excluded minus included", p_value=wilcox.test(a,b,exact=FALSE)$p.value)
}
for (v in c("age_years","score_raw_linear_predictor","project_immune_ssGSEA","project_stromal_ssGSEA","OS_time_months")) tests[[length(tests)+1L]] <- add_cont(v)
pt <- prop.test(c(sum(inc$OS_status),sum(exc$OS_status)),c(nrow(inc),nrow(exc)),correct=FALSE)
tests[[length(tests)+1L]] <- data.frame(variable="OS_status",test="two-sided score test for proportions",included_total_n=nrow(inc),included_observed_n=nrow(inc),included_missing_n=0,
                                       excluded_total_n=nrow(exc),excluded_observed_n=nrow(exc),excluded_missing_n=0,
                                       included_summary=sprintf("%d/%d",sum(inc$OS_status),nrow(inc)),excluded_summary=sprintf("%d/%d",sum(exc$OS_status),nrow(exc)),
                                       effect=mean(exc$OS_status)-mean(inc$OS_status),effect_definition="event proportion: excluded minus included",p_value=pt$p.value)
add_cat <- function(v) {
  x <- addNA(master[[v]], ifany=TRUE); tab <- table(master$m4_included,x); chi <- suppressWarnings(chisq.test(tab,correct=FALSE))
  cramer <- sqrt(as.numeric(chi$statistic)/(sum(tab)*min(nrow(tab)-1,ncol(tab)-1)))
  if(any(chi$expected<5)){tt<-fisher.test(tab);tn<-"two-sided Fisher exact"}else{tt<-chi;tn<-"two-sided Pearson chi-square"}
  data.frame(variable=v,test=tn,included_total_n=nrow(inc),included_observed_n=sum(!is.na(inc[[v]])),included_missing_n=sum(is.na(inc[[v]])),
             excluded_total_n=nrow(exc),excluded_observed_n=sum(!is.na(exc[[v]])),excluded_missing_n=sum(is.na(exc[[v]])),
             included_summary=paste(names(table(x[keep_m4])),as.integer(table(x[keep_m4])),sep="=",collapse=";"),
             excluded_summary=paste(names(table(x[!keep_m4])),as.integer(table(x[!keep_m4])),sep="=",collapse=";"),
             effect=cramer,effect_definition="Cramer's V",p_value=tt$p.value)
}
for(v in c("stage_factor","HPV_p16_3level","CDE3108203_subsite","surgery_factor")) tests[[length(tests)+1L]]<-add_cat(v)
inclusion <- do.call(rbind,tests); inclusion$p_bh_C_MISSINGNESS <- p.adjust(inclusion$p_value,"BH")
write_tsv(inclusion,"M4_INCLUSION_COMPARISON.tsv")
event_diff <- mean(exc$OS_status)-mean(inc$OS_status); event_se <- sqrt(mean(exc$OS_status)*(1-mean(exc$OS_status))/nrow(exc)+mean(inc$OS_status)*(1-mean(inc$OS_status))/nrow(inc))
write_tsv(data.frame(included_n=nrow(inc),included_events=sum(inc$OS_status),excluded_n=nrow(exc),excluded_events=sum(exc$OS_status),
                     excluded_minus_included_event_difference=event_diff,ci_lower_wald=event_diff-1.96*event_se,ci_upper_wald=event_diff+1.96*event_se,
                     limitation_trigger=ifelse(abs(event_diff)>=0.10,"TRIGGER_SELECTION_BIAS_LIMITATION","NO_10_PERCENTAGE_POINT_TRIGGER")),"M4_INCLUSION_EVENT_DIFFERENCE.tsv")

# -----------------------------------------------------------------------------
# PH and functional-form checks on the common roster.
# -----------------------------------------------------------------------------
ph_list <- lapply(names(common_fits), function(id) {
  z <- cox.zph(common_fits[[id]], transform="km", terms=TRUE, singledf=TRUE)
  tab <- as.data.frame(z$table); tab$term <- rownames(tab); rownames(tab)<-NULL
  data.frame(model_id=id,term=tab$term,chisq=tab$chisq,df=tab$df,p_nominal=tab$p,transform="km")
}); names(ph_list)<-names(common_fits)
ph <- do.call(rbind,ph_list); ph$p_bh_within_model_terms <- NA_real_; ph$p_bh_across_model_global <- NA_real_
for(id in names(common_fits)){ii<-ph$model_id==id & ph$term!="GLOBAL";ph$p_bh_within_model_terms[ii]<-p.adjust(ph$p_nominal[ii],"BH")}
gi<-ph$term=="GLOBAL";ph$p_bh_across_model_global[gi]<-p.adjust(ph$p_nominal[gi],"BH")
write_tsv(ph,"PH_DIAGNOSTICS_COMMON_ROSTER.tsv")

tv_rows <- lapply(names(forms),function(id){
  p0<-ph$p_nominal[ph$model_id==id & ph$term=="score_per_global_sd"]
  need(length(p0)==1,"Score PH row missing: %s",id)
  if(p0>=0.05) return(data.frame(model_id=id,trigger_p=p0,status="N/A_WITH_REASON",beta_main=NA,se_main=NA,gamma_log_time=NA,se_gamma=NA,gamma_p=NA,hr_1m=NA,hr_12m=NA,hr_60m=NA,reason="score term nominal cox.zph p>=0.05"))
  f<-update(forms[[id]],. ~ . + tt(score_per_global_sd)); z<-coxph(f,data=common,ties="efron",tt=function(x,t,...)x*log(pmax(t,1/30.44)),x=TRUE)
  sm<-summary(z); vv<-vcov(z); b<-coef(z)["score_per_global_sd"]; g<-coef(z)["tt(score_per_global_sd)"]
  h<-c(1,12,60); lt<-log(h); bh<-b+g*lt
  data.frame(model_id=id,trigger_p=p0,status="COMPUTED_SENSITIVITY",beta_main=b,se_main=sm$coefficients["score_per_global_sd","se(coef)"],
             gamma_log_time=g,se_gamma=sm$coefficients["tt(score_per_global_sd)","se(coef)"],gamma_p=sm$coefficients["tt(score_per_global_sd)","Pr(>|z|)"],
             hr_1m=exp(bh[1]),hr_12m=exp(bh[2]),hr_60m=exp(bh[3]),reason="score PH nominal p<0.05; score-by-log(time) sensitivity")
})
tv <- do.call(rbind,tv_rows); tv$p_bh_C_PH_SCORE_SENSITIVITY_4 <- NA_real_; fi<-is.finite(tv$gamma_p); tv$p_bh_C_PH_SCORE_SENSITIVITY_4[fi]<-p.adjust(tv$gamma_p[fi],"BH",n=4)
write_tsv(tv,"PH_SCORE_TIME_VARYING_SENSITIVITY.tsv")

need(requireNamespace("Hmisc",quietly=TRUE),"Hmisc required for nonlinear check")
nonlin <- lapply(names(forms),function(id){
  x<-common$score_per_global_sd; knots<-as.numeric(quantile(x,c(.1,.5,.9),type=7)); basis<-Hmisc::rcspline.eval(x,knots=knots,inclx=TRUE)
  d<-common;d$.score_nonlin<-basis[,2];ef<-update(forms[[id]],. ~ . + .score_nonlin)
  mm<-model.matrix(ef,d);df_total<-ncol(mm)-1;rank_eff<-qr(mm[,-1,drop=FALSE])$rank;epv<-sum(d$OS_status)/df_total
  base<-data.frame(model_id=id,n=nrow(d),events=sum(d$OS_status),knot_q10=knots[1],knot_q50=knots[2],knot_q90=knots[3],total_df=df_total,effective_df=rank_eff,events_per_df=epv,
                   delta_lr_chisq=NA,nonlinear_p_nominal=NA,status="NOT_FEASIBLE_WITH_EVIDENCE",reason="")
  if(length(unique(knots))<3){base$reason<-"non-distinct knots";return(base)}
  if(df_total!=rank_eff){base$reason<-"rank-deficient expanded design";return(base)}
  if(epv<10){base$reason<-"prespecified events-per-total-df gate <10";return(base)}
  efit<-fit_cox(ef,d);aa<-anova(common_fits[[id]],efit,test="LRT");base$delta_lr_chisq<-aa[2,"Chisq"];base$nonlinear_p_nominal<-aa[2,"Pr(>|Chi|)"];base$status<-"COMPUTED";base$reason<-"three-knot RCS versus linear, 1-df nonlinear component";base
})
nonlin<-do.call(rbind,nonlin);nonlin$nonlinear_p_bh_C_FUNCTIONAL_FORM_4<-NA_real_;ni<-is.finite(nonlin$nonlinear_p_nominal);nonlin$nonlinear_p_bh_C_FUNCTIONAL_FORM_4[ni]<-p.adjust(nonlin$nonlinear_p_nominal[ni],"BH",n=4)
write_tsv(nonlin,"SCORE_NONLINEARITY_COMMON_ROSTER.tsv")

# -----------------------------------------------------------------------------
# Clinical incremental information, non-bootstrap paired C, and risk-set AUC.
# -----------------------------------------------------------------------------
f_cb <- Surv(OS_time_months,OS_status) ~ age_years + stage_factor + HPV_p16_3level + CDE3108203_subsite + surgery_factor
f_cb_score <- update(f_cb,. ~ . + score_per_global_sd)
f_eb <- update(f_cb,. ~ . + project_immune_ssGSEA + project_stromal_ssGSEA)
f_eb_score <- update(f_eb,. ~ . + score_per_global_sd)
inc_forms <- list(CB=f_cb,CB_SCORE=f_cb_score,EB=f_eb,EB_SCORE=f_eb_score)
inc_fits <- lapply(inc_forms,fit_cox,data=common)
inc_metrics <- do.call(rbind,lapply(names(inc_fits),function(id)model_metrics(inc_fits[[id]],common,id,"common M4 roster")))
comp <- data.frame(comparison_id=c("CB_SCORE_vs_CB","EB_SCORE_vs_EB"),base=c("CB","EB"),with_score=c("CB_SCORE","EB_SCORE"))
lrt <- do.call(rbind,lapply(seq_len(nrow(comp)),function(i){a<-anova(inc_fits[[comp$base[i]]],inc_fits[[comp$with_score[i]]],test="LRT");data.frame(comparison_id=comp$comparison_id[i],n=nrow(common),events=sum(common$OS_status),delta_lr_chisq=a[2,"Chisq"],delta_df=a[2,"Df"],p_value=a[2,"Pr(>|Chi|)"],delta_AIC_base_minus_score=AIC(inc_fits[[comp$base[i]]])-AIC(inc_fits[[comp$with_score[i]]]))}))
lrt$p_holm_C_INCREMENTAL_2<-p.adjust(lrt$p_value,"holm");write_tsv(inc_metrics,"CLINICAL_INCREMENTAL_MODEL_METRICS.tsv");write_tsv(lrt,"CLINICAL_INCREMENTAL_LR_AIC.tsv")

need(requireNamespace("compareC",quietly=TRUE),"compareC required for paired non-bootstrap C comparison")
paired_c<-do.call(rbind,lapply(seq_len(nrow(comp)),function(i){
  lp1<-as.numeric(predict(inc_fits[[comp$with_score[i]]],newdata=common,type="lp"));lp0<-as.numeric(predict(inc_fits[[comp$base[i]]],newdata=common,type="lp"))
  z<-compareC::compareC(common$OS_time_months,common$OS_status,-lp1,-lp0)
  c1<-cindex(common$OS_time_months,common$OS_status,lp1)$harrell_c;c0<-cindex(common$OS_time_months,common$OS_status,lp0)$harrell_c
  need(max(abs(c(c1-unname(z$est.c[1]),c0-unname(z$est.c[2]))))<1e-10,"compareC direction mismatch")
  data.frame(comparison_id=comp$comparison_id[i],n=nrow(common),events=sum(common$OS_status),c_with_score=c1,c_base=c0,delta_c=c1-c0,variance_delta=z$est.vardiff_c,z=z$zscore,p_value=z$pval,method="compareC analytic paired covariance; no bootstrap")
}))
paired_c$p_holm_C_INCREMENTAL_C_2<-p.adjust(paired_c$p_value,"holm");write_tsv(paired_c,"CLINICAL_INCREMENTAL_PAIRED_C_NONBOOTSTRAP.tsv")

horizon_gate<-do.call(rbind,lapply(list(SCORE_ONLY_FULL_476=master,CLINICAL_COMMON_320=common),function(d)do.call(rbind,lapply(c(36,60),function(h)data.frame(population=if(nrow(d)==476L)"SCORE_ONLY_FULL_476" else "CLINICAL_COMMON_320",horizon_months=h,n=nrow(d),events_total=sum(d$OS_status),at_risk_at_horizon=sum(d$OS_time_months>=h),events_by_horizon=sum(d$OS_status==1 & d$OS_time_months<=h))))))
horizon_gate$gate_rule<-"at_risk_at_horizon >= 20 AND events_by_horizon >= 20";horizon_gate$gate_status<-ifelse(horizon_gate$at_risk_at_horizon>=20 & horizon_gate$events_by_horizon>=20,"PASS","NOT_FEASIBLE_WITH_EVIDENCE")
write_tsv(horizon_gate,"TIME_AUC_HORIZON_GATES.tsv")
need(requireNamespace("timeROC",quietly=TRUE),"timeROC required for point AUC")
auc_rows<-list();k<-0L
auc_specs<-c(list(SCORE_FULL_476=list(data=master,marker=master$score_per_global_sd,population="SCORE_ONLY_FULL_476")),lapply(names(inc_fits),function(id)list(data=common,marker=as.numeric(predict(inc_fits[[id]],newdata=common,type="lp")),population="CLINICAL_COMMON_320")))
names(auc_specs)<-c("SCORE_FULL_476",names(inc_fits))
for(h in c(36,60))for(id in names(auc_specs)){k<-k+1L;spec<-auc_specs[[id]];dauc<-spec$data;gi<-horizon_gate$population==spec$population & horizon_gate$horizon_months==h;need(sum(gi)==1L,"time-AUC horizon gate key mismatch");gate<-horizon_gate$gate_status[gi];val<-NA_real_;status<-gate
  if(gate=="PASS"){z<-timeROC::timeROC(T=dauc$OS_time_months,delta=dauc$OS_status,marker=spec$marker,cause=1,weighting="marginal",times=h,iid=FALSE);ii<-which(abs(z$times-h)<sqrt(.Machine$double.eps));need(length(ii)==1,"timeROC horizon missing");val<-z$AUC[ii];status<-"COMPUTED_DESCRIPTIVE"}
  auc_rows[[k]]<-data.frame(model_id=id,population=spec$population,horizon_months=h,n=nrow(dauc),events=sum(dauc$OS_status),at_risk_at_horizon=horizon_gate$at_risk_at_horizon[gi],events_by_horizon=horizon_gate$events_by_horizon[gi],apparent_ipcw_auc=val,status=status,multiplicity_family="C10_TIME_AUC_DESCRIPTIVE",claim_rule="DESCRIPTIVE_ONLY; apparent point estimate; no interval or horizon superiority test")}
write_tsv(do.call(rbind,auc_rows),"TIME_AUC_36_60_RISK_SET_GATED.tsv")

# The historical performance branch used resampling-based calibration. New
# resampling is prospectively prohibited, and no alternative calibration
# uncertainty estimand was frozen before inspection. Preserve this as an
# explicit evidence-supported disposition rather than silently reporting an
# apparent curve without uncertainty.
calibration_disposition<-data.frame(
  analysis_id=c("CALIBRATION_36M","CALIBRATION_60M"),
  status="NOT_FEASIBLE_WITH_EVIDENCE",
  reason="Historical method requires new patient resampling, which is prospectively prohibited; no non-resampling calibration uncertainty method was frozen before outcome inspection",
  replacement="none",
  claim_boundary="No calibration, optimism-correction, or clinical-utility claim from this branch",
  stringsAsFactors=FALSE)
write_tsv(calibration_disposition,"CALIBRATION_DISPOSITION.tsv")

# -----------------------------------------------------------------------------
# Existing stage-only MI contract, executed only when the preconditions allow.
# -----------------------------------------------------------------------------
mi_nonstage<-c("OS_time_months","OS_status","score_per_global_sd","project_immune_ssGSEA","age_years","HPV_p16_3level","CDE3108203_subsite","surgery_factor","project_stromal_ssGSEA")
dmi<-droplevels(master[complete.cases(master[,mi_nonstage,drop=FALSE]),c(mi_nonstage,"stage_factor"),drop=FALSE])
mi_warnings <- character()
mi_available <- withCallingHandlers(requireNamespace("mice",quietly=TRUE), warning=function(w){mi_warnings<<-c(mi_warnings,conditionMessage(w));invokeRestart("muffleWarning")})
need(mi_available, "mice is required for the registered stage-only MI sensitivity")
mi_body <- paste(deparse(body(get("mice.impute.logreg",envir=asNamespace("mice")))),collapse="\n")
uses_glm_fit <- grepl("glm.fit",mi_body,fixed=TRUE); invokes_mixed <- grepl("glmer|lmer",mi_body)
matrix_abi <- as.character(Matrix::Matrix.Version()[["abi"]]); lme4_abi <- as.character(get(".Matrix.abi.build.version",envir=asNamespace("lme4"),inherits=FALSE))
backend_ok <- uses_glm_fit && !invokes_mixed
mi_backend <- data.frame(method="logreg",imputation_function="mice::mice.impute.logreg",function_body_sha256=digest::digest(mi_body,algo="sha256",serialize=TRUE),
                         backend="stats::glm.fit with quasibinomial(logit)",uses_glm_fit=uses_glm_fit,invokes_lme4_mixed_model=invokes_mixed,
                         Matrix_ABI_current=matrix_abi,lme4_Matrix_ABI_built=lme4_abi,lme4_Matrix_ABI_match=identical(matrix_abi,lme4_abi),
                         abi_disposition=ifelse(identical(matrix_abi,lme4_abi),"MATCH","MISMATCH_NOT_IN_STAGE_LOGREG_EXECUTION_PATH"),
                         namespace_warnings=ifelse(length(mi_warnings),gsub("[\r\n]+"," ",paste(unique(mi_warnings),collapse=" | ")),"none"),backend_status=ifelse(backend_ok,"PASS_EXACT_LOGREG_BACKEND","FAIL"))
write_tsv(mi_backend,"MI_BACKEND_CONTRACT.tsv")
need(backend_ok,"mice logreg backend is not the registered glm.fit route")
mi_pre_gate<-nrow(dmi)>0 && sum(dmi$OS_status)>0 && sum(is.na(dmi$stage_factor))>0 && length(unique(na.omit(dmi$stage_factor)))==2L && backend_ok
if(mi_pre_gate){
  dmi$nelson_aalen<-mice::nelsonaalen(dmi,OS_time_months,OS_status);method<-mice::make.method(dmi);method[]<-"";method["stage_factor"]<-"logreg"
  pred<-mice::make.predictorMatrix(dmi);pred[,]<-0L;pred["stage_factor",c("age_years","CDE3108203_subsite","score_per_global_sd","OS_status","nelson_aalen")]<-1L
  set.seed(66001L);imp<-mice::mice(dmi,m=20,maxit=20,method=method,predictorMatrix=pred,seed=66001,printFlag=FALSE)
  imputed<-names(imp$imp)[vapply(imp$imp,nrow,integer(1))>0];need(identical(imputed,"stage_factor"),"A variable other than stage was imputed")
  mi_fit<-with(imp,coxph(Surv(OS_time_months,OS_status)~score_per_global_sd+project_immune_ssGSEA+age_years+stage_factor+HPV_p16_3level+CDE3108203_subsite+surgery_factor+project_stromal_ssGSEA))
  pooled<-summary(mice::pool(mi_fit),conf.int=TRUE)
  logs<-if(is.null(imp$loggedEvents))0L else nrow(imp$loggedEvents);conv<-mice::convergence(imp);names(conv)[names(conv)==".it"]<-"iteration";names(conv)[names(conv)=="vrb"]<-"variable";conv<-conv[conv$variable=="stage_factor",]
  last5<-conv$iteration>=16;finite<-all(is.finite(conv$ac[last5]))&&all(is.finite(conv$psrf[last5]));mxac<-if(finite)max(abs(conv$ac[last5])) else Inf;mxpsrf<-if(finite)max(abs(conv$psrf[last5]-1)) else Inf
  chain_mean<-imp$chainMean["stage_factor",,,drop=TRUE];chain_var<-imp$chainVar["stage_factor",,,drop=TRUE]
  gate<-all(is.finite(chain_mean))&&all(is.finite(chain_var))&&finite&&mxac<=.5&&mxpsrf<=.1&&logs==0
  diag<-data.frame(pre_gate="PASS",n=nrow(dmi),events=sum(dmi$OS_status),stage_missing=sum(is.na(dmi$stage_factor)),imputations=20,iterations=20,seed=66001,imputed_variables=paste(imputed,collapse=";"),logged_events=logs,max_abs_autocorrelation_last5=mxac,max_abs_psrf_minus_1_last5=mxpsrf,backend_status=mi_backend$backend_status,lme4_Matrix_ABI_match=mi_backend$lme4_Matrix_ABI_match,abi_disposition=mi_backend$abi_disposition,convergence_gate_status=ifelse(gate,"PASS","NOT_FEASIBLE_WITH_EVIDENCE"),gate_rule="exact glm.fit logreg backend; finite chains; last5 |AC|<=0.50; last5 |PSRF-1|<=0.10; zero logged events")
  pooled_out<-data.frame(term=pooled$term,estimate_log_hr=pooled$estimate,se=pooled$std.error,df=pooled$df,p_value=pooled$p.value,hazard_ratio=exp(pooled$estimate),ci_lower=exp(pooled$`2.5 %`),ci_upper=exp(pooled$`97.5 %`),status=ifelse(gate,"COMPUTED_SENSITIVITY","NOT_FEASIBLE_WITH_EVIDENCE"))
  per_imp<-do.call(rbind,lapply(seq_along(mi_fit$analyses),function(i){s<-summary(mi_fit$analyses[[i]]);z<-s$coefficients["score_per_global_sd",];data.frame(imputation_id=i,n=mi_fit$analyses[[i]]$n,events=mi_fit$analyses[[i]]$nevent,coefficient=z["coef"],se=z["se(coef)"],hazard_ratio=exp(z["coef"]),p_value=z["Pr(>|z|)"],status=ifelse(gate,"COMPUTED_SENSITIVITY","NOT_FEASIBLE_WITH_EVIDENCE"))}))
}else{
  diag<-data.frame(pre_gate="NOT_FEASIBLE_WITH_EVIDENCE",n=nrow(dmi),events=sum(dmi$OS_status),stage_missing=sum(is.na(dmi$stage_factor)),imputations=0,iterations=0,seed=66001,imputed_variables="",logged_events=NA,max_abs_autocorrelation_last5=NA,max_abs_psrf_minus_1_last5=NA,backend_status=mi_backend$backend_status,lme4_Matrix_ABI_match=mi_backend$lme4_Matrix_ABI_match,abi_disposition=mi_backend$abi_disposition,convergence_gate_status="NOT_RUN",gate_rule="exact glm.fit logreg backend; nonstage-complete population; events>0; stage missing>0; two observed stage levels")
  pooled_out<-data.frame(term="score_per_global_sd",estimate_log_hr=NA,se=NA,df=NA,p_value=NA,hazard_ratio=NA,ci_lower=NA,ci_upper=NA,status="NOT_FEASIBLE_WITH_EVIDENCE");per_imp<-data.frame()
}
write_tsv(diag,"STAGE_ONLY_MI_DIAGNOSTICS.tsv");write_tsv(pooled_out,"STAGE_ONLY_MI_POOLED.tsv");if(nrow(per_imp))write_tsv(per_imp,"STAGE_ONLY_MI_PER_IMPUTATION.tsv")

multiplicity<-data.frame(
  family=c("REVISION_PRIMARY_3","C_INCREMENTAL_2","C_INCREMENTAL_C_2","C_FUNCTIONAL_FORM_4","C_PH_WITHIN_MODEL_TERMS","C_PH_GLOBAL_4","C_PH_SCORE_SENSITIVITY_4","C_MISSINGNESS_10","C10_TIME_AUC_DESCRIPTIVE","KM_MEDIAN_DESCRIPTIVE","STAGE_MI_SENSITIVITY"),
  correction=c("Holm pending third component","Holm","Holm","BH with n=4","BH separately within each model","BH across four global tests","BH with n=4","BH","NONE_DESCRIPTIVE","NONE_DESCRIPTIVE","NONE_SINGLE_SENSITIVITY"),
  status=c("INCOMPLETE; nominal P1/P2 retained","COMPLETE","COMPLETE","COMPLETE_WITH_NOT_FEASIBLE_ROWS_RETAINED","COMPLETE","COMPLETE","COMPLETE_WITH_NA_WHEN_NOT_TRIGGERED","COMPLETE","DESCRIPTIVE_NO_PVALUES","DESCRIPTIVE_NOMINAL_ONLY","SENSITIVITY_ONLY"))
write_tsv(multiplicity,"MULTIPLICITY_REGISTRY.tsv")

session <- capture.output(sessionInfo());writeLines(session,file.path(out,"SESSIONINFO.txt"))
receipt_out<-list(schema="V31_HISTORICAL_YELLOW_SPARSE4_TCGA_CLINICAL_BRANCH_V1",status="COMPUTED_PENDING_INDEPENDENT_AUDIT",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),
                  model_id=model_lock$model_id,model_lock_sha256=sha(model_lock_path),score_sha256=sha(score_path),coefficient_sha256=sha(coefficient_path),
                  runtime_input_contract_sha256=sha(file.path(out,"RUNTIME_INPUT_CONTRACT.tsv")),subsite_dictionary_sha256=sha(subsite_dictionary_path),runner_sha256=sha(runner_path),
                  patients=nrow(master),events=sum(master$OS_status),score_gene_count=length(model_lock$genes_and_coefficients),
                  common_m4_n=nrow(common),common_m4_events=sum(common$OS_status),native_sizes=lapply(names(native_data),function(id)list(model=id,n=nrow(native_data[[id]]),events=sum(native_data[[id]]$OS_status))),
                  resampling_executed=FALSE,plots_created=FALSE,model_lock_consumed=TRUE,model_lock_mutated=FALSE,
                  claim_ceiling="same-cohort clinical association and apparent descriptive discrimination only; not external validation, calibration, causality, clinical utility, or prospective performance")
write_json(receipt_out,"ANALYSIS_RECEIPT.json")
cat(toJSON(list(status=receipt_out$status,model_id=receipt_out$model_id,n=receipt_out$patients,events=receipt_out$events,common_n=receipt_out$common_m4_n,common_events=receipt_out$common_m4_events,resampling=FALSE,plots=FALSE),auto_unbox=TRUE),"\n")
