#!/usr/bin/env Rscript

# ==========================================================
# Script   : 11_sparse4_internal_validation.R
# Purpose  : Compute the fold-comparable conditional fixed-OOF discrimination estimand
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: internal validation; Fig. 4H; Supplementary Tables 9 and 14
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/06_oof/run_oof_estimand_repair.R
# ==========================================================

suppressPackageStartupMessages(library(survival))

root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "sparse4_oof_estimand_repair_20261020")
dir.create(file.path(out, "results"), recursive = TRUE, showWarnings = FALSE)

pred_path <- file.path(root, "results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_HELDOUT_PREDICTIONS.tsv")
mean_path <- file.path(root, "results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_PATIENT_MEAN_PREDICTIONS.tsv")
task_path <- file.path(root, "results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_TASK_METADATA.tsv")
stopifnot(file.exists(pred_path), file.exists(mean_path), file.exists(task_path))

pred <- read.delim(pred_path, check.names = FALSE, stringsAsFactors = FALSE)
pm <- read.delim(mean_path, check.names = FALSE, stringsAsFactors = FALSE)
task <- read.delim(task_path, check.names = FALSE, stringsAsFactors = FALSE)

stopifnot(nrow(pred) == 9520L, nrow(pm) == 476L, nrow(task) == 100L)
stopifnot(!anyDuplicated(pred[c("patient_id", "repeat_id")]))
stopifnot(length(unique(pred$patient_id)) == 476L,
          length(unique(pred$repeat_id)) == 20L,
          all(table(pred$repeat_id) == 476L),
          all(table(pred$repeat_id, pred$outer_fold) > 0L))
stopifnot(all(pred$OS_status %in% c(0, 1)), !anyNA(pred$linear_predictor))

fold_counts <- list()
repeat_rows <- list()
for (rid in sort(unique(pred$repeat_id))) {
  z <- pred[pred$repeat_id == rid, , drop = FALSE]
  counts <- setNames(rep(0, 5), c("concordant", "discordant", "tied.x", "tied.y", "tied.xy"))
  for (fid in sort(unique(z$outer_fold))) {
    q <- z[z$outer_fold == fid, , drop = FALSE]
    cc <- concordance(Surv(OS_time, OS_status) ~ linear_predictor,
                      data = q, reverse = TRUE)
    counts <- counts + cc$count
    fold_counts[[length(fold_counts) + 1L]] <- data.frame(
      repeat_id = rid, outer_fold = fid, n = nrow(q), events = sum(q$OS_status),
      concordant = unname(cc$count["concordant"]),
      discordant = unname(cc$count["discordant"]),
      tied_x = unname(cc$count["tied.x"]),
      harrell_c = unname(cc$concordance), stringsAsFactors = FALSE)
  }
  denom <- counts["concordant"] + counts["discordant"] + counts["tied.x"]
  pooled_c <- (counts["concordant"] + 0.5 * counts["tied.x"]) / denom

  z$lp_fold_z <- ave(z$linear_predictor, z$outer_fold,
                     FUN = function(x) as.numeric(scale(x)))
  stopifnot(!anyNA(z$lp_fold_z), all(is.finite(z$lp_fold_z)))
  fit <- coxph(Surv(OS_time, OS_status) ~ lp_fold_z + strata(outer_fold),
               data = z, ties = "efron", x = TRUE)
  sm <- summary(fit)$coefficients["lp_fold_z", ]
  ci <- summary(fit)$conf.int["lp_fold_z", ]
  ph <- cox.zph(fit, transform = "km")$table["lp_fold_z", ]
  repeat_rows[[length(repeat_rows) + 1L]] <- data.frame(
    repeat_id = rid, n = nrow(z), events = sum(z$OS_status),
    fold_stratified_harrell_c = unname(pooled_c),
    log_hr_per_fold_sd = unname(coef(fit)["lp_fold_z"]),
    hr_per_fold_sd = unname(ci["exp(coef)"]),
    ci_low_model_based = unname(ci["lower .95"]),
    ci_high_model_based = unname(ci["upper .95"]),
    nominal_p_descriptive = unname(sm["Pr(>|z|)"]),
    ph_p = unname(ph["p"]), stringsAsFactors = FALSE)
}

fold_df <- do.call(rbind, fold_counts)
rep_df <- do.call(rbind, repeat_rows)

pm$lp_z <- as.numeric(scale(pm$linear_predictor))
pm_c <- concordance(Surv(OS_time, OS_status) ~ linear_predictor,
                    data = pm, reverse = TRUE)
pm_fit <- coxph(Surv(OS_time, OS_status) ~ lp_z, data = pm,
                 ties = "efron", x = TRUE)
pm_sm <- summary(pm_fit)$coefficients["lp_z", ]
pm_ci <- summary(pm_fit)$conf.int["lp_z", ]
pm_ph <- cox.zph(pm_fit, transform = "km")$table["lp_z", ]

summary_df <- data.frame(
  estimand = c(
    "contract_primary_patient_mean_ensemble_harrell_c",
    "supporting_fold_stratified_repeat_harrell_c_mean",
    "supporting_fold_stratified_repeat_harrell_c_sd",
    "supporting_repeat_cox_geometric_mean_hr_per_fold_sd",
    "supporting_repeat_cox_median_hr_per_fold_sd",
    "secondary_patient_mean_ensemble_cox_hr_per_sd"),
  estimate = c(
    unname(pm_c$concordance),
    mean(rep_df$fold_stratified_harrell_c),
    sd(rep_df$fold_stratified_harrell_c),
    exp(mean(rep_df$log_hr_per_fold_sd)),
    median(rep_df$hr_per_fold_sd),
    unname(pm_ci["exp(coef)"])),
  lower = c(
    NA, min(rep_df$fold_stratified_harrell_c), NA,
    min(rep_df$hr_per_fold_sd), min(rep_df$hr_per_fold_sd),
    unname(pm_ci["lower .95"])),
  upper = c(
    NA, max(rep_df$fold_stratified_harrell_c), NA,
    max(rep_df$hr_per_fold_sd), max(rep_df$hr_per_fold_sd),
    unname(pm_ci["upper .95"])),
  interval = c(
    "none; point estimate of repeated-CV patient-mean ensemble",
    "empirical range across 20 dependent repeats",
    "descriptive SD across 20 dependent repeats",
    "empirical HR range across 20 dependent repeats",
    "empirical HR range across 20 dependent repeats",
    "model-based 95% CI; descriptive ensemble association"),
  n = c(476, rep(476, 5)), events = c(220, rep(220, 5)),
  claim_role = c(
    "primary_by_original_sparse_contract_but_ensemble_not_single_model",
    "supporting_scale_valid_repeat_summary", "supporting_repeat_dispersion",
    "supporting_cross_fitted_association_no_combined_p",
    "supporting_cross_fitted_association_no_combined_p",
    "secondary_descriptive_not_confirmatory"),
  stringsAsFactors = FALSE)

ensemble_detail <- data.frame(
  n = nrow(pm), events = sum(pm$OS_status),
  harrell_c = unname(pm_c$concordance),
  hr_per_sd = unname(pm_ci["exp(coef)"]),
  ci_low_model_based = unname(pm_ci["lower .95"]),
  ci_high_model_based = unname(pm_ci["upper .95"]),
  nominal_p_descriptive = unname(pm_sm["Pr(>|z|)"]),
  ph_p = unname(pm_ph["p"]),
  stringsAsFactors = FALSE)

write.table(fold_df, file.path(out, "results/FOLD_LEVEL_CONCORDANCE_COUNTS.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(rep_df, file.path(out, "results/REPEAT_LEVEL_SCALE_VALID_C_AND_HR.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(summary_df, file.path(out, "results/OOF_ESTIMAND_SUMMARY.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
write.table(ensemble_detail, file.path(out, "results/PATIENT_MEAN_ENSEMBLE_DETAIL.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

input_paths <- c(pred_path, mean_path, task_path)
sha256_one <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (length(line) != 1L || !grepl("^[0-9a-f]{64}  ", line))
    stop("sha256sum failed for ", path)
  sub("  .*", "", line)
}
inputs <- data.frame(
  path = input_paths,
  md5 = unname(tools::md5sum(input_paths)),
  sha256 = vapply(input_paths, sha256_one, character(1)),
  stringsAsFactors = FALSE)
write.table(inputs, file.path(out, "results/INPUT_RECEIPT.tsv"), sep = "\t",
            quote = FALSE, row.names = FALSE)

capture.output(sessionInfo(), file = file.path(out, "SESSIONINFO.txt"))
cat("OOF estimand repair completed\n")
print(summary_df)
