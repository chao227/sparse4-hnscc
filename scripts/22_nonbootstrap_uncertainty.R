#!/usr/bin/env Rscript

# ==========================================================
# Script   : 22_nonbootstrap_uncertainty.R
# Purpose  : Compute patient-level non-bootstrap uncertainty summaries
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: uncertainty; Supplementary Table 14
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/10_sensitivities/compute_nonbootstrap_uncertainty.R
# ==========================================================

suppressPackageStartupMessages(library(survival))

root <- normalizePath(getwd(), mustWork = TRUE)
out <- file.path(root, "final_audit_reconciliation_20261022/production/uncertainty")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

oof_path <- file.path(
  root,
  "results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016/YELLOW_SPARSE_NESTED_HELDOUT_PREDICTIONS.tsv"
)
spatial_path <- file.path(
  root,
  "yellow_module_spatial_20261019/production/source_data/YELLOW_SPATIAL_PATIENT_ASSOCIATIONS.tsv"
)
spatial_primary_path <- file.path(
  root,
  "yellow_module_spatial_20261019/production/source_data/YELLOW_SPATIAL_PRIMARY_INFERENCE.tsv"
)
contract_path <- file.path(
  root,
  "final_audit_reconciliation_20261022/ANALYSIS_CONTRACT.md"
)
stopifnot(all(file.exists(c(oof_path, spatial_path, spatial_primary_path, contract_path))))

oof <- read.delim(oof_path, check.names = FALSE, stringsAsFactors = FALSE)
spatial <- read.delim(spatial_path, check.names = FALSE, stringsAsFactors = FALSE)
spatial_primary <- read.delim(spatial_primary_path, check.names = FALSE, stringsAsFactors = FALSE)

required_oof <- c(
  "repeat_id", "outer_fold", "patient_id", "OS_time", "OS_status",
  "linear_predictor"
)
required_spatial <- c("patient_id", "yellow_TLS10_rho")
stopifnot(all(required_oof %in% names(oof)), all(required_spatial %in% names(spatial)))
stopifnot(nrow(oof) == 9520L)
stopifnot(length(unique(oof$repeat_id)) == 20L)
stopifnot(length(unique(oof$patient_id)) == 476L)
stopifnot(!anyDuplicated(oof[c("repeat_id", "patient_id")]))
stopifnot(all(table(oof$repeat_id) == 476L))
stopifnot(all(oof$OS_status %in% c(0, 1)))
stopifnot(!anyNA(oof[required_oof]))
stopifnot(all(is.finite(oof$OS_time)), all(is.finite(oof$linear_predictor)))
stopifnot(all(oof$OS_time > 0))

patient_outcomes <- unique(oof[c("patient_id", "OS_time", "OS_status")])
stopifnot(nrow(patient_outcomes) == 476L)
stopifnot(!anyDuplicated(patient_outcomes$patient_id))
stopifnot(sum(patient_outcomes$OS_status) == 220L)

fold_comparable_repeat_c <- function(dat) {
  by_repeat <- split(dat, dat$repeat_id)
  repeat_c <- vapply(by_repeat, function(z) {
    counts <- setNames(
      rep(0, 5),
      c("concordant", "discordant", "tied.x", "tied.y", "tied.xy")
    )
    by_fold <- split(z, z$outer_fold)
    for (q in by_fold) {
      cc <- concordance(
        Surv(OS_time, OS_status) ~ linear_predictor,
        data = q,
        reverse = TRUE
      )
      counts <- counts + cc$count
    }
    denom <- counts[["concordant"]] + counts[["discordant"]] + counts[["tied.x"]]
    stopifnot(is.finite(denom), denom > 0)
    unname((counts[["concordant"]] + 0.5 * counts[["tied.x"]]) / denom)
  }, numeric(1))
  stopifnot(length(repeat_c) == 20L, all(is.finite(repeat_c)))
  list(estimate = mean(repeat_c), repeat_c = repeat_c)
}

oof_full <- fold_comparable_repeat_c(oof)
expected_oof <- 0.612901390909983
stopifnot(abs(oof_full$estimate - expected_oof) < 1e-12)

patient_ids <- sort(unique(oof$patient_id))
oof_loo <- lapply(patient_ids, function(pid) {
  z <- oof[oof$patient_id != pid, , drop = FALSE]
  est <- fold_comparable_repeat_c(z)
  data.frame(
    omitted_patient_id = pid,
    retained_patients = length(unique(z$patient_id)),
    retained_events = sum(unique(z[c("patient_id", "OS_status")])$OS_status),
    mean_fold_comparable_repeat_harrell_c = est$estimate,
    stringsAsFactors = FALSE
  )
})
oof_loo <- do.call(rbind, oof_loo)
stopifnot(nrow(oof_loo) == 476L, !anyNA(oof_loo))

jackknife_summary <- function(theta_full, theta_delete, n_units) {
  stopifnot(length(theta_delete) == n_units, all(is.finite(theta_delete)))
  theta_bar <- mean(theta_delete)
  variance <- (n_units - 1) / n_units * sum((theta_delete - theta_bar)^2)
  se <- sqrt(variance)
  z <- qnorm(0.975)
  data.frame(
    estimate = theta_full,
    jackknife_delete_mean = theta_bar,
    jackknife_bias_estimate = (n_units - 1) * (theta_bar - theta_full),
    jackknife_se = se,
    ci_low_normal_95 = theta_full - z * se,
    ci_high_normal_95 = theta_full + z * se,
    n_independent_units = n_units,
    stringsAsFactors = FALSE
  )
}

oof_summary <- jackknife_summary(
  oof_full$estimate,
  oof_loo$mean_fold_comparable_repeat_harrell_c,
  length(patient_ids)
)
oof_summary$analysis <- "sparse4_fixed_OOF_fold_comparable_mean_repeat_Harrell_C"
oof_summary$events <- 220L
oof_summary$interval_method <- "delete-one-patient jackknife variance; two-sided normal 95% interval"
oof_summary$estimand <- paste(
  "Arithmetic mean of 20 repeat-level Harrell C estimates; comparable pairs",
  "are pooled only within outer folds"
)
oof_summary$claim_role <- "conditional fixed-OOF uncertainty; no refitting, retuning or full-pipeline claim"
oof_summary <- oof_summary[c(
  "analysis", "estimate", "jackknife_delete_mean", "jackknife_bias_estimate",
  "jackknife_se", "ci_low_normal_95", "ci_high_normal_95",
  "n_independent_units", "events", "interval_method", "estimand", "claim_role"
)]

stopifnot(nrow(spatial) == 10L)
stopifnot(!anyDuplicated(spatial$patient_id))
stopifnot(!anyNA(spatial$yellow_TLS10_rho))
stopifnot(all(is.finite(spatial$yellow_TLS10_rho)))
rho <- spatial$yellow_TLS10_rho
spatial_mean <- mean(rho)
stopifnot(abs(spatial_mean - spatial_primary$mean_patient_rho[1]) < 1e-12)
stopifnot(all(rho > 0))
stopifnot(spatial_primary$exact_signflip_two_sided_p[1] == 0.001953125)

n_spatial <- length(rho)
t_se <- sd(rho) / sqrt(n_spatial)
t_critical <- qt(0.975, df = n_spatial - 1)
t_low <- spatial_mean - t_critical * t_se
t_high <- spatial_mean + t_critical * t_se

spatial_loo <- lapply(seq_along(rho), function(i) {
  data.frame(
    omitted_patient_id = spatial$patient_id[i],
    retained_patients = n_spatial - 1L,
    mean_patient_rho = mean(rho[-i]),
    stringsAsFactors = FALSE
  )
})
spatial_loo <- do.call(rbind, spatial_loo)
spatial_jk <- jackknife_summary(
  spatial_mean,
  spatial_loo$mean_patient_rho,
  n_spatial
)

spatial_summary <- data.frame(
  analysis = "patient_level_Yellow_activity_vs_TLS10_Spearman_rho",
  estimate = spatial_mean,
  median = median(rho),
  minimum = min(rho),
  maximum = max(rho),
  sd_across_patients = sd(rho),
  t_standard_error = t_se,
  t_ci_low_95 = t_low,
  t_ci_high_95 = t_high,
  t_df = n_spatial - 1L,
  jackknife_delete_mean = spatial_jk$jackknife_delete_mean,
  jackknife_bias_estimate = spatial_jk$jackknife_bias_estimate,
  jackknife_se = spatial_jk$jackknife_se,
  jackknife_ci_low_normal_95 = spatial_jk$ci_low_normal_95,
  jackknife_ci_high_normal_95 = spatial_jk$ci_high_normal_95,
  n_independent_patients = n_spatial,
  positive_patients = sum(rho > 0),
  exact_signflip_two_sided_p = spatial_primary$exact_signflip_two_sided_p[1],
  primary_interval_method = "two-sided Student t 95% interval across 10 patient effects",
  sensitivity_interval_method = "delete-one-patient jackknife variance; two-sided normal 95% interval",
  claim_role = paste(
    "patient-level effect uncertainty; sign-flip remains registered directional test;",
    "no spot/section independence, histological-TLS or causal claim"
  ),
  stringsAsFactors = FALSE
)

sha256_one <- function(path) {
  line <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  stopifnot(length(line) == 1L, grepl("^[0-9a-f]{64}  ", line))
  sub("  .*", "", line)
}

input_paths <- c(oof_path, spatial_path, spatial_primary_path, contract_path)
receipt <- data.frame(
  role = c(
    "locked_sparse4_heldout_predictions",
    "locked_spatial_patient_effects",
    "registered_spatial_primary_inference",
    "frozen_analysis_contract"
  ),
  path = normalizePath(input_paths, mustWork = TRUE),
  bytes = file.info(input_paths)$size,
  sha256 = vapply(input_paths, sha256_one, character(1)),
  stringsAsFactors = FALSE
)

write.table(
  oof_loo,
  file.path(out, "SPARSE4_OOF_DELETE_ONE_PATIENT_ESTIMATES.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  oof_summary,
  file.path(out, "SPARSE4_OOF_JACKKNIFE_INTERVAL.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  spatial_loo,
  file.path(out, "SPATIAL_YELLOW_TLS10_DELETE_ONE_PATIENT_ESTIMATES.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  spatial_summary,
  file.path(out, "SPATIAL_YELLOW_TLS10_INTERVALS.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
write.table(
  receipt,
  file.path(out, "INPUT_RECEIPT.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE, na = "NA"
)
capture.output(sessionInfo(), file = file.path(out, "sessionInfo.txt"))

cat("Completed contract A-B non-bootstrap intervals.\n")
print(oof_summary)
print(spatial_summary)
