#!/usr/bin/env Rscript

# ==========================================================
# Script   : 12_objective64_internal_validation.R
# Purpose  : Reproduce the separate lambda-min procedure and its procedure-matched permutation evidence
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: secondary model procedure; Fig. 4H; Supplementary Tables 8-9
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:objective64_evidence_package_20261020/scripts/run_primary_nestedcv_and_permutation.R
# ==========================================================
# Frozen-protocol implementation. Primary fitting never loads legacy score resources.
suppressPackageStartupMessages({ library(survival); library(glmnet) })

base_dir <- normalizePath(getwd(), mustWork = TRUE)
if (!file.exists(file.path(base_dir, "protocol/OBJECTIVE_BIES_REBUILD_PROTOCOL_v1.sha256"))) stop("Frozen protocol SHA256 is absent; refusing real analysis.")
root_dir <- normalizePath(file.path(base_dir, ".."), mustWork = TRUE)
dir.create(file.path(base_dir, "results"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "logs"), showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(base_dir, "checkpoints"), showWarnings = FALSE, recursive = TRUE)

log_file <- file.path(base_dir, "logs/primary_nestedcv_permutation.log")
log_msg <- function(...) { msg <- paste0(...); line <- sprintf("[%s] %s", format(Sys.time(), "%F %T %z"), msg); cat(line, "\n"); cat(line, "\n", file = log_file, append = TRUE) }
worker_cores <- min(64L, max(1L, as.integer(Sys.getenv("OBJECTIVE_BIES_WORKERS", "64"))))

make_stratified_folds <- function(status, k, seed) {
  if (length(unique(status)) != 2L) stop("Both event strata are required for stratified folds.")
  set.seed(seed)
  fold <- integer(length(status))
  for (value in sort(unique(status))) {
    idx <- which(status == value)
    fold[idx] <- sample(rep(seq_len(k), length.out = length(idx)))
  }
  fold
}

cindex_counts <- function(time, status, lp) {
  fit <- survival::concordance(survival::Surv(time, status) ~ lp, reverse = TRUE)
  cnt <- fit$count
  get_count <- function(nm) if (nm %in% names(cnt)) as.numeric(cnt[[nm]]) else 0
  concordant <- get_count("concordant")
  discordant <- get_count("discordant")
  tied_x <- get_count("tied.x")
  tied_y <- get_count("tied.y")
  tied_xy <- get_count("tied.xy")
  comparable <- concordant + discordant + tied_x
  if (comparable <= 0) stop("No comparable survival pairs in held-out fold.")
  data.frame(
    test_cindex = (concordant + 0.5 * tied_x) / comparable,
    concordant_pairs = concordant,
    discordant_pairs = discordant,
    tied_prediction_pairs = tied_x,
    tied_time_pairs = tied_y,
    tied_both_pairs = tied_xy,
    comparable_pairs = comparable,
    stringsAsFactors = FALSE
  )
}

aggregate_cindex <- function(tab) {
  numerator <- sum(tab$concordant_pairs) + 0.5 * sum(tab$tied_prediction_pairs)
  denominator <- sum(tab$comparable_pairs)
  if (denominator <= 0) stop("No comparable pairs after within-fold aggregation.")
  numerator / denominator
}

toy_direction_check <- function() {
  # Higher LP means greater mortality risk; reverse=TRUE must return a high C-index.
  time <- c(10, 9, 8, 7, 6, 5, 4, 3)
  status <- rep(1L, length(time))
  lp <- seq_along(time)
  value <- cindex_counts(time, status, lp)$test_cindex
  if (!is.finite(value) || value <= 0.50) stop("Toy concordance direction check failed.")
  value
}

fit_outer_model <- function(X, time, status, train_idx, test_idx, inner_seed) {
  inner_fold <- make_stratified_folds(status[train_idx], k = 10L, seed = inner_seed)
  cvfit <- glmnet::cv.glmnet(
    x = X[train_idx, , drop = FALSE],
    y = survival::Surv(time[train_idx], status[train_idx]),
    family = "cox", alpha = 0.5, nfolds = 10L, foldid = inner_fold,
    type.measure = "C", standardize = TRUE
  )
  beta <- as.matrix(stats::coef(cvfit, s = "lambda.min"))
  selected <- rownames(beta)[beta[, 1] != 0]
  lp <- as.numeric(stats::predict(cvfit, newx = X[test_idx, , drop = FALSE], s = "lambda.min", type = "link"))
  list(
    inner_fold = inner_fold,
    lambda_min = as.numeric(cvfit$lambda.min),
    beta = beta,
    selected = selected,
    lp = lp,
    cindex = cindex_counts(time[test_idx], status[test_idx], lp)
  )
}

run_nested_cv <- function(X, surv, repeats = 20L, outer_k = 10L) {
  time <- surv$OS_time
  status <- surv$OS_status
  outer_assignments <- lapply(seq_len(repeats), function(repeat_id) {
    outer_fold <- make_stratified_folds(status, k = outer_k, seed = 1000L + repeat_id)
    data.frame(patient_id = rownames(X), repeat_id = repeat_id, outer_fold = outer_fold)
  })
  tasks <- expand.grid(repeat_id = seq_len(repeats), fold_id = seq_len(outer_k), KEEP.OUT.ATTRS = FALSE)
  tasks <- tasks[order(tasks$repeat_id, tasks$fold_id), , drop = FALSE]
  task_results <- parallel::mclapply(seq_len(nrow(tasks)), function(i) {
    repeat_id <- tasks$repeat_id[i]
    fold_id <- tasks$fold_id[i]
    outer_fold <- outer_assignments[[repeat_id]]$outer_fold
    test_idx <- which(outer_fold == fold_id)
    train_idx <- which(outer_fold != fold_id)
    inner_seed <- 100000L + repeat_id * 100L + fold_id
    fit <- fit_outer_model(X, time, status, train_idx, test_idx, inner_seed)
    list(
      fold = cbind(data.frame(repeat_id = repeat_id, fold = fold_id, train_n = length(train_idx), test_n = length(test_idx), train_events = sum(status[train_idx]), test_events = sum(status[test_idx]), lambda_min = fit$lambda_min, n_nonzero = length(fit$selected)), fit$cindex),
      selected = if (length(fit$selected)) data.frame(gene = fit$selected, coefficient = as.numeric(fit$beta[fit$selected, 1]), repeat_id = repeat_id, fold = fold_id) else NULL,
      oof = data.frame(patient_id = rownames(X)[test_idx], repeat_id = repeat_id, outer_fold = fold_id, oof_linear_predictor = fit$lp, time = time[test_idx], event = status[test_idx]),
      inner = data.frame(patient_id = rownames(X)[train_idx], repeat_id = repeat_id, outer_fold = fold_id, inner_fold = fit$inner_fold, inner_seed = inner_seed)
    )
  }, mc.cores = worker_cores, mc.preschedule = FALSE)
  fold_rows <- lapply(task_results, `[[`, "fold")
  selected_rows <- Filter(Negate(is.null), lapply(task_results, `[[`, "selected"))
  list(
    folds = do.call(rbind, fold_rows),
    selected = if (length(selected_rows)) do.call(rbind, selected_rows) else data.frame(gene = character(), coefficient = numeric(), repeat_id = integer(), fold = integer()),
    oof = do.call(rbind, lapply(task_results, `[[`, "oof")),
    outer_assignments = do.call(rbind, outer_assignments),
    inner_assignments = do.call(rbind, lapply(task_results, `[[`, "inner"))
  )
}

run_fixed_fivefold_stat <- function(X, surv, outer_fold, inner_base_seed, permutation_id, save_inner = TRUE) {
  time <- surv$OS_time
  status <- surv$OS_status
  fold_rows <- vector("list", 5L)
  inner_rows <- vector("list", 5L)
  for (fold_id in seq_len(5L)) {
    test_idx <- which(outer_fold == fold_id)
    train_idx <- which(outer_fold != fold_id)
    # Explicit foldid is produced from a deterministic seed for every inner call.
    inner_seed <- inner_base_seed + fold_id
    fit <- fit_outer_model(X, time, status, train_idx, test_idx, inner_seed)
    fold_rows[[fold_id]] <- cbind(data.frame(permutation_id = permutation_id, fold = fold_id, train_n = length(train_idx), test_n = length(test_idx), train_events = sum(status[train_idx]), test_events = sum(status[test_idx]), lambda_min = fit$lambda_min, n_nonzero = length(fit$selected), inner_seed = inner_seed), fit$cindex)
    if (save_inner) inner_rows[[fold_id]] <- data.frame(permutation_id = permutation_id, patient_id = rownames(X)[train_idx], outer_fold = fold_id, inner_fold = fit$inner_fold, inner_seed = inner_seed)
  }
  folds <- do.call(rbind, fold_rows)
  list(statistic = aggregate_cindex(folds), folds = folds, inner = if (save_inner) do.call(rbind, inner_rows) else NULL)
}

write_tsv <- function(x, path, append = FALSE, col.names = !append) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = col.names, append = append)

args <- commandArgs(trailingOnly = TRUE)
toy_value <- toy_direction_check()
log_msg("Toy risk-direction check C-index = ", sprintf("%.3f", toy_value))
if (identical(args, "--toy-only")) quit(save = "no", status = 0L)

feature <- read.delim(file.path(base_dir, "results/OBJECTIVE_FEATURE_UNIVERSE.tsv"), check.names = FALSE)
feature <- feature[, c("gene", "primary_universe")]
genes <- feature$gene[feature$primary_universe]
if (length(genes) != 492L || anyDuplicated(genes)) stop("Feature-universe integrity failure.")
vst <- readRDS(file.path(root_dir, "results/02_WGCNA/vst_top5000.rds"))
surv <- read.csv(file.path(root_dir, "results/04_ML/features/survival_data.csv"), check.names = FALSE, stringsAsFactors = FALSE)
if (anyDuplicated(surv$sample_id) || nrow(surv) != 476L) stop("Survival-ID integrity failure.")
surv <- surv[match(intersect(surv$sample_id, colnames(vst)), surv$sample_id), , drop = FALSE]
if (nrow(surv) != 476L || !all(surv$sample_id %in% colnames(vst))) stop("Locked sample match failure.")
X <- t(vst[genes, surv$sample_id, drop = FALSE])
rownames(X) <- surv$sample_id
if (!all(is.finite(X)) || any(apply(X, 2, var) <= 0)) stop("Primary feature technical eligibility failure.")
if (!all(is.finite(surv$OS_time) & surv$OS_time > 0) || anyNA(surv$OS_status)) stop("Survival-value integrity failure.")
surv$OS_status <- as.integer(surv$OS_status)
if (sum(surv$OS_status) != 220L) stop("Unexpected event count.")

writeLines(capture.output(sessionInfo()), file.path(base_dir, "logs/sessionInfo_primary_permutation.txt"))
log_msg("Loaded locked real data: ", nrow(X), " patients, ", ncol(X), " genes, ", sum(surv$OS_status), " events; deterministic workers = ", worker_cores, ".")

primary_fold_path <- file.path(base_dir, "results/NESTED_CV_FOLD_RESULTS.tsv")
primary_repeat_path <- file.path(base_dir, "results/NESTED_CV_REPEAT_LEVEL_RESULTS.tsv")
if (file.exists(primary_fold_path) && file.exists(primary_repeat_path)) {
  repeat_stat <- read.delim(primary_repeat_path, check.names = FALSE)
  if (nrow(read.delim(primary_fold_path, check.names = FALSE)) != 200L || nrow(repeat_stat) != 20L) stop("Primary checkpoint integrity failure.")
  log_msg("Resuming from validated primary nested-CV checkpoint.")
} else {
  primary <- run_nested_cv(X, surv, repeats = 20L, outer_k = 10L)
  primary_folds_out <- primary$folds; names(primary_folds_out)[names(primary_folds_out) == "repeat_id"] <- "repeat"
  primary_selected_out <- primary$selected; names(primary_selected_out)[names(primary_selected_out) == "repeat_id"] <- "repeat"
  primary_outer_out <- primary$outer_assignments; names(primary_outer_out)[names(primary_outer_out) == "repeat_id"] <- "repeat"
  primary_inner_out <- primary$inner_assignments; names(primary_inner_out)[names(primary_inner_out) == "repeat_id"] <- "repeat"
  primary_oof_out <- primary$oof; names(primary_oof_out)[names(primary_oof_out) == "repeat_id"] <- "repeat"
  write_tsv(primary_folds_out, primary_fold_path)
  write_tsv(primary_selected_out, file.path(base_dir, "results/NESTED_CV_SELECTED_GENES.tsv"))
  write_tsv(primary_oof_out, file.path(base_dir, "results/OOF_PREDICTIONS.tsv"))
  write_tsv(primary_outer_out, file.path(base_dir, "results/PRIMARY_OUTER_FOLD_ASSIGNMENTS.tsv"))
  write_tsv(primary_inner_out, file.path(base_dir, "results/PRIMARY_INNER_FOLD_ASSIGNMENTS.tsv"))
  repeat_stat <- do.call(rbind, lapply(split(primary$folds, primary$folds$repeat_id), function(z) data.frame(repeat_id = unique(z$repeat_id), repeat_cindex = aggregate_cindex(z), concordant_pairs = sum(z$concordant_pairs), discordant_pairs = sum(z$discordant_pairs), tied_prediction_pairs = sum(z$tied_prediction_pairs), comparable_pairs = sum(z$comparable_pairs))))
  names(repeat_stat)[names(repeat_stat) == "repeat_id"] <- "repeat"
  write_tsv(repeat_stat, primary_repeat_path)
  log_msg("Primary nested CV complete; mean repeat-level C-index = ", sprintf("%.5f", mean(repeat_stat$repeat_cindex)))
}

perm_outer_path <- file.path(base_dir, "results/PERMUTATION_OUTER_FOLD_ASSIGNMENTS.tsv")
if (file.exists(perm_outer_path)) {
  perm_outer <- read.delim(perm_outer_path, check.names = FALSE)$outer_fold
  if (length(perm_outer) != nrow(X)) stop("Permutation outer-fold checkpoint integrity failure.")
} else {
  perm_outer <- make_stratified_folds(surv$OS_status, k = 5L, seed = 400000L)
  write_tsv(data.frame(patient_id = rownames(X), outer_fold = perm_outer, outer_seed = 400000L), perm_outer_path)
}
observed_path <- file.path(base_dir, "results/PERMUTATION_OBSERVED_STATISTIC.tsv")
if (file.exists(observed_path)) {
  observed_statistic <- read.delim(observed_path, check.names = FALSE)$mean_outer_test_cindex[1]
  log_msg("Resuming from validated observed permutation statistic checkpoint.")
} else {
  observed <- run_fixed_fivefold_stat(X, surv, perm_outer, inner_base_seed = 400000L, permutation_id = 0L, save_inner = TRUE)
  write_tsv(observed$folds, file.path(base_dir, "results/PERMUTATION_OBSERVED_FOLD_RESULTS.tsv"))
  write_tsv(observed$inner, file.path(base_dir, "results/PERMUTATION_INNER_FOLD_ASSIGNMENTS.tsv"))
  write_tsv(data.frame(permutation_id = 0L, mean_outer_test_cindex = observed$statistic), observed_path)
  observed_statistic <- observed$statistic
  log_msg("Fixed-fivefold observed permutation statistic = ", sprintf("%.5f", observed_statistic))
}

null_path <- file.path(base_dir, "results/PERMUTATION_NULL.tsv")
progress_path <- file.path(base_dir, "checkpoints/permutation_progress.tsv")
if (file.exists(null_path)) {
  null <- read.delim(null_path, check.names = FALSE)
  done <- as.integer(null$permutation_id)
  if (anyDuplicated(done) || any(!done %in% seq_len(499L))) stop("Permutation resume file integrity failure.")
} else {
  done <- integer()
}
pending <- setdiff(seq_len(499L), done)
for (batch_start in seq(1L, length(pending), by = worker_cores)) {
  batch <- pending[batch_start:min(batch_start + worker_cores - 1L, length(pending))]
  batch_results <- parallel::mclapply(batch, function(b) {
    set.seed(500000L + b)
    perm_idx <- sample.int(nrow(X))
    X_perm <- X[perm_idx, , drop = FALSE]
    rownames(X_perm) <- rownames(X)
    one <- run_fixed_fivefold_stat(X_perm, surv, perm_outer, inner_base_seed = 500000L + b * 100L, permutation_id = b, save_inner = TRUE)
    list(permutation_id = b, statistic = one$statistic, inner = one$inner)
  }, mc.cores = min(worker_cores, length(batch)), mc.preschedule = FALSE)
  for (one in batch_results) {
    write_tsv(data.frame(permutation_id = one$permutation_id, mean_outer_test_cindex = one$statistic), null_path, append = file.exists(null_path))
    write_tsv(one$inner, file.path(base_dir, "results/PERMUTATION_INNER_FOLD_ASSIGNMENTS.tsv"), append = TRUE)
    write_tsv(data.frame(last_completed_permutation = one$permutation_id, timestamp = format(Sys.time(), "%FT%T%z")), progress_path)
  }
  log_msg("Permutation batch complete through ID ", max(batch), "/499; completed = ", sum(batch %in% batch), ".")
}
null <- read.delim(null_path, check.names = FALSE)
if (nrow(null) != 499L || !identical(sort(as.integer(null$permutation_id)), seq_len(499L))) stop("All 499 unique permutations are required.")
empirical_p <- (sum(null$mean_outer_test_cindex >= observed_statistic) + 1) / 500
summary <- data.frame(observed_statistic = observed_statistic, null_mean = mean(null$mean_outer_test_cindex), null_sd = sd(null$mean_outer_test_cindex), null_95_percentile = as.numeric(quantile(null$mean_outer_test_cindex, 0.95, names = FALSE)), null_99_percentile = as.numeric(quantile(null$mean_outer_test_cindex, 0.99, names = FALSE)), empirical_p = empirical_p, n_permutations = 499L)
write_tsv(summary, file.path(base_dir, "results/PERMUTATION_SUMMARY.tsv"))
write_tsv(data.frame(metric = c("n_patients", "n_events", "n_eligible_yellow_genes", "p_over_n", "p_over_events"), value = c(nrow(X), sum(surv$OS_status), ncol(X), ncol(X) / nrow(X), ncol(X) / sum(surv$OS_status))), file.path(base_dir, "results/FEATURE_UNIVERSE_DIMENSIONS.tsv"))
log_msg("Permutation complete; empirical P = ", sprintf("%.5f", empirical_p))
