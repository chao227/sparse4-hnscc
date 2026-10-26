#!/usr/bin/env Rscript

# ==========================================================
# Script   : 09_sparse4_construction.R
# Purpose  : Reproduce conditional elastic-net selection from eligible Yellow genes
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: sparse4 development; Fig. 4A-B; Supplementary Tables 1-2
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/05_sparse4/run_v31_historical_yellow_sparse_enet.R
# ==========================================================
suppressPackageStartupMessages({library(glmnet); library(survival); library(jsonlite)})
options(stringsAsFactors = FALSE)

need <- function(ok, ...) if (!isTRUE(ok)) stop(sprintf(...), call. = FALSE)
sha256sum <- function(path) {
  z <- system2("sha256sum", shQuote(path), stdout = TRUE, stderr = TRUE)
  need(length(z) >= 1L, "sha256sum failed for %s", path)
  sub("[[:space:]].*$", "", z[[1L]])
}
atomic_tsv <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  write.table(x, tmp, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
  need(file.rename(tmp, path), "cannot install %s", path)
}
atomic_json <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  jsonlite::write_json(x, tmp, auto_unbox = TRUE, pretty = TRUE, null = "null", digits = 16)
  need(file.rename(tmp, path), "cannot install %s", path)
}
atomic_rds <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(path, ".tmp-", Sys.getpid())
  saveRDS(x, tmp)
  need(file.rename(tmp, path), "cannot install %s", path)
}

args <- commandArgs(trailingOnly = TRUE)
need(length(args) >= 2L, "usage: full|task|aggregate CONTRACT [TASK_ID]")
mode <- args[[1L]]
contract_path <- normalizePath(args[[2L]], mustWork = TRUE)
need(mode %in% c("full", "task", "aggregate"), "unknown mode")
contract <- jsonlite::fromJSON(contract_path, simplifyVector = FALSE)
root <- normalizePath(getwd(), mustWork = TRUE)
need(!grepl("(^|/)reback2(/|$)", root), "reback2 is prohibited")
need(identical(contract$status, "FROZEN_BEFORE_SPARSE_MODEL_FIT"), "contract status drift")
for (nm in names(contract$inputs)) {
  p <- contract$inputs[[nm]]$path
  need(file.exists(p), "missing input %s", p)
  need(identical(sha256sum(p), contract$inputs[[nm]]$sha256), "input hash drift: %s", nm)
}

out <- file.path(root, contract$output_root)
workers <- file.path(root, "agent_work", "yellow_sparse_enet_nested_20261016", "workers")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(workers, recursive = TRUE, showWarnings = FALSE)

vst <- readRDS(contract$inputs$vst_matrix$path)
wgcna <- readRDS(contract$inputs$historical_wgcna$path)
surv <- read.csv(contract$inputs$survival$path, check.names = FALSE)
need(is.matrix(vst) && !anyDuplicated(rownames(vst)) && !anyDuplicated(colnames(vst)), "invalid locked VST matrix")
need(length(wgcna$gene_names) == length(wgcna$module_colors), "historical WGCNA membership alignment failure")
yellow_all <- sort(unique(as.character(wgcna$gene_names[tolower(as.character(wgcna$module_colors)) == "yellow"])))
yellow <- sort(intersect(yellow_all, rownames(vst)))
need(length(yellow_all) == contract$candidate_rule$historical_members_expected, "historical Yellow size drift")
need(length(yellow) == contract$candidate_rule$eligible_members_expected, "eligible Yellow size drift")
need(!anyDuplicated(surv$sample_id) && all(surv$sample_id %in% colnames(vst)), "survival/VST sample mapping failure")
need(nrow(surv) == contract$roster$patients && sum(surv$OS_status == 1L) == contract$roster$events, "roster drift")
need(all(is.finite(surv$OS_time) & surv$OS_time > 0) && all(surv$OS_status %in% 0:1), "invalid survival endpoint")
X <- t(vst[yellow, surv$sample_id, drop = FALSE])
Y <- Surv(surv$OS_time, surv$OS_status)
need(nrow(X) == contract$roster$patients && ncol(X) == contract$candidate_rule$eligible_members_expected && !anyNA(X) && all(is.finite(X)), "design matrix contract failure")

stratified_folds <- function(status, k, seed) {
  set.seed(seed)
  z <- integer(length(status))
  for (s in sort(unique(status))) {
    ix <- which(status == s)
    z[ix] <- sample(rep(seq_len(k), length.out = length(ix)))
  }
  need(all(tabulate(z, nbins = k) > 0L), "empty fold")
  z
}
task_split <- function(task) {
  m <- regexec("^R([0-9]{2})F([1-5])$", task)
  q <- regmatches(task, m)[[1L]]
  need(length(q) == 3L, "invalid task id %s", task)
  r <- as.integer(q[[2L]])
  k <- as.integer(q[[3L]])
  need(r >= 1L && r <= contract$internal_validation$outer_repeats, "repeat out of range")
  fold <- stratified_folds(surv$OS_status, contract$internal_validation$outer_folds,
                           contract$internal_validation$fold_seed + r)
  list(repeat_index = r, fold_index = k, train = which(fold != k), test = which(fold == k))
}
harrell_c <- function(time, event, score) {
  z <- survival::concordance(Surv(time, event) ~ score, reverse = TRUE)$concordance
  need(is.finite(z), "concordance unavailable")
  unname(z)
}
fit_sparse <- function(train, inner_seed) {
  inner <- stratified_folds(surv$OS_status[train], contract$elastic_net$inner_folds, inner_seed)
  alphas <- as.numeric(unlist(contract$elastic_net$alpha_grid))
  need(length(alphas) >= 1L && all(is.finite(alphas)) && all(alphas >= 0 & alphas <= 1), "invalid alpha grid")
  cvs <- vector("list", length(alphas))
  paths <- vector("list", length(alphas))
  for (i in seq_along(alphas)) {
    cvs[[i]] <- cv.glmnet(X[train, , drop = FALSE], Y[train], family = "cox", alpha = alphas[[i]],
                          foldid = inner, type.measure = "C", standardize = TRUE, parallel = FALSE)
    paths[[i]] <- data.frame(alpha = alphas[[i]], lambda_index = seq_along(cvs[[i]]$lambda),
                             lambda = cvs[[i]]$lambda, cvm = cvs[[i]]$cvm, cvsd = cvs[[i]]$cvsd,
                             cvup = cvs[[i]]$cvup, cvlo = cvs[[i]]$cvlo,
                             nzero = as.integer(cvs[[i]]$nzero), stringsAsFactors = FALSE)
  }
  surface <- do.call(rbind, paths)
  ok_best <- which(is.finite(surface$cvm) & is.finite(surface$cvsd))
  need(length(ok_best) > 0L, "no finite alpha/lambda candidate")
  best_cvm <- max(surface$cvm[ok_best])
  best_rows <- ok_best[surface$cvm[ok_best] == best_cvm]
  best_row <- best_rows[order(surface$cvsd[best_rows], -surface$lambda[best_rows], -surface$alpha[best_rows])][[1L]]
  threshold <- best_cvm - surface$cvsd[[best_row]]
  eligible <- which(is.finite(surface$cvm) & surface$cvm >= threshold & surface$nzero > 0L)
  need(length(eligible) > 0L, "no nonzero model within global one-SE envelope")
  ord <- order(surface$nzero[eligible], -surface$lambda[eligible], -surface$alpha[eligible], -surface$cvm[eligible])
  selected_row <- eligible[ord[[1L]]]
  surface$global_best <- FALSE
  surface$global_best[[best_row]] <- TRUE
  surface$within_global_one_se <- is.finite(surface$cvm) & surface$cvm >= threshold
  surface$selected <- FALSE
  surface$selected[[selected_row]] <- TRUE
  selected_alpha <- surface$alpha[[selected_row]]
  selected_lambda <- surface$lambda[[selected_row]]
  selected_cv <- cvs[[which(alphas == selected_alpha)[[1L]]]]
  beta <- as.matrix(coef(selected_cv, s = selected_lambda))[, 1L]
  beta <- beta[beta != 0]
  need(length(beta) == surface$nzero[[selected_row]] && length(beta) > 0L, "coefficient count mismatch")
  list(cv = selected_cv, all_cv = cvs, surface = surface, lambda = selected_lambda,
       alpha = selected_alpha, beta = beta, inner_fold = inner,
       best_cvm = best_cvm, threshold = threshold)
}

if (mode == "full") {
  candidate <- data.frame(historical_module = "Yellow", gene_symbol = yellow_all,
                          in_locked_vst = yellow_all %in% rownames(vst),
                          eligible_for_model = yellow_all %in% yellow,
                          exclusion_reason = ifelse(yellow_all %in% yellow, "", "ABSENT_FROM_LOCKED_VST_TOP5000"))
  atomic_tsv(candidate, file.path(out, "YELLOW_SPARSE_CANDIDATE_PROVENANCE.tsv"))
  fit <- fit_sparse(seq_len(nrow(X)), 86013L)
  pred <- as.numeric(predict(fit$cv, newx = X, s = fit$lambda, type = "link"))
  co <- data.frame(gene_symbol = names(fit$beta), coefficient = unname(fit$beta), stringsAsFactors = FALSE)
  sc <- data.frame(sample_id = surv$sample_id, patient_id = substr(surv$sample_id, 1L, 12L),
                   OS_time = surv$OS_time, OS_status = surv$OS_status, linear_predictor = pred)
  scale_tab <- data.frame(gene_symbol = yellow, training_mean = colMeans(X), training_sd = apply(X, 2L, sd),
                          selected = yellow %in% names(fit$beta), coefficient = unname(fit$beta[yellow]))
  scale_tab$coefficient[is.na(scale_tab$coefficient)] <- 0
  scale_tab$external_z_weight <- scale_tab$coefficient * scale_tab$training_sd
  atomic_tsv(co, file.path(out, "YELLOW_SPARSE_FULL_COEFFICIENTS.tsv"))
  atomic_tsv(sc, file.path(out, "YELLOW_SPARSE_FULL_SCORE.tsv"))
  atomic_tsv(scale_tab, file.path(out, "YELLOW_SPARSE_TRAINING_SCALE.tsv"))
  atomic_tsv(fit$surface, file.path(out, "YELLOW_SPARSE_FULL_SELECTION_SURFACE.tsv"))
  atomic_tsv(data.frame(sample_id = surv$sample_id, inner_fold = fit$inner_fold),
             file.path(out, "YELLOW_SPARSE_FULL_INNER_FOLDS.tsv"))
  atomic_rds(list(selected_cv = fit$cv, all_cv = fit$all_cv, lambda = fit$lambda, alpha = fit$alpha,
                  coefficients = fit$beta, candidate_genes = yellow, historical_yellow = yellow_all,
                  contract_sha256 = sha256sum(contract_path)), file.path(out, "YELLOW_SPARSE_FULL_MODEL.rds"))
  receipt <- list(schema = "V31_HISTORICAL_YELLOW_SPARSE_FULL_FIT_V1",
                  status = "FROZEN_BEFORE_EXTERNAL_EVALUATION_PENDING_INDEPENDENT_AUDIT",
                  created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
                  contract_sha256 = sha256sum(contract_path), patients = nrow(X), events = sum(surv$OS_status),
                  historical_yellow_n = length(yellow_all), candidate_n = length(yellow), selected_n = nrow(co),
                  alpha_grid = as.numeric(unlist(contract$elastic_net$alpha_grid)), selected_alpha = fit$alpha,
                  selected_lambda = fit$lambda, global_best_inner_c = fit$best_cvm,
                  global_one_se_threshold = fit$threshold,
                  apparent_harrell_c = harrell_c(surv$OS_time, surv$OS_status, pred),
                  claim_boundary = "DESCRIPTIVE_FULL_COHORT_FIT_NOT_INTERNAL_OR_EXTERNAL_VALIDATION",
                  prohibitions = c("NO_TARGET_GENE_COUNT", "NO_OUTCOME_GUIDED_REFIT", "NO_EXTERNAL_GUIDED_REFIT"))
  atomic_json(receipt, file.path(out, "YELLOW_SPARSE_FULL_FIT_RECEIPT.json"))
  cat(jsonlite::toJSON(receipt, auto_unbox = TRUE, pretty = TRUE), "\n")
  quit(save = "no", status = 0L)
}

if (mode == "task") {
  need(length(args) == 3L, "task mode requires TASK_ID")
  task <- args[[3L]]
  sp <- task_split(task)
  wout <- file.path(workers, task)
  dir.create(wout, recursive = TRUE, showWarnings = FALSE)
  fit <- fit_sparse(sp$train, 86013L + 100L * sp$repeat_index + sp$fold_index)
  pred <- as.numeric(predict(fit$cv, newx = X[sp$test, , drop = FALSE], s = fit$lambda, type = "link"))
  co <- data.frame(task_id = task, gene_symbol = names(fit$beta), coefficient = unname(fit$beta))
  pr <- data.frame(task_id = task, repeat_id = sprintf("R%02d", sp$repeat_index), outer_fold = sp$fold_index,
                   sample_id = surv$sample_id[sp$test], patient_id = substr(surv$sample_id[sp$test], 1L, 12L),
                   OS_time = surv$OS_time[sp$test], OS_status = surv$OS_status[sp$test], linear_predictor = pred)
  meta <- data.frame(task_id = task, repeat_id = sp$repeat_index, outer_fold = sp$fold_index,
                     train_n = length(sp$train), test_n = length(sp$test), train_events = sum(surv$OS_status[sp$train]),
                     test_events = sum(surv$OS_status[sp$test]), candidate_n = ncol(X), selected_n = length(fit$beta),
                     selected_alpha = fit$alpha, selected_lambda = fit$lambda,
                     global_best_inner_c = fit$best_cvm, global_one_se_threshold = fit$threshold,
                     heldout_harrell_c = harrell_c(surv$OS_time[sp$test], surv$OS_status[sp$test], pred))
  atomic_tsv(pr, file.path(wout, "HELDOUT_PREDICTIONS.tsv"))
  atomic_tsv(co, file.path(wout, "COEFFICIENTS.tsv"))
  atomic_tsv(meta, file.path(wout, "METADATA.tsv"))
  atomic_tsv(fit$surface, file.path(wout, "SELECTION_SURFACE.tsv"))
  receipt <- list(schema = "V31_HISTORICAL_YELLOW_SPARSE_NESTED_TASK_V1", status = "PASS", task_id = task,
                  contract_sha256 = sha256sum(contract_path), candidate_n = ncol(X), selected_n = length(fit$beta),
                  selected_alpha = fit$alpha, selected_lambda = fit$lambda, prediction_rows = nrow(pr),
                  artifact_sha256 = list(predictions = sha256sum(file.path(wout, "HELDOUT_PREDICTIONS.tsv")),
                    coefficients = sha256sum(file.path(wout, "COEFFICIENTS.tsv")),
                    metadata = sha256sum(file.path(wout, "METADATA.tsv")),
                    selection_surface = sha256sum(file.path(wout, "SELECTION_SURFACE.tsv"))))
  atomic_json(receipt, file.path(wout, "TASK_RECEIPT.json"))
  quit(save = "no", status = 0L)
}

task_ids <- unlist(lapply(seq_len(contract$internal_validation$outer_repeats),
                          function(r) sprintf("R%02dF%d", r, seq_len(contract$internal_validation$outer_folds))))
all_pr <- all_co <- all_meta <- vector("list", length(task_ids))
for (i in seq_along(task_ids)) {
  d <- file.path(workers, task_ids[[i]])
  rp <- file.path(d, "TASK_RECEIPT.json")
  need(file.exists(rp), "missing receipt %s", task_ids[[i]])
  rr <- jsonlite::fromJSON(rp, simplifyVector = FALSE)
  need(identical(rr$status, "PASS") && identical(rr$task_id, task_ids[[i]]), "non-PASS receipt %s", task_ids[[i]])
  files <- c(predictions = "HELDOUT_PREDICTIONS.tsv", coefficients = "COEFFICIENTS.tsv",
             metadata = "METADATA.tsv", selection_surface = "SELECTION_SURFACE.tsv")
  for (nm in names(files)) {
    need(identical(sha256sum(file.path(d, files[[nm]])), rr$artifact_sha256[[nm]]),
         "hash mismatch %s %s", task_ids[[i]], nm)
  }
  all_pr[[i]] <- read.delim(file.path(d, files[["predictions"]]), check.names = FALSE)
  all_co[[i]] <- read.delim(file.path(d, files[["coefficients"]]), check.names = FALSE)
  all_meta[[i]] <- read.delim(file.path(d, files[["metadata"]]), check.names = FALSE)
}
pr <- do.call(rbind, all_pr)
co <- do.call(rbind, all_co)
meta <- do.call(rbind, all_meta)
expected_predictions <- nrow(X) * contract$internal_validation$outer_repeats
need(nrow(pr) == expected_predictions && !anyDuplicated(paste(pr$task_id, pr$sample_id)), "prediction coverage failure")
repeat_c <- do.call(rbind, lapply(split(pr, pr$repeat_id), function(z) {
  data.frame(repeat_id = z$repeat_id[[1L]], n = nrow(z), events = sum(z$OS_status),
             heldout_harrell_c = harrell_c(z$OS_time, z$OS_status, z$linear_predictor))
}))
mean_pred <- aggregate(linear_predictor ~ sample_id + patient_id + OS_time + OS_status, pr, mean)
mean_c <- harrell_c(mean_pred$OS_time, mean_pred$OS_status, mean_pred$linear_predictor)
freq_count <- table(factor(co$gene_symbol, levels = yellow))
freq <- data.frame(gene_symbol = yellow, selected_outer_tasks = as.integer(freq_count),
                   selection_frequency = as.integer(freq_count) / length(task_ids))
freq <- freq[order(-freq$selection_frequency, freq$gene_symbol), ]
alpha_freq <- as.data.frame(table(meta$selected_alpha), stringsAsFactors = FALSE)
names(alpha_freq) <- c("selected_alpha", "outer_tasks")
alpha_freq$selection_frequency <- alpha_freq$outer_tasks / length(task_ids)
atomic_tsv(pr, file.path(out, "YELLOW_SPARSE_NESTED_HELDOUT_PREDICTIONS.tsv"))
atomic_tsv(mean_pred, file.path(out, "YELLOW_SPARSE_NESTED_PATIENT_MEAN_PREDICTIONS.tsv"))
atomic_tsv(repeat_c, file.path(out, "YELLOW_SPARSE_NESTED_REPEAT_CONCORDANCE.tsv"))
atomic_tsv(co, file.path(out, "YELLOW_SPARSE_NESTED_COEFFICIENTS.tsv"))
atomic_tsv(meta, file.path(out, "YELLOW_SPARSE_NESTED_TASK_METADATA.tsv"))
atomic_tsv(freq, file.path(out, "YELLOW_SPARSE_NESTED_SELECTION_FREQUENCY.tsv"))
atomic_tsv(alpha_freq, file.path(out, "YELLOW_SPARSE_NESTED_ALPHA_FREQUENCY.tsv"))
summary <- list(schema = "V31_HISTORICAL_YELLOW_SPARSE_NESTED_SUMMARY_V1", status = "PENDING_INDEPENDENT_AUDIT",
                created_utc = format(Sys.time(), tz = "UTC", usetz = TRUE), contract_sha256 = sha256sum(contract_path),
                outer_tasks = length(task_ids), outer_repeats = contract$internal_validation$outer_repeats,
                outer_folds = contract$internal_validation$outer_folds, patients = nrow(mean_pred),
                events = sum(mean_pred$OS_status), candidate_n = ncol(X), heldout_predictions = nrow(pr),
                patient_mean_heldout_harrell_c = mean_c,
                repeat_mean_harrell_c = mean(repeat_c$heldout_harrell_c),
                repeat_sd_harrell_c = sd(repeat_c$heldout_harrell_c),
                repeat_min_harrell_c = min(repeat_c$heldout_harrell_c),
                repeat_max_harrell_c = max(repeat_c$heldout_harrell_c),
                selected_n_median = median(meta$selected_n), selected_n_min = min(meta$selected_n),
                selected_n_max = max(meta$selected_n), claim_boundary = contract$internal_validation$claim_boundary)
atomic_json(summary, file.path(out, "YELLOW_SPARSE_NESTED_SUMMARY.json"))
cat(jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE), "\n")
