#!/usr/bin/env Rscript

# ==========================================================
# Script   : 14_calibration_nonbootstrap.R
# Purpose  : Compute non-bootstrap apparent calibration summaries
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: calibration; Supplementary Fig. 2I-J
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/07_clinical/run_nonbootstrap_calibration.R
# ==========================================================

suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(jsonlite)
  library(digest)
})

options(stringsAsFactors = FALSE)

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop("Missing value after ", flag, call. = FALSE)
  args[[i + 1L]]
}
root <- normalizePath(arg_value("--root", "."), mustWork = TRUE)
module <- normalizePath(arg_value("--out", file.path(root, "calibration_nonbootstrap_20261020")), mustWork = TRUE)
stopifnot(basename(root) == "reback", startsWith(module, paste0(root, .Platform$file.sep)))

dirs <- file.path(module, c("results", "source_data", "figures", "logs", "qc", "reports"))
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

need <- function(ok, msg) if (!isTRUE(ok)) stop(msg, call. = FALSE)
sha256 <- function(path) digest::digest(file = path, algo = "sha256", serialize = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

paths <- c(
  clinical_master = file.path(root, "agent_work/sparse4_tcga_clinical_20261016/CLINICAL_MASTER_476.tsv"),
  independent_master = file.path(root, "agent_work/sparse4_tcga_clinical_independent_audit_20261016/INDEPENDENT_CLINICAL_MASTER_476.tsv"),
  independent_verdict = file.path(root, "agent_work/sparse4_tcga_clinical_independent_audit_20261016/INDEPENDENT_VERDICT.json"),
  model_lock = file.path(root, "results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/MODEL_LOCK.json"),
  score_table = file.path(root, "results/v3_1/P1_02_YELLOW_SPARSE4_MODEL_LOCK_20261016/SCORE_TABLE.tsv"),
  contract = file.path(module, "CALIBRATION_ANALYSIS_CONTRACT.md")
)
need(all(file.exists(paths)), "One or more frozen calibration inputs are missing")
need(sha256(paths[["clinical_master"]]) == sha256(paths[["independent_master"]]),
     "Production and independently recomputed clinical masters are not byte-identical")

verdict <- jsonlite::fromJSON(paths[["independent_verdict"]])
lock <- jsonlite::fromJSON(paths[["model_lock"]])
need(identical(verdict$status, "PASS"), "Current clinical source did not pass different-role audit")
need(identical(lock$status, "VALIDATED_MODEL_LOCK"), "Sparse4 MODEL_LOCK is not validated")
need(length(lock$genes_and_coefficients$gene_symbol) == 4L, "MODEL_LOCK is not sparse4")

master <- read.delim(paths[["clinical_master"]], check.names = FALSE, stringsAsFactors = FALSE)
score <- read.delim(paths[["score_table"]], check.names = FALSE, stringsAsFactors = FALSE)
need(nrow(master) == 476L && !anyDuplicated(master$patient_id) && sum(master$OS_status) == 220L,
     "Clinical master patient/event contract failed")
need(nrow(score) == 476L && !anyDuplicated(score$patient_id), "Locked score roster contract failed")
mi <- match(master$patient_id, score$patient_id)
need(!anyNA(mi), "Clinical master to locked score join is incomplete")
need(max(abs(master$score_raw_linear_predictor - score$linear_predictor[mi])) <= 1e-12,
     "Clinical score differs from locked sparse4 score")
need(max(abs(master$OS_time_months - score$OS_time[mi])) <= 1e-12 &&
       identical(as.integer(master$OS_status), as.integer(score$OS_status[mi])),
     "Clinical endpoint differs from locked score endpoint")

master$stage_factor <- factor(master$stage_factor, levels = c("I_II", "III_IV"))
master$HPV_p16_3level <- factor(master$HPV_p16_3level, levels = c("Negative", "Positive", "Unknown"))
master$CDE3108203_subsite <- factor(master$CDE3108203_subsite,
                                    levels = c("OralCavity", "Oropharynx", "Larynx_Hypopharynx"))
master$surgery_factor <- factor(master$surgery_factor, levels = c("No", "Yes"))

required <- c("OS_time_months", "OS_status", "age_years", "stage_factor", "HPV_p16_3level",
              "CDE3108203_subsite", "surgery_factor", "score_per_global_sd")
keep <- complete.cases(master[, required, drop = FALSE])
d <- droplevels(master[keep, , drop = FALSE])
need(nrow(d) == 320L && sum(d$OS_status) == 106L, "Derived complete-case roster is not 320 patients/106 events")
need(!anyDuplicated(d$patient_id) && all(d$OS_time_months > 0) && all(d$OS_status %in% c(0L, 1L)),
     "Complete-case patient/endpoint contract failed")

ind <- read.delim(paths[["independent_master"]], check.names = FALSE, stringsAsFactors = FALSE)
ind_keep <- complete.cases(ind[, required, drop = FALSE])
need(identical(sort(d$patient_id), sort(ind$patient_id[ind_keep])),
     "Derived roster differs from the independently recomputed roster")
missingness_flow <- data.frame(
  source_n = nrow(master), retained_n = nrow(d), excluded_n = sum(!keep),
  source_events = sum(master$OS_status), retained_events = sum(d$OS_status),
  excluded_events = sum(master$OS_status[!keep]),
  rule = paste("complete cases for", paste(required, collapse = ";")),
  stringsAsFactors = FALSE
)

forms <- list(
  CB = Surv(OS_time_months, OS_status) ~ age_years + stage_factor + HPV_p16_3level + CDE3108203_subsite + surgery_factor,
  CB_SCORE = Surv(OS_time_months, OS_status) ~ age_years + stage_factor + HPV_p16_3level + CDE3108203_subsite + surgery_factor + score_per_global_sd
)
fits <- lapply(forms, function(f) coxph(f, data = d, ties = "efron", x = TRUE, y = TRUE, model = TRUE))
need(all(vapply(fits, function(z) is.null(z$fail) && all(is.finite(coef(z))), logical(1))), "Cox fit failed")

predict_risk <- function(fit, newdata, horizon) {
  bh <- basehaz(fit, centered = FALSE)
  H0 <- if (any(bh$time <= horizon)) max(bh$hazard[bh$time <= horizon]) else 0
  lp <- as.numeric(predict(fit, newdata = newdata, type = "lp", reference = "zero"))
  risk <- 1 - exp(-H0 * exp(lp))
  need(all(is.finite(risk)) && all(risk >= 0 & risk <= 1), "Predicted absolute risk is invalid")
  risk
}

km_at <- function(dat, horizon) {
  fit <- survfit(Surv(OS_time_months, OS_status) ~ 1, data = dat, conf.type = "log-log")
  sm <- summary(fit, times = horizon, extend = TRUE)
  need(length(sm$surv) == 1L && is.finite(sm$surv), "Kaplan-Meier estimate is unavailable")
  data.frame(
    observed_risk = 1 - sm$surv,
    observed_risk_lower = max(0, 1 - sm$upper),
    observed_risk_upper = min(1, 1 - sm$lower),
    n_risk_at_horizon = as.integer(sm$n.risk),
    events_by_horizon = sum(dat$OS_status == 1L & dat$OS_time_months <= horizon),
    stringsAsFactors = FALSE
  )
}

horizons <- c(36L, 60L)
full_gates <- do.call(rbind, lapply(horizons, function(h) data.frame(
  horizon_months = h,
  n = nrow(d),
  events_total = sum(d$OS_status),
  n_risk_at_horizon = sum(d$OS_time_months >= h),
  events_by_horizon = sum(d$OS_status == 1L & d$OS_time_months <= h),
  gate_rule = "n_risk_at_horizon >= 20 and events_by_horizon >= 20",
  status = if (sum(d$OS_time_months >= h) >= 20L && sum(d$OS_status == 1L & d$OS_time_months <= h) >= 20L) "PASS" else "NOT_FEASIBLE_WITH_EVIDENCE"
)))
need(all(full_gates$status == "PASS"), "At least one frozen calibration horizon failed feasibility")

prediction_rows <- list()
group_rows <- list()
summary_rows <- list()
k <- 0L
for (h in horizons) {
  full_km <- km_at(d, h)
  for (id in names(fits)) {
    risk <- predict_risk(fits[[id]], d, h)
    ord <- order(risk, d$patient_id)
    group <- integer(nrow(d))
    group[ord] <- rep(seq_len(5L), each = nrow(d) / 5L)
    need(identical(as.integer(table(group)), rep(64L, 5L)), "Calibration groups are not five equal groups")
    pred <- data.frame(patient_id = d$patient_id, horizon_months = h, model_id = id,
                       predicted_event_risk = risk, calibration_group = group,
                       OS_time_months = d$OS_time_months, OS_status = d$OS_status,
                       stringsAsFactors = FALSE)
    prediction_rows[[length(prediction_rows) + 1L]] <- pred
    one_groups <- do.call(rbind, lapply(seq_len(5L), function(g) {
      z <- pred[pred$calibration_group == g, , drop = FALSE]
      km <- km_at(z, h)
      data.frame(horizon_months = h, model_id = id, calibration_group = g,
                 n = nrow(z), mean_predicted_risk = mean(z$predicted_event_risk),
                 median_predicted_risk = median(z$predicted_event_risk), km,
                 absolute_error = abs(km$observed_risk - mean(z$predicted_event_risk)),
                 interval_method = "pointwise 95% log-log Greenwood CI for 1-KM",
                 evidence_layer = "apparent same-cohort grouped calibration",
                 stringsAsFactors = FALSE)
    }))
    need(all(one_groups$n == 64L) && all(one_groups$n_risk_at_horizon >= 1L),
         "A calibration group has no patient at risk at the horizon")
    group_rows[[length(group_rows) + 1L]] <- one_groups
    summary_rows[[length(summary_rows) + 1L]] <- data.frame(
      horizon_months = h, model_id = id, n = nrow(d), events_total = sum(d$OS_status),
      n_risk_at_horizon = full_km$n_risk_at_horizon, events_by_horizon = full_km$events_by_horizon,
      mean_predicted_risk = mean(risk), observed_risk = full_km$observed_risk,
      observed_risk_lower = full_km$observed_risk_lower, observed_risk_upper = full_km$observed_risk_upper,
      absolute_overall_error = abs(full_km$observed_risk - mean(risk)),
      weighted_mean_absolute_group_error = weighted.mean(one_groups$absolute_error, one_groups$n),
      maximum_absolute_group_error = max(one_groups$absolute_error),
      evaluation_layer = "apparent same-cohort calibration; no optimism correction",
      uncertainty_scope = "Greenwood uncertainty for observed risk only; fitted-prediction and model-development uncertainty excluded",
      stringsAsFactors = FALSE
    )
    k <- k + 1L
  }
}

predictions <- do.call(rbind, prediction_rows)
groups <- do.call(rbind, group_rows)
summaries <- do.call(rbind, summary_rows)
coef_rows <- do.call(rbind, lapply(names(fits), function(id) {
  s <- summary(fits[[id]])
  data.frame(model_id = id, term = rownames(s$coefficients), coefficient = s$coefficients[, "coef"],
             se = s$coefficients[, "se(coef)"], hazard_ratio = s$conf.int[, "exp(coef)"],
             ci_lower = s$conf.int[, "lower .95"], ci_upper = s$conf.int[, "upper .95"],
             p_value = s$coefficients[, "Pr(>|z|)"], row.names = NULL, stringsAsFactors = FALSE)
}))

write_tsv(d[, c("patient_id", "sample_id", "OS_time_months", "OS_status")], file.path(module, "source_data/CALIBRATION_COMMON_ROSTER_320.tsv"))
write_tsv(predictions, file.path(module, "source_data/CALIBRATION_PATIENT_PREDICTIONS.tsv"))
write_tsv(groups, file.path(module, "source_data/CALIBRATION_GROUPED_KM.tsv"))
write_tsv(summaries, file.path(module, "results/CALIBRATION_APPARENT_SUMMARY.tsv"))
write_tsv(coef_rows, file.path(module, "results/CALIBRATION_MODEL_COEFFICIENTS.tsv"))
write_tsv(full_gates, file.path(module, "results/CALIBRATION_HORIZON_GATES.tsv"))
write_tsv(missingness_flow, file.path(module, "results/CALIBRATION_COMPLETE_CASE_FLOW.tsv"))
write_tsv(data.frame(input_id = names(paths), absolute_path = unname(paths), sha256 = vapply(paths, sha256, character(1))),
          file.path(module, "results/INPUT_RECEIPTS.tsv"))

labels <- c(CB = "Clinical base", CB_SCORE = "Clinical base + sparse4")
colors <- c(CB = "#4C78A8", CB_SCORE = "#C95D52")
shapes <- c(CB = 16, CB_SCORE = 17)
axis_max <- 1

theme_set(theme_classic(base_size = 7.5, base_family = "Arial") +
            theme(axis.line = element_line(linewidth = 0.35), axis.ticks = element_line(linewidth = 0.35),
                  axis.title = element_text(size = 7.5), axis.text = element_text(size = 7.5),
                  plot.title = element_text(size = 8.2, face = "bold", hjust = 0.5),
                  legend.position = "bottom", legend.title = element_blank(),
                  legend.text = element_text(size = 7.5), legend.key.width = grid::unit(5, "mm"),
                  legend.spacing.x = grid::unit(1.5, "mm"),
                  plot.margin = margin(3, 3, 2, 3, unit = "mm")))

for (h in horizons) {
  z <- groups[groups$horizon_months == h, , drop = FALSE]
  z$model_label <- factor(labels[z$model_id], levels = unname(labels))
  p <- ggplot(z, aes(mean_predicted_risk, observed_risk, color = model_id, shape = model_id, group = model_id)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "#777777", linewidth = 0.45) +
    geom_errorbar(aes(ymin = observed_risk_lower, ymax = observed_risk_upper), width = 0.012, linewidth = 0.45) +
    geom_line(linewidth = 0.55) + geom_point(size = 1.8, stroke = 0.25) +
    scale_color_manual(values = colors, breaks = names(labels), labels = c("Clinical", "Clinical + sparse4")) +
    scale_shape_manual(values = shapes, breaks = names(labels), labels = c("Clinical", "Clinical + sparse4")) +
    scale_x_continuous(limits = c(0, axis_max), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, axis_max), breaks = seq(0, 1, by = 0.25), expand = c(0, 0)) +
    coord_equal(clip = "off") +
    labs(title = sprintf("%d-month calibration\n(n = 320)", h),
         x = "Mean predicted risk", y = "Observed risk (1-KM)")
  panel_id <- if (h == 36L) "SFig2I_sparse4_calibration_36m" else "SFig2J_sparse4_calibration_60m"
  stem <- file.path(module, "figures", panel_id)
  ggsave(paste0(stem, ".pdf"), p, width = 60, height = 60, units = "mm", device = cairo_pdf)
  ggsave(paste0(stem, ".png"), p, width = 60, height = 60, units = "mm", dpi = 600, device = "png", bg = "white")
}

writeLines(capture.output(sessionInfo()), file.path(module, "logs/sessionInfo.txt"), useBytes = TRUE)

contract_record <- data.frame(
  field = c("independent_unit", "horizons_months", "roster", "models", "risk_definition", "grouping",
            "observed_risk", "uncertainty", "multiplicity", "bootstrap", "time_dependent_auc", "claim_ceiling"),
  value = c("patient", "36;60", "derived complete-case n=320; events=106", "CB;CB_SCORE",
            "1-S0(t)^exp(LP) from same fitted Cox model", "five deterministic equal-frequency groups per model; n=64 each",
            "1-Kaplan-Meier at horizon", "pointwise log-log Greenwood 95% CI for observed group risk only",
            "no hypothesis tests; no multiplicity adjustment", "PROHIBITED_NOT_RUN", "NOT_COMPUTED_NOT_DISPLAYED",
            "apparent same-cohort descriptive calibration only; no optimism correction, external validation, transportability, clinical utility or formal improvement claim"),
  stringsAsFactors = FALSE
)
write_tsv(contract_record, file.path(module, "results/EXECUTED_ESTIMAND_CONTRACT.tsv"))

cat(jsonlite::toJSON(list(status = "PRODUCTION_COMPLETE_PENDING_DIFFERENT_ROLE_AUDIT",
                          n = nrow(d), events = sum(d$OS_status), horizons = horizons,
                          models = names(fits), bootstrap = FALSE, time_dependent_auc = FALSE,
                          figures_count = as.integer(length(horizons))), auto_unbox = TRUE, pretty = TRUE), "\n")
