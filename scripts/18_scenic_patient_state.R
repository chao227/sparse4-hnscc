#!/usr/bin/env Rscript

# ==========================================================
# Script   : 18_scenic_patient_state.R
# Purpose  : Fit the patient-by-state IKZF1 regulon mixed model
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: SCENIC patient-state analysis; Fig. 5A
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:scan_1019_remediation_20261019/scripts/08_f3_patient_state_independent_recompute.R
# ==========================================================

suppressPackageStartupMessages({
  library(lme4)
  library(lmerTest)
  library(emmeans)
})
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
out <- file.path(root, "results", "P3_07_F3_patient_state")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

source_path <- file.path(root, "inputs", "IKZF1_patient_aware_AUCell.tsv")
legacy_model_path <- file.path(root, "inputs", "IKZF1_patient_aware_model.tsv")
d <- read_tsv(source_path)
legacy <- read_tsv(legacy_model_path)
state_levels <- c("Naive/Stem-like", "Progenitor-Tex", "Intermediate-Tex", "Terminal-Tex", "Cytotoxic-Eff", "Memory")
d$patient <- factor(d$patient)
d$CD8_state <- factor(d$CD8_state, levels = state_levels)
stopifnot(nrow(d) == 108L, length(unique(d$patient)) == 18L,
          all(table(d$patient) == 6L), all(table(d$CD8_state) == 18L), !anyDuplicated(d[c("patient", "CD8_state")]))

full <- lmerTest::lmer(IKZF1_AUCell_median ~ CD8_state + (1 | patient), data = d, REML = FALSE)
null <- lmerTest::lmer(IKZF1_AUCell_median ~ 1 + (1 | patient), data = d, REML = FALSE)
lrt <- anova(null, full)
conv_msg <- full@optinfo$conv$lme4$messages
if (is.null(conv_msg)) conv_msg <- ""
lrt_out <- data.frame(n_patient_state = nrow(d), patients = nlevels(d$patient), states = nlevels(d$CD8_state),
                      null_logLik = as.numeric(logLik(null)), full_logLik = as.numeric(logLik(full)),
                      lr_chisq = lrt$Chisq[2], df = lrt$Df[2], p_value = lrt$`Pr(>Chisq)`[2],
                      legacy_p = legacy$p_value[legacy$term == "global_CD8_state_LRT"],
                      abs_p_difference = abs(lrt$`Pr(>Chisq)`[2] - legacy$p_value[legacy$term == "global_CD8_state_LRT"]),
                      singular = isSingular(full), convergence_messages = paste(conv_msg, collapse = ";"))
write_tsv(lrt_out, file.path(out, "P3_07_F3_GLOBAL_LRT.tsv"))

cf <- summary(full)$coefficients
coef_out <- data.frame(term = rownames(cf), estimate = cf[, "Estimate"], se = cf[, "Std. Error"],
                       df = cf[, "df"], t = cf[, "t value"], p_value = cf[, "Pr(>|t|)"], row.names = NULL)
write_tsv(coef_out, file.path(out, "P3_07_F3_FIXED_EFFECTS.tsv"))

emm <- emmeans(full, ~ CD8_state)
pw <- as.data.frame(pairs(emm, adjust = "tukey"))
write_tsv(pw, file.path(out, "P3_07_F3_TUKEY_PAIRWISE.tsv"))

sens_full <- lm(IKZF1_AUCell_median ~ patient + CD8_state, data = d)
sens_null <- lm(IKZF1_AUCell_median ~ patient, data = d)
sa <- anova(sens_null, sens_full)
sens <- data.frame(n = nrow(d), patients = nlevels(d$patient), states = nlevels(d$CD8_state),
                   numerator_df = sa$Df[2], denominator_df = df.residual(sens_full), F = sa$F[2], p_value = sa$`Pr(>F)`[2])
write_tsv(sens, file.path(out, "P3_07_F3_FIXED_PATIENT_SENSITIVITY.tsv"))

coverage <- aggregate(n_cells ~ CD8_state, data = d, FUN = function(x) paste0("min=", min(x), ";median=", median(x), ";max=", max(x), ";lt20=", sum(x < 20)))
write_tsv(coverage, file.path(out, "P3_07_F3_CELL_COUNT_COVERAGE.tsv"))

legacy_pair <- legacy[legacy$model_component == "patient_summary_emmeans_Tukey_pairwise", ]
legacy_pair$key <- gsub("[()]", "", legacy_pair$term)
pw$key <- gsub("[()]", "", gsub(" - ", " - ", pw$contrast))
pair_cmp <- merge(pw[, c("key", "estimate", "SE", "df", "p.value")],
                  legacy_pair[, c("key", "estimate", "SE", "df", "p_value")], by = "key", suffixes = c("_scan", "_legacy"), all = TRUE)
pair_cmp$max_numeric_difference <- apply(pair_cmp[, c("estimate_scan", "estimate_legacy", "SE_scan", "SE_legacy", "df_scan", "df_legacy", "p.value", "p_value")], 1,
  function(x) max(abs(c(x[1]-x[2], x[3]-x[4], x[5]-x[6], x[7]-x[8])), na.rm = TRUE))
write_tsv(pair_cmp, file.path(out, "P3_07_F3_LEGACY_PAIRWISE_RECONCILIATION.tsv"))

qc <- data.frame(check_id = c("COMPLETE_18x6", "GLOBAL_P_MATCH", "PAIRWISE_15", "PAIRWISE_MAX_DIFF", "FIXED_PATIENT_DF", "NO_BOOTSTRAP"),
                 observed = c(nrow(d), lrt_out$abs_p_difference, nrow(pw), max(pair_cmp$max_numeric_difference), sens$denominator_df, TRUE),
                 expected = c(108, "<=1e-10", 15, "<=1e-8", 85, TRUE),
                 pass = c(nrow(d) == 108, lrt_out$abs_p_difference <= 1e-10, nrow(pw) == 15,
                          max(pair_cmp$max_numeric_difference) <= 1e-8, sens$denominator_df == 85, TRUE))
write_tsv(qc, file.path(out, "P3_07_F3_QC.tsv")); stopifnot(all(qc$pass))
writeLines(capture.output(sessionInfo()), file.path(out, "SESSIONINFO.txt"))
cat("F3 patient-state independent recomputation completed\n")
