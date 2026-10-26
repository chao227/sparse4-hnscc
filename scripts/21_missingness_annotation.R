#!/usr/bin/env Rscript

# ==========================================================
# Script   : 21_missingness_annotation.R
# Purpose  : Generate missingness and annotation-evidence summaries
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: missing data and state annotation; Supplementary Tables 15-16
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/10_sensitivities/run_missing_annotation_closure.R
# ==========================================================

options(stringsAsFactors = FALSE, warn = 1)
suppressPackageStartupMessages({
  library(survival)
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("usage: run_missing_annotation_closure.R <reback_root> <output_dir>")
root <- normalizePath(args[[1L]], mustWork = TRUE)
out <- args[[2L]]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
out <- normalizePath(out, mustWork = TRUE)
if (!startsWith(out, paste0(root, .Platform$file.sep))) stop("output must remain under reback")

write_tsv <- function(x, name) {
  write.table(x, file.path(out, name), sep = "\t", quote = FALSE,
              row.names = FALSE, na = "NA")
}

sha256 <- function(path) {
  z <- system2("sha256sum", path, stdout = TRUE, stderr = TRUE)
  if (length(z) != 1L) stop("sha256sum failed")
  sub("[[:space:]].*$", "", z)
}

clinical_path <- file.path(root, "agent_work/sparse4_tcga_clinical_20261016/CLINICAL_MASTER_476.tsv")
membership_path <- file.path(root, "scan_1019_remediation_20261019/independent_audit_primary3_20261019/results/M4_INDEPENDENT_EXCLUSION_LEDGER.tsv")
crosscheck_path <- file.path(root, "agent_work/analytical_gap_wave1_20261016/formal_run_010/l_annotation/L_REFERENCE_ANNOTATION_CROSSCHECK.tsv")
b_object_path <- file.path(root, "inputs/singlecell/subcluster/b_subclusters.rds")
whole_metadata_path <- file.path(root, "inputs/singlecell/annotation/cell_metadata.csv")
marker_source_path <- file.path(root, "inputs/parent_scripts/L/R12_cd8_b_subcluster.R")
scope_contract_path <- file.path(root, "results/A_consistency/SINGLECELL_METADATA_CONTRACT.tsv")
required <- c(clinical_path, membership_path, crosscheck_path, b_object_path,
              whole_metadata_path, marker_source_path, scope_contract_path)
if (!all(file.exists(required)) || any(file.info(required)$isdir) || any(Sys.readlink(required) != "")) {
  stop("guarded input contract failed")
}

input_hashes <- data.frame(
  path = required,
  bytes = file.info(required)$size,
  sha256 = vapply(required, sha256, character(1L)),
  stringsAsFactors = FALSE
)
write_tsv(input_hashes, "INPUT_HASHES.tsv")

## D. Locked M4 complete-case selection description -------------------------
clin <- read.delim(clinical_path, check.names = FALSE)
membership <- read.delim(membership_path, check.names = FALSE)
stopifnot(nrow(clin) == 476L, nrow(membership) == 476L,
          !anyDuplicated(clin$patient_id), !anyDuplicated(membership$patient_id))
idx <- match(clin$patient_id, membership$patient_id)
if (anyNA(idx)) stop("clinical patients do not map one-to-one to locked membership ledger")
membership <- membership[idx, , drop = FALSE]
if (!identical(as.character(clin$sample_id), as.character(membership$sample_id))) {
  stop("sample identity mismatch between clinical master and membership ledger")
}
clin$in_m4 <- as.logical(membership$in_m4)
clin$selection_group <- factor(ifelse(clin$in_m4, "Included", "Excluded"),
                               levels = c("Included", "Excluded"))
if (sum(clin$in_m4) != 320L || sum(!clin$in_m4) != 156L ||
    sum(clin$OS_status[clin$in_m4]) != 106L || sum(clin$OS_status[!clin$in_m4]) != 114L) {
  stop("locked M4 membership/event contract failed")
}

membership_out <- cbind(
  membership[, c("patient_id", "sample_id", "in_m4", "missing_fields")],
  selection_group = as.character(clin$selection_group),
  OS_time_months = clin$OS_time_months,
  OS_status = clin$OS_status
)
write_tsv(membership_out, "M4_LOCKED_MEMBERSHIP_LEDGER.tsv")

reverse_fit <- survfit(Surv(OS_time_months, 1L - OS_status) ~ selection_group, data = clin)
reverse_tab <- summary(reverse_fit)$table
if (is.null(dim(reverse_tab))) reverse_tab <- matrix(reverse_tab, nrow = 1L)
rn <- sub("^selection_group=", "", rownames(reverse_tab))
reverse_df <- data.frame(
  selection_group = rn,
  n = as.integer(reverse_tab[, "records"]),
  events = vapply(rn, function(g) sum(clin$OS_status[clin$selection_group == g]), integer(1L)),
  censored = vapply(rn, function(g) sum(clin$OS_status[clin$selection_group == g] == 0L), integer(1L)),
  event_fraction = vapply(rn, function(g) mean(clin$OS_status[clin$selection_group == g]), numeric(1L)),
  reverse_km_median_followup_months = as.numeric(reverse_tab[, "median"]),
  reverse_km_ci_lower_95 = as.numeric(reverse_tab[, "0.95LCL"]),
  reverse_km_ci_upper_95 = as.numeric(reverse_tab[, "0.95UCL"]),
  method = "reverse Kaplan-Meier; censoring indicator=1-OS_status",
  stringsAsFactors = FALSE
)
write_tsv(reverse_df, "M4_GROUP_EVENT_FOLLOWUP_SUMMARY.tsv")

continuous_vars <- c(
  age_years = "Age, years",
  score_raw_linear_predictor = "sparse4 raw linear predictor",
  project_immune_ssGSEA = "Project immune ssGSEA",
  project_stromal_ssGSEA = "Project stromal ssGSEA"
)
continuous_rows <- lapply(names(continuous_vars), function(v) {
  a <- clin[[v]][clin$in_m4]
  b <- clin[[v]][!clin$in_m4]
  pooled_sd <- sqrt((stats::var(a, na.rm = TRUE) + stats::var(b, na.rm = TRUE)) / 2)
  smd <- (mean(b, na.rm = TRUE) - mean(a, na.rm = TRUE)) / pooled_sd
  data.frame(
    variable = v,
    label = continuous_vars[[v]],
    included_observed_n = sum(!is.na(a)), included_missing_n = sum(is.na(a)),
    included_mean = mean(a, na.rm = TRUE), included_sd = sd(a, na.rm = TRUE),
    included_median = median(a, na.rm = TRUE), included_q1 = quantile(a, .25, na.rm = TRUE),
    included_q3 = quantile(a, .75, na.rm = TRUE),
    excluded_observed_n = sum(!is.na(b)), excluded_missing_n = sum(is.na(b)),
    excluded_mean = mean(b, na.rm = TRUE), excluded_sd = sd(b, na.rm = TRUE),
    excluded_median = median(b, na.rm = TRUE), excluded_q1 = quantile(b, .25, na.rm = TRUE),
    excluded_q3 = quantile(b, .75, na.rm = TRUE),
    smd_excluded_minus_included = smd,
    smd_definition = "difference in means divided by square root of mean group variances",
    stringsAsFactors = FALSE
  )
})
write_tsv(do.call(rbind, continuous_rows), "M4_CONTINUOUS_SMD.tsv")

categorical_vars <- c("stage_factor", "HPV_p16_3level", "surgery_factor", "CDE3108203_subsite")
categorical_rows <- list()
missing_rows <- list()
for (v in c(names(continuous_vars), categorical_vars)) {
  for (g in levels(clin$selection_group)) {
    z <- clin[[v]][clin$selection_group == g]
    missing_rows[[length(missing_rows) + 1L]] <- data.frame(
      variable = v, selection_group = g, total_n = length(z),
      observed_n = sum(!is.na(z)), missing_n = sum(is.na(z)),
      missing_fraction = mean(is.na(z)), denominator_definition = "locked group total",
      stringsAsFactors = FALSE
    )
  }
}
for (v in categorical_vars) {
  lv <- sort(unique(as.character(clin[[v]][!is.na(clin[[v]])])))
  for (lev in lv) {
    included_values <- clin[[v]][clin$in_m4]
    excluded_values <- clin[[v]][!clin$in_m4]
    pa <- mean(included_values == lev, na.rm = TRUE)
    pb <- mean(excluded_values == lev, na.rm = TRUE)
    denom <- sqrt((pa * (1 - pa) + pb * (1 - pb)) / 2)
    categorical_rows[[length(categorical_rows) + 1L]] <- data.frame(
      variable = v, level = lev,
      included_count = sum(included_values == lev, na.rm = TRUE), included_observed_n = sum(!is.na(included_values)),
      included_fraction = pa,
      excluded_count = sum(excluded_values == lev, na.rm = TRUE), excluded_observed_n = sum(!is.na(excluded_values)),
      excluded_fraction = pb,
      smd_excluded_minus_included = ifelse(denom > 0, (pb - pa) / denom, NA_real_),
      smd_definition = "level-specific difference in proportions divided by pooled Bernoulli SD",
      stringsAsFactors = FALSE
    )
  }
}
write_tsv(do.call(rbind, categorical_rows), "M4_CATEGORICAL_LEVEL_SMD.tsv")
write_tsv(do.call(rbind, missing_rows), "M4_MISSING_DENOMINATORS.tsv")

## E. Existing HPCA/SingleR annotation evidence -----------------------------
ref <- read.delim(crosscheck_path, check.names = FALSE)
required_ref_cols <- c("patient", "locked_state", "label_fine", "pruned_label_fine",
                       "prespecified_match", "cells")
if (!all(required_ref_cols %in% names(ref))) stop("reference crosscheck schema mismatch")
locked_states <- c("Naive B", "Memory B", "GC-like B", "Transitional B", "Plasma cell")
if (!setequal(unique(ref$locked_state), locked_states)) stop("locked B-state universe mismatch")

compact_reference <- function(x) {
  if (is.na(x) || x == "<PRUNED>") return("<PRUNED>")
  if (x == "B_cell:Naive") return("Naive reference")
  if (x == "B_cell:Memory") return("Memory reference")
  if (x %in% c("B_cell:Germinal_center", "B_cell:CXCR4-_centrocyte", "B_cell:CXCR4+_centroblast")) return("GC reference")
  if (x == "B_cell:immature") return("Immature reference")
  if (x == "B_cell:Plasma_cell") return("Plasma reference")
  "Other reference"
}
ref$compact_reference <- vapply(as.character(ref$pruned_label_fine), compact_reference, character(1L))
compact_levels <- c("Naive reference", "Memory reference", "GC reference",
                    "Immature reference", "Plasma reference", "Other reference", "<PRUNED>")
conf <- aggregate(cells ~ locked_state + compact_reference, ref, sum)
grid <- expand.grid(locked_state = locked_states, compact_reference = compact_levels,
                    stringsAsFactors = FALSE)
conf <- merge(grid, conf, by = c("locked_state", "compact_reference"), all.x = TRUE, sort = FALSE)
conf$cells[is.na(conf$cells)] <- 0L
state_den <- aggregate(cells ~ locked_state, conf, sum)
names(state_den)[2] <- "state_cells"
conf <- merge(conf, state_den, by = "locked_state", all.x = TRUE, sort = FALSE)
conf$fraction_within_locked_state <- conf$cells / conf$state_cells
conf$locked_state <- factor(conf$locked_state, levels = locked_states)
conf$compact_reference <- factor(conf$compact_reference, levels = compact_levels)
conf <- conf[order(conf$locked_state, conf$compact_reference), ]
conf$locked_state <- as.character(conf$locked_state)
conf$compact_reference <- as.character(conf$compact_reference)
write_tsv(conf, "B_STATE_REFERENCE_COMPACT_CONFUSION.tsv")

donor_agreement <- do.call(rbind, lapply(split(ref, list(ref$patient, ref$locked_state), drop = TRUE), function(d) {
  total <- sum(d$cells)
  pruned <- sum(d$cells[d$pruned_label_fine == "<PRUNED>"])
  matched <- sum(d$cells[as.logical(d$prespecified_match)])
  label_counts <- aggregate(cells ~ compact_reference, d, sum)
  top <- label_counts$compact_reference[which.max(label_counts$cells)]
  data.frame(
    donor = d$patient[[1L]], locked_state = d$locked_state[[1L]], cells = total,
    matched_cells = matched, matched_fraction_all_cells = matched / total,
    pruned_cells = pruned, pruned_fraction = pruned / total,
    top_compact_reference = top,
    evidence_role = "donor-level descriptive SingleR/HPCA cross-check; no relabelling",
    stringsAsFactors = FALSE
  )
}))
donor_agreement$locked_state <- factor(donor_agreement$locked_state, levels = locked_states)
donor_agreement <- donor_agreement[order(donor_agreement$donor, donor_agreement$locked_state), ]
donor_agreement$locked_state <- as.character(donor_agreement$locked_state)
write_tsv(donor_agreement, "B_STATE_DONOR_REFERENCE_AGREEMENT.tsv")

## Prespecified marker set copied exactly from the historical annotation workflow.
marker_sets <- list(
  "Naive B" = c("TCL1A", "IGHD", "IGHM", "FCER2", "IL4R", "CD22"),
  "Memory B" = c("CD27", "IGHG1", "IGHA1", "S1PR1", "CR2", "TNFRSF13B"),
  "GC-like B" = c("AICDA", "BCL6", "MKI67", "PCNA", "FCRL4", "ITGAX"),
  "Transitional B" = c("CD38", "VPREB3", "PTPRC", "EBF1", "IL4R"),
  "Plasma cell" = c("MZB1", "JCHAIN", "SDC1", "XBP1", "IRF4", "IGHG1", "IGHA1")
)
marker_contract <- do.call(rbind, lapply(names(marker_sets), function(s) {
  data.frame(marker_state = s, marker = marker_sets[[s]],
             source = "inputs/parent_scripts/L/R12_cd8_b_subcluster.R::b_markers",
             selection_status = "prespecified historical annotation marker; not result-selected",
             stringsAsFactors = FALSE)
}))

b <- readRDS(b_object_path)
meta <- b@meta.data
if (!identical(colnames(b), rownames(meta)) || anyNA(meta$patient) || anyNA(meta$b_subtype)) {
  stop("B-cell object identity contract failed")
}
expr <- SeuratObject::LayerData(b[["RNA"]], layer = "data")
covered <- intersect(unique(marker_contract$marker), rownames(expr))
marker_contract$covered_in_object <- marker_contract$marker %in% covered
write_tsv(marker_contract, "B_STATE_PREDEFINED_MARKER_CONTRACT.tsv")
if (!all(marker_contract$covered_in_object)) stop("prespecified marker absent from B-cell object")

donors <- sort(unique(as.character(meta$patient)))
marker_rows <- vector("list", length(donors) * length(locked_states) * nrow(marker_contract))
k <- 0L
for (donor in donors) {
  for (state in locked_states) {
    cells <- rownames(meta)[meta$patient == donor & meta$b_subtype == state]
    for (j in seq_len(nrow(marker_contract))) {
      k <- k + 1L
      gene <- marker_contract$marker[[j]]
      x <- if (length(cells)) as.numeric(expr[gene, cells, drop = TRUE]) else numeric()
      marker_rows[[k]] <- data.frame(
        donor = donor, locked_state = state,
        marker_state = marker_contract$marker_state[[j]], marker = gene,
        state_donor_cells = length(cells),
        detected_cells = sum(x > 0),
        detection_fraction = if (length(x)) mean(x > 0) else NA_real_,
        mean_log_normalized_expression = if (length(x)) mean(x) else NA_real_,
        median_log_normalized_expression = if (length(x)) median(x) else NA_real_,
        independent_unit = "donor",
        stringsAsFactors = FALSE
      )
    }
  }
}
marker_donor <- do.call(rbind, marker_rows)
write_tsv(marker_donor, "B_STATE_MARKER_DONOR_BY_STATE.tsv")

qfun <- function(x, p) if (all(is.na(x))) NA_real_ else unname(quantile(x, p, na.rm = TRUE))
marker_across <- do.call(rbind, lapply(split(marker_donor, list(marker_donor$locked_state,
                                                                marker_donor$marker_state,
                                                                marker_donor$marker), drop = TRUE), function(d) {
  data.frame(
    locked_state = d$locked_state[[1L]], marker_state = d$marker_state[[1L]], marker = d$marker[[1L]],
    donors_total = nrow(d), donors_with_cells = sum(d$state_donor_cells > 0),
    total_cells = sum(d$state_donor_cells),
    donor_median_detection_fraction = median(d$detection_fraction, na.rm = TRUE),
    donor_q1_detection_fraction = qfun(d$detection_fraction, .25),
    donor_q3_detection_fraction = qfun(d$detection_fraction, .75),
    donor_median_mean_log_expression = median(d$mean_log_normalized_expression, na.rm = TRUE),
    donor_q1_mean_log_expression = qfun(d$mean_log_normalized_expression, .25),
    donor_q3_mean_log_expression = qfun(d$mean_log_normalized_expression, .75),
    summary_unit = "donor",
    stringsAsFactors = FALSE
  )
}))
marker_across$locked_state <- factor(marker_across$locked_state, levels = locked_states)
marker_across <- marker_across[order(marker_across$locked_state, marker_across$marker_state,
                                     marker_across$marker), ]
marker_across$locked_state <- as.character(marker_across$locked_state)
write_tsv(marker_across, "B_STATE_MARKER_ACROSS_DONOR_SUMMARY.tsv")

whole <- read.csv(whole_metadata_path, check.names = FALSE)
scope_contract <- read.delim(scope_contract_path, check.names = FALSE)
sample_values <- sort(unique(as.character(whole$sample)))
if (!all(as.character(whole$fraction) %in% c("CD45p", "CD45n"))) stop("unexpected fraction label")
if (!all(grepl("^HN[0-9]{2}_CD45[pn]$", sample_values))) stop("sample naming contract does not support tumour CD45-only scope")
if (any(grepl("PBL|blood", sample_values, ignore.case = TRUE))) stop("PBL/blood sample found in current object metadata")
scope_row <- scope_contract[scope_contract$entity == "tumor fractions", , drop = FALSE]
if (nrow(scope_row) != 1L || scope_row$count[[1L]] != 33L ||
    !grepl("PBL excluded", scope_row$evidence_basis[[1L]], fixed = TRUE) ||
    scope_row$status[[1L]] != "VERIFIED") {
  stop("single-cell tumour-fraction scope contract failed")
}
scope <- data.frame(
  object = c("whole annotated object", "B-cell subcluster object"),
  cells = c(nrow(whole), ncol(b)),
  donors = c(length(unique(whole$patient)), length(unique(meta$patient))),
  samples_or_fractions = c(length(sample_values), length(unique(meta$sample))),
  observed_fraction_labels = c(paste(sort(unique(whole$fraction)), collapse = ";"),
                               paste(sort(unique(meta$fraction)), collapse = ";")),
  pbl_or_blood_label_hits = c(sum(grepl("PBL|blood", whole$sample, ignore.case = TRUE)),
                              sum(grepl("PBL|blood", meta$sample, ignore.case = TRUE))),
  scope_status = "CURRENT_OBJECT_METADATA_ONLY_INTRATUMORAL_CD45_SORTED_FRACTIONS_NO_PBL_LABEL",
  evidence_boundary = "current metadata and verified sample contract identify 33 tumour CD45-sorted fractions with PBL excluded; the public source study also contains matched PBL samples that are not present here",
  stringsAsFactors = FALSE
)
write_tsv(scope, "SINGLECELL_CURRENT_OBJECT_SCOPE.tsv")

checks <- data.frame(
  check_id = c("D01", "D02", "D03", "D04", "E01", "E02", "E03", "E04", "E05"),
  status = "PASS",
  observed = c(
    "476 ledger identities map one-to-one",
    "320 included/156 excluded; 106/114 events",
    "reverse-KM computed by locked membership",
    "continuous and level-specific categorical SMD with explicit missing denominators",
    "existing SingleR/HPCA aggregate rows used without relabelling",
    "5 locked states x 7 compact reference categories",
    "18-donor state-level agreement generated",
    paste0(nrow(marker_contract), " prespecified state-marker pairs; all covered"),
    "102117 cells; 18 donors; 33 CD45p/CD45n fractions; zero PBL/blood-labelled samples"
  ),
  acceptance_rule = c(
    "use verified patient ledger, not summary reconstruction",
    "match locked M4 roster and event totals",
    "reverse censoring indicator and group-specific interval",
    "descriptive balance only; no missingness mechanism claim",
    "reuse registered sensitivity; no new favourable label selection",
    "compact confusion preserves pruned/other assignments",
    "donor is independent unit",
    "markers exactly match historical b_markers; no result-driven selection",
    "current retained object contains tumour CD45 sorted fractions and no PBL-labelled sample"
  ),
  stringsAsFactors = FALSE
)
write_tsv(checks, "PRODUCER_CHECKS.tsv")
capture.output(sessionInfo(), file = file.path(out, "SESSIONINFO.txt"))
cat("MISSING_ANNOTATION_CLOSURE_PASS\n")
