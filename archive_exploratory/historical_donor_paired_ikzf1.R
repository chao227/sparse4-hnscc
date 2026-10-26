#!/usr/bin/env Rscript

# ==========================================================
# Script   : 19_donor_paired_ikzf1.R
# Purpose  : Reproduce the donor-paired IKZF1 expression analysis retained as historical uncertainty evidence
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: donor-paired expression; Fig. 7B
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:scripts/modules/F/run_f4_cd8_pseudobulk.R
# ==========================================================

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(ggplot2)
})

script_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
if (length(script_arg) != 1L) stop("Cannot resolve script path")
script_path <- normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
if (!file.exists(file.path(root, "config", "params.yml"))) stop("Repository marker config/params.yml is missing", call. = FALSE)

input <- file.path(root, "inputs", "singlecell", "subcluster", "cd8_subclusters.rds")
out_dir <- file.path(root, "results", "F_singlecell")
log_dir <- file.path(root, "logs", "F_singlecell")
figure_dir <- file.path(root, "figures", "F_singlecell")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv <- function(x, path) {
  write.table(x, file = path, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")
}

write_tsv_gz <- function(x, path) {
  connection <- gzfile(path, open = "wt")
  on.exit(close(connection))
  write.table(x, file = connection, sep = "\t", quote = FALSE, row.names = FALSE,
              col.names = TRUE, na = "")
}

genes <- c("GZMK", "GNLY", "NKG7", "PRF1", "GZMB", "IFNG")
minimum_cells <- 20L
minimum_paired_donors <- 5L
bootstrap_reps <- 10000L
bootstrap_seed <- 70002L

object <- readRDS(input)
if (!inherits(object, "Seurat")) stop("Input is not a Seurat object")
if (!"patient" %in% colnames(object@meta.data)) stop("Missing patient metadata")
if (length(unique(object$patient)) != 18L) stop("Expected 18 donors")
if (!all(genes %in% rownames(object))) stop("Missing one or more prespecified genes")

normalized <- GetAssayData(object, assay = "RNA", layer = "data")
counts <- GetAssayData(object, assay = "RNA", layer = "counts")
if (!"IKZF1" %in% rownames(normalized)) stop("IKZF1 absent")
ikzf1 <- as.numeric(normalized["IKZF1", ])
names(ikzf1) <- colnames(object)
threshold <- unname(quantile(ikzf1, probs = 2 / 3, names = FALSE, type = 7))
group <- rep("MID_EXCLUDED", length(ikzf1))
names(group) <- names(ikzf1)
group[ikzf1 == 0] <- "IKZF1_low_zero"
group[ikzf1 >= threshold & ikzf1 > 0] <- "IKZF1_high_top_third"

meta <- object@meta.data
meta$cell_id <- rownames(meta)
meta$IKZF1_group <- group[meta$cell_id]
meta$IKZF1_raw_count <- as.numeric(counts["IKZF1", meta$cell_id])
meta$IKZF1_log_normalized <- ikzf1[meta$cell_id]
assignment_columns <- c("cell_id", "patient")
for (optional_column in c("cell_type", "cd8_subtype")) {
  if (optional_column %in% colnames(meta)) assignment_columns <- c(assignment_columns, optional_column)
}
cell_assignment <- meta[, assignment_columns, drop = FALSE]
cell_assignment$IKZF1_raw_count <- meta$IKZF1_raw_count
cell_assignment$IKZF1_log_normalized <- meta$IKZF1_log_normalized
cell_assignment$IKZF1_group <- meta$IKZF1_group
cell_assignment$included_in_pseudobulk <- meta$IKZF1_group %in% c(
  "IKZF1_high_top_third", "IKZF1_low_zero"
)
cell_assignment$threshold_probability <- 2 / 3
cell_assignment$threshold_type <- 7L
cell_assignment$threshold_global_log_normalized <- threshold
write_tsv_gz(
  cell_assignment,
  file.path(out_dir, "IKZF1_CD8_CELL_GROUP_ASSIGNMENTS.tsv.gz")
)
gate_counts <- as.data.frame.matrix(xtabs(~ patient + IKZF1_group, data = meta))
for (needed in c("IKZF1_high_top_third", "IKZF1_low_zero", "MID_EXCLUDED")) {
  if (!needed %in% colnames(gate_counts)) gate_counts[[needed]] <- 0L
}
gate_counts$patient <- rownames(gate_counts)
gate_counts$passes_high <- gate_counts$IKZF1_high_top_third >= minimum_cells
gate_counts$passes_low <- gate_counts$IKZF1_low_zero >= minimum_cells
gate_counts$paired_eligible <- gate_counts$passes_high & gate_counts$passes_low
gate_counts$minimum_cells_per_group <- minimum_cells
gate_counts$high_threshold_global_log_normalized <- threshold
gate_counts <- gate_counts[, c(
  "patient", "IKZF1_high_top_third", "IKZF1_low_zero", "MID_EXCLUDED",
  "passes_high", "passes_low", "paired_eligible", "minimum_cells_per_group",
  "high_threshold_global_log_normalized"
)]
write_tsv(gate_counts, file.path(out_dir, "IKZF1_CD8_PSEUDOBULK_GATE.tsv"))

eligible <- gate_counts$patient[gate_counts$paired_eligible]
gate_summary <- data.frame(
  input = input,
  n_cd8_cells = ncol(object),
  n_donors_total = length(unique(meta$patient)),
  high_definition = "global IKZF1 normalized expression >= 66.7th percentile and >0",
  low_definition = "IKZF1 normalized expression == 0",
  mid_definition = "nonzero cells below the global 66.7th percentile; excluded",
  high_threshold = threshold,
  minimum_cells_per_group = minimum_cells,
  minimum_paired_donors = minimum_paired_donors,
  observed_paired_donors = length(eligible),
  eligible_donor_ids = paste(eligible, collapse = ";"),
  gate = ifelse(length(eligible) >= minimum_paired_donors, "PASS", "FAIL"),
  old_random_subsampling_reused = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(gate_summary, file.path(out_dir, "IKZF1_CD8_PSEUDOBULK_GATE_SUMMARY.tsv"))

if (length(eligible) < minimum_paired_donors) {
  result <- data.frame(
    status = "NOT_FEASIBLE_WITH_EVIDENCE",
    reason = sprintf("Only %d donors pass; protocol requires >=%d", length(eligible), minimum_paired_donors),
    evidence = file.path(out_dir, "IKZF1_CD8_PSEUDOBULK_GATE.tsv")
  )
  write_tsv(result, file.path(out_dir, "IKZF1_CD8_PSEUDOBULK.tsv"))
  cat("F4_GATE_FAIL\n")
  quit(status = 0L)
}

sample_defs <- do.call(rbind, lapply(eligible, function(donor) {
  data.frame(
    patient = donor,
    group = c("IKZF1_high_top_third", "IKZF1_low_zero"),
    sample_id = paste(donor, c("high", "low"), sep = "__"),
    stringsAsFactors = FALSE
  )
}))

aggregate_columns <- lapply(seq_len(nrow(sample_defs)), function(i) {
  cells <- meta$cell_id[
    meta$patient == sample_defs$patient[i] & meta$IKZF1_group == sample_defs$group[i]
  ]
  if (length(cells) < minimum_cells) stop("Internal gate inconsistency")
  Matrix::rowSums(counts[, cells, drop = FALSE])
})
pseudobulk_counts <- do.call(cbind, aggregate_columns)
rownames(pseudobulk_counts) <- rownames(counts)
colnames(pseudobulk_counts) <- sample_defs$sample_id
if (any(colSums(pseudobulk_counts) <= 0)) stop("Zero pseudobulk library")

dge <- DGEList(counts = pseudobulk_counts)
keep_for_expression_gate <- filterByExpr(dge, group = factor(sample_defs$group))
if (sum(keep_for_expression_gate) < 100L) stop("Insufficient expressed genes for TMM normalization")
dge <- calcNormFactors(dge, method = "TMM")
log_cpm <- cpm(dge, log = TRUE, prior.count = 2)[genes, , drop = FALSE]

sample_values <- do.call(rbind, lapply(seq_len(nrow(sample_defs)), function(i) {
  data.frame(
    sample_defs[i, ],
    library_size = sum(pseudobulk_counts[, i]),
    norm_factor = dge$samples$norm.factors[i],
    t(log_cpm[, i, drop = FALSE]),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
}))
write_tsv(sample_values, file.path(out_dir, "IKZF1_CD8_PSEUDOBULK_DONOR_LOGCPM.tsv"))

raw_target_counts <- data.frame(
  sample_defs,
  cells_in_group = vapply(seq_len(nrow(sample_defs)), function(i) {
    sum(meta$patient == sample_defs$patient[i] & meta$IKZF1_group == sample_defs$group[i])
  }, integer(1)),
  library_size = colSums(pseudobulk_counts),
  t(as.matrix(pseudobulk_counts[genes, , drop = FALSE])),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
write_tsv(
  raw_target_counts,
  file.path(out_dir, "IKZF1_CD8_TARGET_RAW_PSEUDOBULK_COUNTS.tsv")
)

set.seed(bootstrap_seed)
result_rows <- lapply(genes, function(gene) {
  high_ids <- paste0(eligible, "__high")
  low_ids <- paste0(eligible, "__low")
  high <- as.numeric(log_cpm[gene, high_ids])
  low <- as.numeric(log_cpm[gene, low_ids])
  difference <- high - low
  test <- suppressWarnings(wilcox.test(high, low, paired = TRUE, exact = FALSE,
                                       alternative = "two.sided"))
  boot_medians <- replicate(
    bootstrap_reps,
    median(sample(difference, length(difference), replace = TRUE))
  )
  data.frame(
    gene = gene,
    paired_donors = length(difference),
    median_high_logCPM = median(high),
    median_low_logCPM = median(low),
    median_paired_difference_high_minus_low = median(difference),
    paired_difference_bootstrap_CI_lower = unname(quantile(boot_medians, 0.025, type = 6)),
    paired_difference_bootstrap_CI_upper = unname(quantile(boot_medians, 0.975, type = 6)),
    wilcoxon_statistic = unname(test$statistic),
    p_value = test$p.value,
    stringsAsFactors = FALSE
  )
})
results <- do.call(rbind, result_rows)
results$p_adjust_BH <- p.adjust(results$p_value, method = "BH")
results$bootstrap_reps <- bootstrap_reps
results$bootstrap_seed <- bootstrap_seed
results$analysis_unit <- "paired donor"
results$normalization <- "raw-count pseudobulk; TMM factors; edgeR logCPM prior.count=2"
results$claim_boundary <- "Exploratory donor-paired association; IKZF1-defined groups are observational"
results$status <- "COMPUTED_PENDING_AUDIT"
write_tsv(results, file.path(out_dir, "IKZF1_CD8_PSEUDOBULK.tsv"))

figure_source <- results[, c(
  "gene", "paired_donors", "median_paired_difference_high_minus_low",
  "paired_difference_bootstrap_CI_lower", "paired_difference_bootstrap_CI_upper",
  "p_value", "p_adjust_BH", "bootstrap_reps", "bootstrap_seed"
)]
figure_source$BH_status <- ifelse(figure_source$p_adjust_BH < 0.05, "BH q < 0.05", "BH q >= 0.05")
figure_source$display_order <- seq_len(nrow(figure_source))
write_tsv(figure_source, file.path(out_dir, "IKZF1_CD8_FIG7B_SOURCE_DATA.tsv"))

plot_data <- figure_source
plot_data$gene <- factor(plot_data$gene, levels = rev(genes))
plot_data$q_label <- sprintf("q = %.3g", plot_data$p_adjust_BH)
plot_data$BH_status <- factor(plot_data$BH_status, levels = c("BH q < 0.05", "BH q >= 0.05"))

fig7b <- ggplot(
  plot_data,
  aes(
    x = median_paired_difference_high_minus_low,
    y = gene,
    colour = BH_status
  )
) +
  geom_vline(xintercept = 0, linewidth = 0.35, linetype = "dashed", colour = "#606060") +
  geom_errorbar(
    aes(
      xmin = paired_difference_bootstrap_CI_lower,
      xmax = paired_difference_bootstrap_CI_upper
    ),
    orientation = "y",
    width = 0,
    linewidth = 0.65
  ) +
  geom_point(size = 2.2) +
  geom_text(
    aes(x = 0.075, label = q_label),
    hjust = 0,
    colour = "#202020",
    family = "Arial",
    size = 2.50
  ) +
  scale_colour_manual(values = c("BH q < 0.05" = "#2F6F9F", "BH q >= 0.05" = "#9AA0A6")) +
  coord_cartesian(xlim = c(-0.44, 0.19), clip = "off") +
  labs(
    x = "Median paired difference (log2 CPM)\nIKZF1-detected/high - IKZF1-undetected/zero",
    y = NULL,
    colour = NULL,
    subtitle = "n = 18 paired donors; 95% paired-donor bootstrap CI"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line = element_line(linewidth = 0.35, colour = "black"),
    axis.ticks = element_line(linewidth = 0.35, colour = "black"),
    axis.text = element_text(size = 7, colour = "black"),
    axis.text.y = element_text(face = "italic"),
    axis.title.x = element_text(size = 7.2, margin = margin(t = 5)),
    plot.subtitle = element_text(size = 7.2, margin = margin(b = 6)),
    legend.position = "bottom",
    legend.text = element_text(size = 7.0),
    legend.key.width = grid::unit(7, "pt"),
    plot.margin = margin(6, 18, 6, 6)
  ) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

figure_base <- file.path(figure_dir, "Fig7B_IKZF1_CD8_DONOR_PSEUDOBULK")
figure_width_in <- 89 / 25.4
figure_height_in <- 82 / 25.4
svglite::svglite(
  paste0(figure_base, ".svg"),
  width = figure_width_in,
  height = figure_height_in
)
print(fig7b)
grDevices::dev.off()
grDevices::cairo_pdf(
  paste0(figure_base, ".pdf"),
  width = figure_width_in,
  height = figure_height_in,
  family = "Arial",
  onefile = TRUE
)
print(fig7b)
grDevices::dev.off()
ragg::agg_png(
  paste0(figure_base, ".png"),
  width = figure_width_in,
  height = figure_height_in,
  units = "in",
  res = 600,
  background = "white"
)
print(fig7b)
grDevices::dev.off()

figure_contract <- data.frame(
  field = c(
    "core_conclusion", "archetype", "analysis_unit", "estimand", "interval",
    "multiplicity", "claim_ceiling", "source_data", "final_size_mm", "exports"
  ),
  value = c(
    "Within the 18 observed donors, IKZF1-detected/high and IKZF1-undetected/zero CD8 pseudobulks differ for four of six prespecified cytotoxic genes after BH correction.",
    "quantitative effect-size plot",
    "paired donor",
    "median across donors of within-donor TMM-logCPM high-minus-zero contrasts",
    "95% percentile interval from 10000 paired-donor bootstrap resamples; seed 70002; quantile type 6",
    "one Benjamini-Hochberg family across six prespecified genes",
    "exploratory observational association; no causal, cell-intrinsic, protein-validation or external-validation claim",
    "results/F_singlecell/IKZF1_CD8_FIG7B_SOURCE_DATA.tsv",
    "89 x 82",
    "600-dpi PNG, Cairo PDF and editable SVG with Arial text"
  ),
  stringsAsFactors = FALSE
)
write_tsv(figure_contract, file.path(out_dir, "IKZF1_CD8_FIG7B_CONTRACT.tsv"))

capture.output(sessionInfo(), file = file.path(log_dir, "F4_sessionInfo.txt"))
cat(sprintf("F4_PASS paired_donors=%d\n", length(eligible)))
