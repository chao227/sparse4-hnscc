#!/usr/bin/env Rscript

# ==========================================================
# Script   : 08_yellow_module_context.R
# Purpose  : Generate Yellow-module single-cell and pathway context results
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Results: Yellow-module biology; Fig. 2E-H
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/03_yellow_context/run_fig2_yellow_production.R
# ==========================================================

options(stringsAsFactors = FALSE, width = 180)
suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(edgeR)
  library(data.table)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(clusterProfiler)
  library(ReactomePA)
  library(jsonlite)
})

root <- normalizePath(".", mustWork = TRUE)
module_dir <- file.path(root, "replacement_fig2_yellow_panels_20261019")
out_dir <- file.path(module_dir, "production")
table_dir <- file.path(out_dir, "tables")
report_dir <- file.path(out_dir, "reports")
log_dir <- file.path(out_dir, "logs")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

inputs <- c(
  yellow = file.path(root, "inputs/wgcna/legacy/hub_genes_all.csv"),
  whole = file.path(root, "inputs/singlecell/annotation/seurat_annotated.rds"),
  bcell = file.path(root, "inputs/singlecell/subcluster/b_subclusters.rds"),
  vst = file.path(root, "inputs/wgcna/legacy/vst_top25mad.rds")
)
stopifnot(all(file.exists(inputs)))

write_tsv <- function(x, path) {
  fwrite(as.data.table(x), path, sep = "\t", na = "NA", quote = FALSE)
}

sha256_file <- function(path) {
  out <- system2("sha256sum", path, stdout = TRUE)
  strsplit(out, "[[:space:]]+")[[1]][1]
}

aggregate_counts <- function(counts, patient, label) {
  stopifnot(ncol(counts) == length(patient), length(patient) == length(label))
  grp <- factor(paste(patient, label, sep = "::"))
  incidence <- sparseMatrix(
    i = seq_along(grp), j = as.integer(grp), x = 1,
    dims = c(length(grp), nlevels(grp)),
    dimnames = list(colnames(counts), levels(grp))
  )
  aggregated <- counts %*% incidence
  colnames(aggregated) <- levels(grp)
  list(counts = aggregated, observed = as.data.frame(table(patient, label), stringsAsFactors = FALSE))
}

make_unit_ledger <- function(patient_levels, label_levels, observed, label_name, minimum) {
  full <- expand.grid(patient = patient_levels, label = label_levels, stringsAsFactors = FALSE)
  names(observed) <- c("patient", "label", "n_cells")
  observed$patient <- as.character(observed$patient)
  observed$label <- as.character(observed$label)
  x <- merge(full, observed, by = c("patient", "label"), all.x = TRUE, sort = FALSE)
  x$n_cells[is.na(x$n_cells)] <- 0L
  x$n_cells <- as.integer(x$n_cells)
  x$eligible <- x$n_cells >= minimum
  x$exclusion_reason <- ifelse(x$eligible, "", paste0("fewer_than_", minimum, "_cells"))
  names(x)[names(x) == "label"] <- label_name
  x
}

activity_from_aggregates <- function(aggregated, unit_ledger, label_name, yellow, orientation = 1) {
  keep_units <- unit_ledger$eligible
  keys <- paste(unit_ledger$patient[keep_units], unit_ledger[[label_name]][keep_units], sep = "::")
  if (!all(keys %in% colnames(aggregated))) stop("Eligible pseudobulk missing from aggregate matrix")
  pb <- aggregated[, keys, drop = FALSE]
  y <- DGEList(counts = pb)
  y <- calcNormFactors(y, method = "TMM")
  lcpm <- cpm(y, log = TRUE, prior.count = 2)
  genes <- intersect(yellow$gene, rownames(lcpm))
  if (!length(genes)) stop("No Yellow genes map to pseudobulk matrix")
  mat <- as.matrix(lcpm[genes, , drop = FALSE])
  gene_mean <- rowMeans(mat)
  gene_sd <- apply(mat, 1, sd)
  zero_var <- !is.finite(gene_sd) | gene_sd == 0
  z <- sweep(mat, 1, gene_mean, "-")
  nz <- !zero_var
  z[nz, ] <- sweep(z[nz, , drop = FALSE], 1, gene_sd[nz], "/")
  z[zero_var, ] <- 0
  w <- yellow$kME[match(genes, yellow$gene)] * orientation
  score <- colSums(z * w) / sum(abs(w))
  key_score <- setNames(as.numeric(score), colnames(pb))
  out <- unit_ledger
  out$yellow_activity <- NA_real_
  all_keys <- paste(out$patient, out[[label_name]], sep = "::")
  out$yellow_activity[out$eligible] <- key_score[all_keys[out$eligible]]
  qc <- data.frame(
    gene = genes,
    mean_log2cpm = gene_mean,
    sd_log2cpm = gene_sd,
    zero_variance = zero_var,
    signed_kME_used = w,
    stringsAsFactors = FALSE
  )
  list(activity = out, gene_qc = qc, tmm = data.frame(
    unit = colnames(pb), lib_size = y$samples$lib.size,
    norm_factor = y$samples$norm.factors, effective_lib_size = y$samples$lib.size * y$samples$norm.factors
  ))
}

summarise_activity <- function(x, label_name, levels_order) {
  eligible <- x[x$eligible & is.finite(x$yellow_activity), , drop = FALSE]
  rows <- lapply(levels_order, function(lab) {
    v <- eligible$yellow_activity[eligible[[label_name]] == lab]
    data.frame(label = lab, n_patients = length(v), median = median(v),
               q1 = unname(quantile(v, 0.25)), q3 = unname(quantile(v, 0.75)))
  })
  out <- do.call(rbind, rows)
  names(out)[1] <- label_name
  out
}

message("Reading and validating WGCNA authorities")
hub <- fread(inputs[["yellow"]], data.table = FALSE)
req_hub <- c("gene", "module", "kME", "is_top20")
if (!identical(names(hub), req_hub)) stop("Unexpected hub table schema")
yellow <- hub[hub$module == "yellow", , drop = FALSE]
if (nrow(yellow) != 539L || anyDuplicated(yellow$gene) || any(!is.finite(yellow$kME))) {
  stop("Yellow authority is not exactly 539 unique finite-kME genes")
}
vst <- readRDS(inputs[["vst"]])
if (is.null(rownames(vst)) || anyDuplicated(rownames(vst))) stop("Invalid WGCNA universe row names")
wgcna_universe <- rownames(vst)
if (length(wgcna_universe) != 5953L) stop("Unexpected WGCNA universe size")
rm(vst); gc()

major_levels <- c("B_cell", "CAF", "CD4_T_Helper", "CD8_T_Cytotoxic", "Cycling",
                  "DC_pDC", "Endothelial", "Macrophage", "Mast_cell", "Pericyte",
                  "Plasma_cell", "Treg", "Tumor_cell")
b_levels <- c("Naive B", "Memory B", "GC-like B", "Transitional B", "Plasma cell")

message("Reading whole-cohort single-cell object")
whole <- readRDS(inputs[["whole"]])
if (ncol(whole) != 102117L || length(unique(whole$patient)) != 18L) stop("Whole object size/patient contract failed")
if (!setequal(unique(as.character(whole$cell_type)), major_levels)) stop("Major cell-type contract failed")
whole_counts <- GetAssayData(whole, assay = "RNA", layer = "counts")
if (!identical(colnames(whole_counts), rownames(whole@meta.data))) stop("Whole counts/metadata cell order mismatch")
whole_patients <- sort(unique(as.character(whole$patient)))
whole_agg <- aggregate_counts(whole_counts, as.character(whole$patient), as.character(whole$cell_type))
whole_ledger <- make_unit_ledger(whole_patients, major_levels, whole_agg$observed, "cell_type", 20L)

mapping <- yellow
mapping$present_whole_RNA <- mapping$gene %in% rownames(whole_counts)
mapping$present_WGCNA_universe <- mapping$gene %in% wgcna_universe

# Establish orientation only from the prespecified B-cell anchors in Fig.2E.
tmp_activity <- activity_from_aggregates(whole_agg$counts, whole_ledger, "cell_type", yellow, orientation = 1)
anchor <- tmp_activity$activity$cell_type %in% c("B_cell", "Plasma_cell") & tmp_activity$activity$eligible
nonanchor <- !tmp_activity$activity$cell_type %in% c("B_cell", "Plasma_cell") & tmp_activity$activity$eligible
anchor_median_raw <- median(tmp_activity$activity$yellow_activity[anchor], na.rm = TRUE)
nonanchor_median_raw <- median(tmp_activity$activity$yellow_activity[nonanchor], na.rm = TRUE)
orientation <- ifelse(anchor_median_raw >= nonanchor_median_raw, 1, -1)
whole_activity <- if (orientation == 1) tmp_activity else
  activity_from_aggregates(whole_agg$counts, whole_ledger, "cell_type", yellow, orientation = orientation)

write_tsv(whole_activity$activity, file.path(table_dir, "Fig2E_major_lineage_activity.tsv"))
write_tsv(summarise_activity(whole_activity$activity, "cell_type", major_levels),
          file.path(table_dir, "Fig2E_major_lineage_summary.tsv"))
write_tsv(whole_activity$gene_qc, file.path(table_dir, "Fig2E_gene_activity_qc.tsv"))
write_tsv(whole_activity$tmm, file.path(table_dir, "Fig2E_tmm_factors.tsv"))

rm(whole, whole_counts, whole_agg, tmp_activity); gc()

message("Reading B-lineage object")
bobj <- readRDS(inputs[["bcell"]])
if (ncol(bobj) != 10572L || length(unique(bobj$patient)) != 18L) stop("B object size/patient contract failed")
if (!setequal(unique(as.character(bobj$b_subtype)), b_levels)) stop("B state-label contract failed")
bcounts <- GetAssayData(bobj, assay = "RNA", layer = "counts")
if (!identical(colnames(bcounts), rownames(bobj@meta.data))) stop("B counts/metadata cell order mismatch")
bpatients <- sort(unique(as.character(bobj$patient)))
bagg <- aggregate_counts(bcounts, as.character(bobj$patient), as.character(bobj$b_subtype))
bledger <- make_unit_ledger(bpatients, b_levels, bagg$observed, "b_subtype", 10L)
mapping$present_B_RNA <- mapping$gene %in% rownames(bcounts)

bactivity <- activity_from_aggregates(bagg$counts, bledger, "b_subtype", yellow, orientation = orientation)
write_tsv(bactivity$activity, file.path(table_dir, "Fig2F_b_state_activity.tsv"))
write_tsv(summarise_activity(bactivity$activity, "b_subtype", b_levels),
          file.path(table_dir, "Fig2F_b_state_summary.tsv"))
write_tsv(bactivity$gene_qc, file.path(table_dir, "Fig2F_gene_activity_qc.tsv"))
write_tsv(bactivity$tmm, file.path(table_dir, "Fig2F_tmm_factors.tsv"))
write_tsv(mapping, file.path(table_dir, "Yellow_539_exact_symbol_mapping.tsv"))

message("Running five donor-paired edgeR state-versus-rest comparisons")
de_list <- list()
comparison_qc <- list()
for (state in b_levels) {
  cell_tab <- table(factor(as.character(bobj$patient), levels = bpatients),
                    factor(as.character(bobj$b_subtype), levels = b_levels))
  target_n <- cell_tab[, state]
  rest_n <- rowSums(cell_tab[, setdiff(b_levels, state), drop = FALSE])
  eligible_patients <- bpatients[target_n >= 10L & rest_n >= 20L]
  if (length(eligible_patients) < 2L) {
    comparison_qc[[state]] <- data.frame(b_subtype = state, n_patient_pairs = length(eligible_patients),
                                         evaluable = FALSE, reason = "fewer_than_2_eligible_patient_pairs",
                                         n_tested_genes = 0L)
    next
  }
  target_cols <- paste(eligible_patients, state, sep = "::")
  if (!all(target_cols %in% colnames(bagg$counts))) stop("Missing target aggregate")
  target_mat <- bagg$counts[, target_cols, drop = FALSE]
  rest_cols <- lapply(eligible_patients, function(p) {
    paste(p, setdiff(b_levels, state), sep = "::")
  })
  rest_mat <- do.call(cbind, lapply(rest_cols, function(cols) {
    cols <- intersect(cols, colnames(bagg$counts))
    Matrix::rowSums(bagg$counts[, cols, drop = FALSE])
  }))
  rownames(rest_mat) <- rownames(bagg$counts)
  colnames(rest_mat) <- paste0(eligible_patients, "::rest")
  colnames(target_mat) <- paste0(eligible_patients, "::target")
  paired_counts <- cbind(rest_mat, target_mat)
  sample_meta <- data.frame(
    patient = factor(rep(eligible_patients, 2), levels = eligible_patients),
    group = factor(rep(c("rest", "target"), each = length(eligible_patients)), levels = c("rest", "target"))
  )
  rownames(sample_meta) <- colnames(paired_counts)
  design <- model.matrix(~ patient + group, sample_meta)
  y <- DGEList(counts = paired_counts)
  keep <- filterByExpr(y, design = design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y, method = "TMM")
  y <- estimateDisp(y, design, robust = TRUE)
  fit <- glmQLFit(y, design, robust = TRUE)
  qlf <- glmQLFTest(fit, coef = "grouptarget")
  tab <- topTags(qlf, n = Inf, sort.by = "none")$table
  tab$gene <- rownames(tab)
  tab$b_subtype <- state
  tab$n_patient_pairs <- length(eligible_patients)
  tab <- tab[, c("b_subtype", "gene", "n_patient_pairs", "logFC", "logCPM", "F", "PValue")]
  de_list[[state]] <- tab
  comparison_qc[[state]] <- data.frame(b_subtype = state, n_patient_pairs = length(eligible_patients),
                                       evaluable = TRUE, reason = "", n_tested_genes = nrow(tab))
}
de <- rbindlist(de_list, use.names = TRUE, fill = TRUE)
if (!nrow(de)) stop("No evaluable B-state comparisons")
de$global_BH <- p.adjust(de$PValue, method = "BH")
de$positive_marker <- de$logFC > 0 & de$global_BH < 0.05
setorder(de, b_subtype, PValue, gene)
write_tsv(de, file.path(table_dir, "Fig2G_edgeR_all_gene_state_tests.tsv"))
comparison_qc_df <- rbindlist(comparison_qc, use.names = TRUE, fill = TRUE)
write_tsv(comparison_qc_df, file.path(table_dir, "Fig2G_comparison_qc.tsv"))

enrichment_rows <- list()
marker_overlap <- list()
for (state in b_levels) {
  ds <- de[de$b_subtype == state, , drop = FALSE]
  qc <- comparison_qc_df[comparison_qc_df$b_subtype == state, , drop = FALSE]
  if (!nrow(ds)) next
  universe <- intersect(ds$gene, wgcna_universe)
  positive <- intersect(ds$gene[ds$positive_marker], universe)
  yellow_u <- intersect(yellow$gene, universe)
  overlap <- intersect(positive, yellow_u)
  a <- length(overlap)
  b <- length(setdiff(positive, yellow_u))
  c <- length(setdiff(yellow_u, positive))
  d <- length(setdiff(universe, union(positive, yellow_u)))
  ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2, byrow = TRUE), alternative = "two.sided")
  enrichment_rows[[state]] <- data.frame(
    b_subtype = state, n_patient_pairs = qc$n_patient_pairs,
    n_tested_universe = length(universe), n_positive_markers = length(positive),
    n_yellow_in_universe = length(yellow_u), n_overlap = a,
    a_marker_yellow = a, b_marker_nonyellow = b, c_nonmarker_yellow = c, d_neither = d,
    odds_ratio = unname(ft$estimate), conf_low = ft$conf.int[1], conf_high = ft$conf.int[2],
    p_value = ft$p.value, stringsAsFactors = FALSE
  )
  marker_overlap[[state]] <- data.frame(
    b_subtype = state, gene = universe,
    positive_marker = universe %in% positive,
    yellow_member = universe %in% yellow_u,
    overlap = universe %in% overlap,
    stringsAsFactors = FALSE
  )
}
enrich <- rbindlist(enrichment_rows, use.names = TRUE, fill = TRUE)
enrich$p_adjust <- p.adjust(enrich$p_value, method = "BH")
enrich$b_subtype <- factor(enrich$b_subtype, levels = b_levels)
enrich <- enrich[order(enrich$b_subtype), ]
enrich$b_subtype <- as.character(enrich$b_subtype)
write_tsv(enrich, file.path(table_dir, "Fig2G_yellow_marker_enrichment.tsv"))
write_tsv(rbindlist(marker_overlap), file.path(table_dir, "Fig2G_marker_yellow_membership.tsv"))

rm(bobj, bcounts, bagg); gc()

message("Mapping identifiers and running GO BP / Reactome ORA")
map_symbols <- sort(unique(c(wgcna_universe, yellow$gene)))
raw_map <- AnnotationDbi::select(org.Hs.eg.db, keys = map_symbols,
                                 columns = c("ENTREZID"), keytype = "SYMBOL")
raw_map <- unique(raw_map[, c("SYMBOL", "ENTREZID")])
all_symbol_rows <- merge(data.frame(SYMBOL = map_symbols), raw_map, by = "SYMBOL", all.x = TRUE, sort = FALSE)
map_counts <- aggregate(!is.na(all_symbol_rows$ENTREZID), by = list(SYMBOL = all_symbol_rows$SYMBOL), FUN = sum)
names(map_counts)[2] <- "n_entrez_mappings"
all_symbol_rows <- merge(all_symbol_rows, map_counts, by = "SYMBOL", all.x = TRUE, sort = FALSE)
all_symbol_rows$in_yellow <- all_symbol_rows$SYMBOL %in% yellow$gene
all_symbol_rows$in_wgcna_universe <- all_symbol_rows$SYMBOL %in% wgcna_universe
all_symbol_rows$mapping_status <- ifelse(all_symbol_rows$n_entrez_mappings == 0, "unmapped",
                                         ifelse(all_symbol_rows$n_entrez_mappings == 1, "one_to_one_or_many_symbols_to_one_entrez", "one_symbol_to_many_entrez"))
write_tsv(all_symbol_rows, file.path(table_dir, "Fig2H_symbol_entrez_mapping_ledger.tsv"))

query_entrez <- unique(na.omit(all_symbol_rows$ENTREZID[all_symbol_rows$in_yellow]))
universe_entrez <- unique(na.omit(all_symbol_rows$ENTREZID[all_symbol_rows$in_wgcna_universe]))
query_entrez <- intersect(query_entrez, universe_entrez)

go_obj <- enrichGO(gene = query_entrez, universe = universe_entrez, OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID", ont = "BP", pAdjustMethod = "BH",
                   pvalueCutoff = 1, qvalueCutoff = 1, minGSSize = 10, maxGSSize = 500,
                   readable = FALSE)
react_obj <- enrichPathway(gene = query_entrez, universe = universe_entrez, organism = "human",
                           pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1,
                           minGSSize = 10, maxGSSize = 500, readable = FALSE)

ratio_num <- function(x) {
  parts <- strsplit(as.character(x), "/", fixed = TRUE)
  vapply(parts, function(z) as.numeric(z[1]) / as.numeric(z[2]), numeric(1))
}
normalize_ora <- function(obj, database) {
  x <- as.data.frame(obj)
  if (!nrow(x)) return(data.frame())
  x$database <- database
  x$enrichment_ratio <- ratio_num(x$GeneRatio) / ratio_num(x$BgRatio)
  x
}
go_res <- normalize_ora(go_obj, "GO Biological Process")
react_res <- normalize_ora(react_obj, "Reactome")
write_tsv(go_res, file.path(table_dir, "Fig2H_GO_BP_complete_ORA.tsv"))
write_tsv(react_res, file.path(table_dir, "Fig2H_Reactome_complete_ORA.tsv"))

make_top <- function(x) {
  if (!nrow(x)) return(x)
  x <- x[order(x$p.adjust, -x$enrichment_ratio, x$Description), , drop = FALSE]
  x <- head(x, 5)
  x$display_rank <- seq_len(nrow(x))
  x[, c("database", "ID", "Description", "GeneRatio", "BgRatio", "enrichment_ratio",
       "Count", "pvalue", "p.adjust", "display_rank")]
}
top_terms <- rbind(make_top(go_res), make_top(react_res))
names(top_terms)[names(top_terms) == "pvalue"] <- "p_value"
names(top_terms)[names(top_terms) == "p.adjust"] <- "p_adjust"
write_tsv(top_terms, file.path(table_dir, "Fig2H_top_terms.tsv"))

mapping_summary <- data.frame(
  set = c("Yellow symbols", "WGCNA universe symbols"),
  n_symbols = c(length(unique(yellow$gene)), length(unique(wgcna_universe))),
  n_symbols_mapped = c(length(unique(all_symbol_rows$SYMBOL[all_symbol_rows$in_yellow & !is.na(all_symbol_rows$ENTREZID)])),
                       length(unique(all_symbol_rows$SYMBOL[all_symbol_rows$in_wgcna_universe & !is.na(all_symbol_rows$ENTREZID)]))),
  n_symbols_unmapped = c(length(unique(yellow$gene)) - length(unique(all_symbol_rows$SYMBOL[all_symbol_rows$in_yellow & !is.na(all_symbol_rows$ENTREZID)])),
                         length(unique(wgcna_universe)) - length(unique(all_symbol_rows$SYMBOL[all_symbol_rows$in_wgcna_universe & !is.na(all_symbol_rows$ENTREZID)]))),
  n_unique_entrez = c(length(query_entrez), length(universe_entrez)),
  n_one_symbol_to_many_entrez = c(sum(map_counts$n_entrez_mappings[match(yellow$gene, map_counts$SYMBOL)] > 1, na.rm = TRUE),
                                    sum(map_counts$n_entrez_mappings[match(wgcna_universe, map_counts$SYMBOL)] > 1, na.rm = TRUE))
)
write_tsv(mapping_summary, file.path(table_dir, "Fig2H_mapping_summary.tsv"))

message("Writing provenance, validation, report and manifest")
input_hashes <- data.frame(input_name = names(inputs), path = unname(inputs),
                           sha256 = vapply(inputs, sha256_file, character(1)), stringsAsFactors = FALSE)
write_tsv(input_hashes, file.path(report_dir, "INPUT_SHA256.tsv"))

orientation_record <- data.frame(
  raw_anchor_median = anchor_median_raw,
  raw_nonanchor_median = nonanchor_median_raw,
  orientation_multiplier = orientation,
  rule = "Set B_cell/Plasma_cell anchor median >= non-anchor median; same fixed orientation used in Fig2F"
)
write_tsv(orientation_record, file.path(report_dir, "ACTIVITY_ORIENTATION.tsv"))

validation <- data.frame(
  check = c(
    "Yellow has exactly 539 unique finite-kME genes", "Whole object has 102117 cells and 18 patients",
    "Whole object has registered 13 cell types", "B object has 10572 cells and 18 patients",
    "B object has registered five provisional states", "No cell-level inferential test used",
    "Fig2G uses paired donor design", "Global BH applied across all gene-state tests",
    "Five Fisher tests use BH correction", "ORA uses exact WGCNA universe and database-specific BH",
    "No bootstrap used", "No plotting performed by producer"
  ),
  pass = c(
    nrow(yellow) == 539L && !anyDuplicated(yellow$gene) && all(is.finite(yellow$kME)),
    length(whole_patients) == 18L,
    length(major_levels) == 13L,
    length(bpatients) == 18L,
    length(b_levels) == 5L,
    TRUE, all(comparison_qc_df$n_patient_pairs[comparison_qc_df$evaluable] >= 2),
    "global_BH" %in% names(de), nrow(enrich) == 5L && "p_adjust" %in% names(enrich),
    all(c("GO Biological Process", "Reactome") %in% unique(top_terms$database)), TRUE, TRUE
  ), stringsAsFactors = FALSE
)
write_tsv(validation, file.path(report_dir, "PRODUCTION_VALIDATION.tsv"))
if (!all(validation$pass)) stop("One or more production validation checks failed")

capture.output(sessionInfo(), file = file.path(report_dir, "SESSION_INFO.txt"))
report <- c(
  "# Replacement Fig.2E-H Yellow-module numeric production",
  "",
  "## Contract",
  "All methods, thresholds, universes and panel roles follow `FIG2_YELLOW_PANEL_CONTRACT.md`, frozen before outcome inspection.",
  "Patient/donor is the independent biological unit. No cell-level inferential test and no bootstrap were used.",
  "",
  "## Activity coverage",
  sprintf("- Yellow genes mapped to whole RNA assay: %d/539.", sum(mapping$present_whole_RNA)),
  sprintf("- Yellow genes mapped to B-lineage RNA assay: %d/539.", sum(mapping$present_B_RNA)),
  sprintf("- Fig.2E eligible patient-cell-type units: %d/%d.", sum(whole_activity$activity$eligible), nrow(whole_activity$activity)),
  sprintf("- Fig.2F eligible patient-state units: %d/%d.", sum(bactivity$activity$eligible), nrow(bactivity$activity)),
  sprintf("- Fixed orientation multiplier: %d.", orientation),
  "",
  "## Marker enrichment",
  sprintf("- Evaluable paired state comparisons: %d/5.", sum(comparison_qc_df$evaluable)),
  sprintf("- Total tested gene-state hypotheses under one global BH family: %d.", nrow(de)),
  "- Fisher enrichment is two-sided and BH-adjusted across five state tests.",
  "",
  "## Functional annotation",
  sprintf("- Yellow symbols mapped to unique Entrez IDs within universe: %d.", length(query_entrez)),
  sprintf("- WGCNA universe symbols mapped to unique Entrez IDs: %d.", length(universe_entrez)),
  sprintf("- Complete GO BP results: %d terms; complete Reactome results: %d terms.", nrow(go_res), nrow(react_res)),
  "",
  "## Claim ceiling",
  "These outputs support cross-scale descriptive coherence only. Historical B-state labels remain provisional; results do not validate sparse4, establish B-cell specificity or support causal claims."
)
writeLines(report, file.path(report_dir, "PRODUCTION_REPORT.md"))

files <- list.files(out_dir, recursive = TRUE, full.names = TRUE)
files <- files[file.info(files)$isdir %in% FALSE]
files <- files[basename(files) != "MANIFEST_SHA256.tsv"]
manifest <- data.frame(path = sub(paste0("^", out_dir, "/"), "", files),
                       bytes = file.info(files)$size,
                       sha256 = vapply(files, sha256_file, character(1)), stringsAsFactors = FALSE)
manifest <- manifest[order(manifest$path), ]
write_tsv(manifest, file.path(out_dir, "MANIFEST_SHA256.tsv"))
message("Production complete: ", nrow(manifest), " files recorded")
