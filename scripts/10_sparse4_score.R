#!/usr/bin/env Rscript

# ==========================================================
# Script   : 10_sparse4_score.R
# Purpose  : Score compatible expression matrices with the locked four-gene coefficients
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: scoring specification; Supplementary Table 2
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/scoring/score_sparse4.R
# ==========================================================
# Score the locked four-gene sparse4 model on TCGA-compatible DESeq2-VST values.
coef <- c(CAMK2N1 = 0.0123106633921438,
          MASP1 = -0.0175920351043370,
          `MIR9-3HG` = -0.0039450538526462,
          MS4A1 = -0.0218626023395592)

score_sparse4_vst <- function(samples_by_genes) {
  missing <- setdiff(names(coef), colnames(samples_by_genes))
  if (length(missing)) stop("FULL sparse4 score requires all four genes; missing: ", paste(missing, collapse = ", "))
  x <- as.matrix(samples_by_genes[, names(coef), drop = FALSE])
  storage.mode(x) <- "double"
  drop(x %*% coef)
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  if (length(args) != 2L) stop("usage: score_sparse4.R INPUT.tsv OUTPUT.tsv")
  dat <- read.delim(args[[1L]], check.names = FALSE)
  if (!"sample_id" %in% names(dat)) stop("input requires sample_id")
  out <- data.frame(sample_id = dat$sample_id, sparse4_linear_predictor = score_sparse4_vst(dat))
  dir.create(dirname(args[[2L]]), recursive = TRUE, showWarnings = FALSE)
  write.table(out, args[[2L]], sep = "\t", quote = FALSE, row.names = FALSE)
}
