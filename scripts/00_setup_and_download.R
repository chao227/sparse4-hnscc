#!/usr/bin/env Rscript
# ==========================================================
# Script   : 00_setup_and_download.R
# Purpose  : Validate public accessions and prepare local data directories
# Inputs   : data/accessions.csv
# Outputs  : data/raw/, data/processed/, logs/00_setup_and_download_session_info.txt
# Paper    : Methods > Data sources
# Seed     : not applicable
# Runtime  : less than one minute without data transfer
# Last run : 2026-10-26
# ==========================================================

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
manifest <- read.csv(file.path(root, "data", "accessions.csv"), check.names = FALSE)
required <- c("dataset", "accession", "platform", "n_patients", "role", "source_url", "download_date", "version_note")
missing <- setdiff(required, names(manifest))
if (length(missing)) stop("Accessions file is missing: ", paste(missing, collapse = ", "), call. = FALSE)
for (dataset in manifest$dataset) {
  dir.create(file.path(root, "data", "raw", dataset), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(root, "data", "processed", dataset), recursive = TRUE, showWarnings = FALSE)
}
if ("--download" %in% args) {
  stop("Automatic bulk download is disabled because repository landing pages do not uniquely identify all deposited files. Follow data/README.md and record checksums after retrieval.", call. = FALSE)
}
source(file.path(root, "R", "utils_logging.R"))
write_session_info("00_setup_and_download.R", root)
message("Prepared directories for ", nrow(manifest), " public datasets. No data were downloaded.")
