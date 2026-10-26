#!/usr/bin/env Rscript
# ==========================================================
# Script   : 23_make_figures.R
# Purpose  : Verify released panel assets and report their source generators
# Inputs   : docs/figure_table_map.csv; results/figures/
# Outputs  : logs/23_make_figures_verification.tsv
# Paper    : Figures 1-7 and Supplementary Figures 1-4
# Seed     : not applicable
# Runtime  : less than one minute in verification mode
# Last run : 2026-10-26
# ==========================================================

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
mapping <- read.csv(file.path(root, "docs", "figure_table_map.csv"), check.names = FALSE)
figure_rows <- mapping[mapping$item_type == "figure_panel", , drop = FALSE]
exists <- vapply(figure_rows$output_file, function(x) {
  if (!nzchar(x)) return(FALSE)
  all(file.exists(file.path(root, strsplit(x, ";", fixed = TRUE)[[1]])))
}, logical(1))
report <- data.frame(paper_item = figure_rows$paper_item, released_asset_present = exists,
                     generator = figure_rows$script, note = figure_rows$description,
                     stringsAsFactors = FALSE)
dir.create(file.path(root, "logs"), showWarnings = FALSE)
write.table(report, file.path(root, "logs", "23_make_figures_verification.tsv"),
            sep = "	", quote = FALSE, row.names = FALSE)
message(sum(exists), " of ", nrow(report), " mapped panels have released binary assets. See the map for omitted large or author-edited panels.")
