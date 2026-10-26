#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_render_fig7c_traceable.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/render_fig7c_traceable.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
out <- file.path(root, "panel_layout_revision_20261020")
figdir <- file.path(out, "figures")
srcdir <- file.path(out, "source_data")
dir.create(figdir, recursive = TRUE, showWarnings = FALSE)
dir.create(srcdir, recursive = TRUE, showWarnings = FALSE)

vst <- readRDS(file.path(root, "inputs/tcga/expression/vst_matrix.rds"))
tme <- read.csv(file.path(root, "inputs/tcga/features/tme_scores.csv"),
                check.names = FALSE, stringsAsFactors = FALSE)
surv <- read.csv(file.path(root, "inputs/tcga/clinical/survival_data.csv"),
                 check.names = FALSE, stringsAsFactors = FALSE)

stopifnot("IKZF1" %in% rownames(vst), nrow(surv) == 476L,
          !anyDuplicated(surv$sample_id), !anyDuplicated(tme$sample_id))
ix <- match(surv$sample_id, colnames(vst))
it <- match(substr(surv$sample_id, 1, 15), tme$sample_id)
stopifnot(!anyNA(ix), !anyNA(it))

d <- data.frame(
  sample_id = surv$sample_id,
  patient_sample_id = substr(surv$sample_id, 1, 15),
  IKZF1_VST = as.numeric(vst["IKZF1", ix]),
  ImmuneScore = tme$ImmuneScore[it],
  stringsAsFactors = FALSE
)
rho <- unname(cor(d$IKZF1_VST, d$ImmuneScore, method = "spearman"))
pval <- suppressWarnings(cor.test(d$IKZF1_VST, d$ImmuneScore,
                                  method = "spearman", exact = FALSE)$p.value)
stopifnot(nrow(d) == 476L, isTRUE(all.equal(rho, 0.8919581, tolerance = 1e-7)))

p <- ggplot(d, aes(IKZF1_VST, ImmuneScore)) +
  geom_point(size = 0.60, alpha = 0.52, colour = "#4C91C6") +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE,
              linewidth = 0.50, colour = "#17365D") +
  annotate("text", x = 5.0, y = 0.77, hjust = 0, vjust = 1,
           label = sprintf("n = %d\nSpearman r = %.3f\nP < 0.001", nrow(d), rho),
           family = "Arial", size = 2.50, lineheight = 1.05) +
  labs(title = "IKZF1 expression and\nproject-defined immune\nssGSEA score",
       x = "IKZF1 expression\n(VST, TCGA)",
       y = "Project immune ssGSEA score") +
  theme_classic(base_family = "Arial", base_size = 7.5) +
  theme(plot.title = element_text(size = 7.2, face = "bold", hjust = 0.5,
                                  lineheight = 0.95),
        axis.title = element_text(size = 7.5),
        axis.text = element_text(size = 7.0, colour = "black"),
        plot.margin = margin(5, 5, 4, 5))

pdf_file <- file.path(figdir, "Fig7C_IKZF1_ImmuneScore_traceable.pdf")
png_file <- file.path(figdir, "Fig7C_IKZF1_ImmuneScore_traceable.png")
grDevices::cairo_pdf(pdf_file, width = 50 / 25.4, height = 70 / 25.4,
                     family = "Arial", onefile = FALSE)
print(p)
dev.off()
ragg::agg_png(png_file, width = 50, height = 70, units = "mm", res = 600,
              background = "white")
print(p)
dev.off()

write.table(d, file.path(srcdir, "Fig7C_IKZF1_ImmuneScore_source.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(n = nrow(d), spearman_rho = rho, p_value = pval,
                       independent_unit = "TCGA patient/sample"),
            file.path(srcdir, "Fig7C_IKZF1_ImmuneScore_statistics.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
capture.output(sessionInfo(), file = file.path(out, "logs/Fig7C_sessionInfo.txt"))
