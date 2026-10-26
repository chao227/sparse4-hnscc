#!/usr/bin/env Rscript

# ==========================================================
# Script   : fig2c_production_render_fig2c.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_repairs/Fig2C_production/render_fig2c.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(WGCNA)
  library(ragg)
  library(digest)
})

options(stringsAsFactors = FALSE)
need <- function(ok, msg, ...) if (!isTRUE(ok)) stop(sprintf(msg, ...), call. = FALSE)
write_tsv <- function(x, path) write.table(x, path, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag) {
  i <- match(flag, args)
  need(!is.na(i) && i < length(args), "Missing argument: %s", flag)
  args[[i + 1L]]
}
root <- normalizePath(get_arg("--root"), mustWork = TRUE)
out <- normalizePath(get_arg("--out"), mustWork = TRUE)
need(basename(root) == "reback", "Root must be the reback workspace")
need(startsWith(paste0(out, .Platform$file.sep), paste0(root, .Platform$file.sep)), "Output escapes reback")

input_path <- file.path(root, "inputs/wgcna/legacy/wgcna_results.rds")
need(file.exists(input_path) && !file.info(input_path)$isdir, "Locked WGCNA result object is missing")
input_sha <- digest(file = input_path, algo = "sha256", serialize = FALSE)
wg <- readRDS(input_path)
need(all(c("moduleTraitCor", "moduleTraitP") %in% names(wg)), "WGCNA result lacks module-trait matrices")
cor_mat <- wg$moduleTraitCor
p_mat <- wg$moduleTraitP
need(is.matrix(cor_mat) && is.matrix(p_mat) && identical(dim(cor_mat), c(10L, 26L)), "Expected 10 x 26 matrices")
need(identical(dimnames(cor_mat), dimnames(p_mat)), "Correlation and P-value dimnames differ")
need(all(is.finite(cor_mat)) && all(is.finite(p_mat)), "Nonfinite matrix value detected")
need(all(cor_mat >= -1 & cor_mat <= 1) && all(p_mat >= 0 & p_mat <= 1), "Correlation/P-value range failed")

module_order <- rownames(cor_mat)
trait_order <- colnames(cor_mat)
idx <- expand.grid(module_index = seq_len(nrow(cor_mat)), trait_index = seq_len(ncol(cor_mat)))
dat <- data.frame(
  module = module_order[idx$module_index],
  trait = trait_order[idx$trait_index],
  module_order = idx$module_index,
  trait_order = idx$trait_index,
  correlation = cor_mat[cbind(idx$module_index, idx$trait_index)],
  p_value = p_mat[cbind(idx$module_index, idx$trait_index)],
  stringsAsFactors = FALSE
)
compact_leading_zero <- function(x) sub("^(-?)0\\.", "\\1.", x)
# The exact full-precision values and the original significant-digit rendering
# are retained in source data. Fixed two-decimal r labels prevent rare near-zero
# scientific notation from colliding at the immutable 130-mm panel width.
dat$original_display_correlation <- as.character(signif(dat$correlation, 2))
dat$original_display_p_value <- as.character(signif(dat$p_value, 1))
dat$display_correlation <- compact_leading_zero(sprintf("%.2f", dat$correlation))
dat$display_correlation[dat$display_correlation %in% c("-.00", ".00")] <- ".00"
dat$display_p_value <- compact_leading_zero(dat$original_display_p_value)
dat$display_label <- paste0(dat$display_correlation, "\n(", dat$display_p_value, ")")
dat$text_colour <- ifelse(abs(dat$correlation) >= 0.62, "white", "black")
dat$module <- factor(dat$module, levels = module_order)
dat$trait <- factor(dat$trait, levels = rev(trait_order))
need(nrow(dat) == 260L && !anyDuplicated(paste(dat$module, dat$trait)), "Long-form cell contract failed")

source_out <- transform(dat,
  module = as.character(module), trait = as.character(trait),
  correlation = sprintf("%.17g", correlation), p_value = sprintf("%.17g", p_value)
)
source_out$display_label <- gsub("\n", "\\\\n", source_out$display_label, fixed = TRUE)
source_out <- source_out[, c("module", "trait", "module_order", "trait_order", "correlation", "p_value", "original_display_correlation", "original_display_p_value", "display_correlation", "display_p_value", "display_label")]
write_tsv(source_out, file.path(out, "Fig2C_source_data.tsv"))

palette <- WGCNA::blueWhiteRed(50)
pt_to_mm <- function(pt) pt / (72 / 25.4)
p <- ggplot(dat, aes(x = module, y = trait, fill = correlation)) +
  geom_tile(colour = "grey88", linewidth = 0.16) +
  geom_text(aes(label = display_label, colour = text_colour),
            family = "Arial", size = pt_to_mm(7), lineheight = 0.61,
            fontface = "plain", show.legend = FALSE) +
  scale_colour_identity() +
  scale_fill_gradientn(
    colours = palette, limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1),
    name = "Correlation", guide = guide_colourbar(
      title.position = "top", title.hjust = 0.5, barwidth = unit(3.4, "mm"),
      barheight = unit(58, "mm"), ticks.colour = "black", frame.colour = "black"
    )
  ) +
  scale_x_discrete(position = "bottom", expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  coord_cartesian(clip = "off") +
  labs(title = "Module–Trait Relationships", x = NULL, y = NULL) +
  theme_minimal(base_family = "Arial", base_size = 7) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1, colour = "black"),
    axis.text.y = element_text(size = 7, colour = "black", margin = margin(r = 1.4, unit = "mm")),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 7, face = "plain", hjust = 0.5),
    legend.text = element_text(size = 7, colour = "black"),
    legend.margin = margin(l = 0.5, unit = "mm"),
    legend.box.margin = margin(0, 0, 0, 0),
    plot.title = element_text(size = 9, face = "bold", hjust = 0.5, margin = margin(b = 1.2, unit = "mm")),
    plot.margin = margin(t = 5.0, r = 1.5, b = 1.5, l = 1.5, unit = "mm")
  )

width_mm <- 130
height_mm <- 120
pdf_path <- file.path(out, "Fig2C_module_trait_heatmap.pdf")
png_path <- file.path(out, "Fig2C_module_trait_heatmap.png")
svg_path <- file.path(out, "Fig2C_module_trait_heatmap.svg")
grDevices::cairo_pdf(pdf_path, width = width_mm / 25.4, height = height_mm / 25.4, family = "Arial")
print(p)
dev.off()
ragg::agg_png(png_path, width = width_mm, height = height_mm, units = "mm", res = 600, background = "white")
print(p)
dev.off()
svglite::svglite(svg_path, width = width_mm / 25.4, height = height_mm / 25.4, bg = "white")
print(p)
dev.off()

receipt <- data.frame(
  item = c("input_wgcna_rds", "matrix_shape", "cells", "correlation_range", "p_value_range", "width_mm", "height_mm", "cell_text_pt", "axis_text_pt", "legend_text_pt", "title_pt", "backend", "estimand_change", "row_exclusion"),
  observed = c(input_sha, paste(dim(cor_mat), collapse = "x"), nrow(dat), paste(range(cor_mat), collapse = ";"), paste(range(p_mat), collapse = ";"), width_mm, height_mm, 7, 7, 7, 9, "R/ggplot2 + Cairo PDF + ragg PNG", "none", "none"),
  acceptance_rule = c("recorded SHA256", "10x26", "260", "within [-1,1]", "within [0,1]", "130", "120", ">=7", ">=7", ">=7", ">=7", "R only", "unchanged", "all cells retained"),
  status = "PASS", stringsAsFactors = FALSE
)
write_tsv(receipt, file.path(out, "Fig2C_production_receipt.tsv"))
writeLines(capture.output(sessionInfo()), file.path(out, "SESSIONINFO.txt"))
cat("Fig2C rendered from 260 locked cells at 130 x 120 mm\n")
