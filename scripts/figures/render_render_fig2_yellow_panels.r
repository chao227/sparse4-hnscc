#!/usr/bin/env Rscript

# ==========================================================
# Script   : render_render_fig2_yellow_panels.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:replacement_fig2_yellow_panels_20261019/render/scripts/render_fig2_yellow_panels.R
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(ragg)
})

options(stringsAsFactors = FALSE)

root <- normalizePath("replacement_fig2_yellow_panels_20261019", mustWork = TRUE)
table_dir <- file.path(root, "production", "tables")
figure_dir <- file.path(root, "render", "figures")
report_dir <- file.path(root, "render", "reports")
log_dir <- file.path(root, "render", "logs")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(report_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

read_tsv_strict <- function(name, required) {
  path <- file.path(table_dir, name)
  if (!file.exists(path)) stop("Missing production source table: ", path)
  x <- read.delim(path, check.names = FALSE, quote = "", comment.char = "")
  missing <- setdiff(required, names(x))
  if (length(missing)) stop(name, " is missing columns: ", paste(missing, collapse = ", "))
  x
}

assert_finite <- function(x, columns, name) {
  for (column in columns) {
    if (any(!is.finite(x[[column]]))) stop(name, " has non-finite ", column)
  }
}

theme_panel <- function() {
  theme_classic(base_size = 7.5, base_family = "Arial") +
    theme(
      axis.line = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks = element_line(linewidth = 0.3, colour = "black"),
      axis.ticks.length = grid::unit(1.2, "mm"),
      axis.title = element_text(size = 7.5, colour = "black"),
      axis.text = element_text(size = 7.5, colour = "black"),
      plot.title = element_text(size = 8.2, face = "bold", colour = "black", hjust = 0.5),
      strip.background = element_blank(),
      strip.text = element_text(size = 7.5, face = "bold", colour = "black", hjust = 0),
      legend.title = element_text(size = 7.5),
      legend.text = element_text(size = 7.5),
      legend.key.height = grid::unit(2.5, "mm"),
      legend.key.width = grid::unit(2.8, "mm"),
      panel.grid = element_blank(),
      plot.title.position = "plot",
      plot.margin = margin(2.0, 2.0, 1.5, 1.5, unit = "mm")
    )
}

export_pair <- function(plot, stem, width_mm, height_mm) {
  pdf_path <- file.path(figure_dir, paste0(stem, ".pdf"))
  png_path <- file.path(figure_dir, paste0(stem, ".png"))
  grDevices::cairo_pdf(
    filename = pdf_path,
    width = width_mm / 25.4,
    height = height_mm / 25.4,
    family = "Arial",
    onefile = TRUE
  )
  print(plot)
  grDevices::dev.off()
  ragg::agg_png(
    filename = png_path,
    width = width_mm,
    height = height_mm,
    units = "mm",
    res = 600,
    background = "white"
  )
  print(plot)
  grDevices::dev.off()
  invisible(c(pdf_path, png_path))
}

fmt_q <- function(x) {
  ifelse(x < 0.001, format(x, scientific = TRUE, digits = 2), sprintf("%.3f", x))
}

# Fig.2E: eligible donor x lineage pseudobulks only; no inferential test.
e_points <- read_tsv_strict(
  "Fig2E_major_lineage_activity.tsv",
  c("patient", "cell_type", "n_cells", "eligible", "exclusion_reason", "yellow_activity")
)
e_summary <- read_tsv_strict(
  "Fig2E_major_lineage_summary.tsv",
  c("cell_type", "n_patients", "median", "q1", "q3")
)
e_points <- e_points[e_points$eligible %in% c(TRUE, "TRUE", 1, "1"), , drop = FALSE]
if (!nrow(e_points)) stop("Fig.2E has no eligible pseudobulks")
assert_finite(e_points, c("yellow_activity"), "Fig.2E activity")
assert_finite(e_summary, c("n_patients", "median", "q1", "q3"), "Fig.2E summary")
if (nrow(e_summary) != 13L) stop("Fig.2E summary must contain exactly 13 registered cell types")
e_order <- as.character(e_summary$cell_type)
e_display <- gsub("_", " ", e_summary$cell_type, fixed = TRUE)
e_display[e_summary$cell_type == "DC_pDC"] <- "DC/pDC"
e_labels <- setNames(sprintf("%s  (n=%d)", e_display, e_summary$n_patients), e_summary$cell_type)
e_points$cell_type <- factor(e_points$cell_type, levels = rev(e_order))
e_summary$cell_type <- factor(e_summary$cell_type, levels = rev(e_order))
p_e <- ggplot(e_points, aes(y = cell_type, x = yellow_activity)) +
  geom_point(
    position = position_jitter(width = 0, height = 0.10, seed = 20261019),
    shape = 16, size = 0.75, alpha = 0.55, colour = "#6B9AC4"
  ) +
  geom_errorbar(
    data = e_summary,
    aes(y = cell_type, xmin = q1, xmax = q3),
    inherit.aes = FALSE, orientation = "y", width = 0, linewidth = 0.7, colour = "#1D3557"
  ) +
  geom_point(
    data = e_summary,
    aes(y = cell_type, x = median),
    inherit.aes = FALSE, shape = 21, size = 1.55, stroke = 0.4,
    fill = "white", colour = "#1D3557"
  ) +
  scale_y_discrete(labels = e_labels) +
  labs(
    title = "Yellow activity: lineages",
    x = "Signed-kME activity", y = NULL
  ) +
  theme_panel() +
  theme(axis.text.y = element_text(size = 7.5), plot.margin = margin(2, 2, 1.5, 1, unit = "mm"))

# Fig.2F: eligible donor x provisional B-state pseudobulks only.
f_points <- read_tsv_strict(
  "Fig2F_b_state_activity.tsv",
  c("patient", "b_subtype", "n_cells", "eligible", "exclusion_reason", "yellow_activity")
)
f_summary <- read_tsv_strict(
  "Fig2F_b_state_summary.tsv",
  c("b_subtype", "n_patients", "median", "q1", "q3")
)
f_points <- f_points[f_points$eligible %in% c(TRUE, "TRUE", 1, "1"), , drop = FALSE]
if (!nrow(f_points)) stop("Fig.2F has no eligible pseudobulks")
assert_finite(f_points, c("yellow_activity"), "Fig.2F activity")
assert_finite(f_summary, c("n_patients", "median", "q1", "q3"), "Fig.2F summary")
b_order <- c("Naive B", "Memory B", "GC-like B", "Transitional B", "Plasma cell")
if (!setequal(as.character(f_summary$b_subtype), b_order)) stop("Fig.2F state labels differ from the five frozen provisional states")
f_labels <- setNames(sprintf("%s  (n=%d)", f_summary$b_subtype, f_summary$n_patients), f_summary$b_subtype)
f_points$b_subtype <- factor(f_points$b_subtype, levels = rev(b_order))
f_summary$b_subtype <- factor(f_summary$b_subtype, levels = rev(b_order))
p_f <- ggplot(f_points, aes(y = b_subtype, x = yellow_activity)) +
  geom_point(
    position = position_jitter(width = 0, height = 0.11, seed = 20261019),
    shape = 16, size = 0.85, alpha = 0.6, colour = "#D08C60"
  ) +
  geom_errorbar(
    data = f_summary,
    aes(y = b_subtype, xmin = q1, xmax = q3),
    inherit.aes = FALSE, orientation = "y", width = 0, linewidth = 0.8, colour = "#7F5539"
  ) +
  geom_point(
    data = f_summary,
    aes(y = b_subtype, x = median),
    inherit.aes = FALSE, shape = 21, size = 1.7, stroke = 0.45,
    fill = "white", colour = "#7F5539"
  ) +
  scale_y_discrete(labels = f_labels) +
  labs(
    title = "Yellow activity in B states",
    x = "Signed-kME activity", y = NULL
  ) +
  theme_panel() +
  theme(axis.text.y = element_text(size = 7.5), plot.margin = margin(2, 2, 1.5, 1, unit = "mm"))

# Fig.2G: producer-supplied Fisher estimates and five-test global BH q values.
g <- read_tsv_strict(
  "Fig2G_yellow_marker_enrichment.tsv",
  c("b_subtype", "n_patient_pairs", "n_tested_universe", "n_positive_markers",
    "n_yellow_in_universe", "n_overlap", "a_marker_yellow", "b_marker_nonyellow",
    "c_nonmarker_yellow", "d_neither", "odds_ratio", "conf_low", "conf_high",
    "p_value", "p_adjust")
)
if (!setequal(as.character(g$b_subtype), b_order)) stop("Fig.2G state labels differ from the five frozen provisional states")
assert_finite(g, c("odds_ratio", "conf_low", "conf_high", "p_value", "p_adjust"), "Fig.2G")
if (any(g$odds_ratio <= 0 | g$conf_low <= 0 | g$conf_high <= 0)) stop("Fig.2G log-scale estimates must be positive")
if (any(g$conf_low > g$odds_ratio | g$conf_high < g$odds_ratio)) stop("Fig.2G estimate lies outside its confidence interval")
g$display_label <- sprintf("%s\n%d/%d; q=%s", g$b_subtype, g$n_overlap, g$n_positive_markers, fmt_q(g$p_adjust))
g$display_label <- factor(g$display_label, levels = rev(g$display_label[match(b_order, g$b_subtype)]))
p_g <- ggplot(g, aes(y = display_label, x = odds_ratio, xmin = conf_low, xmax = conf_high)) +
  geom_vline(xintercept = 1, linewidth = 0.35, linetype = 2, colour = "#777777") +
  geom_errorbar(orientation = "y", width = 0, linewidth = 0.7, colour = "#264653") +
  geom_point(shape = 21, size = 1.8, stroke = 0.45, fill = "#6AB0A3", colour = "#264653") +
  scale_x_log10() +
  labs(
    title = "Marker overlap",
    x = "Odds ratio (95% CI; log scale)", y = NULL
  ) +
  theme_panel() +
  theme(axis.text.y = element_text(size = 7.5, lineheight = 0.88), plot.margin = margin(2, 2, 1.5, 1, unit = "mm"))

# Fig.2H: supplied top-five-per-database table; no enrichment recomputation.
h <- read_tsv_strict(
  "Fig2H_top_terms.tsv",
  c("database", "ID", "Description", "GeneRatio", "BgRatio", "enrichment_ratio",
    "Count", "p_value", "p_adjust", "display_rank")
)
if (!nrow(h)) stop("Fig.2H top-term table is empty")
if (any(table(h$database) > 5L)) stop("Fig.2H contains more than five displayed terms per database")
assert_finite(h, c("enrichment_ratio", "Count", "p_value", "p_adjust", "display_rank"), "Fig.2H")
if (any(h$enrichment_ratio <= 0 | h$Count <= 0 | h$p_adjust < 0 | h$p_adjust > 1)) stop("Fig.2H contains invalid enrichment values")
short_term <- function(x) {
  replacements <- c(
    "production of molecular mediator of immune response" = "Immune mediator production",
    "B cell mediated immunity" = "B-cell immunity",
    "immunoglobulin mediated immune response" = "Ig-mediated immunity",
    "immunoglobulin production" = "Ig production",
    "adaptive immune response" = "Adaptive immunity",
    "Antigen activates B Cell Receptor (BCR) leading to generation of second messengers" = "BCR second messengers",
    "Signaling by the B Cell Receptor (BCR)" = "BCR signaling",
    "Disorders of transmembrane transporters" = "TM transporter disorders",
    "ABC transporter disorders" = "ABC transporter disorders",
    "Adrenaline,noradrenaline inhibits insulin secretion" = "Adrenergic insulin inhibition"
  )
  out <- x
  hit <- match(x, names(replacements))
  out[!is.na(hit)] <- unname(replacements[hit[!is.na(hit)]])
  out
}
h <- h[order(h$database, h$display_rank, h$p_adjust, -h$enrichment_ratio), , drop = FALSE]
h$term_key <- paste(h$database, h$ID, sep = "___")
h$display_description <- short_term(h$Description)
h$display_wrapped <- h$display_description
h$term_axis <- factor(
  h$term_key,
  levels = rev(h$term_key),
  labels = rev(h$display_wrapped)
)
h$minus_log10_q <- -log10(pmax(h$p_adjust, .Machine$double.xmin))
h$database_display <- h$database
reactome_min_q <- min(h$p_adjust[h$database == "Reactome"], na.rm = TRUE)
h$database_display[h$database %in% c("GO Biological Process", "GO_BP")] <- "GO BP"
h$database_display[h$database == "Reactome"] <- "Reactome"
h$database_display <- factor(h$database_display, levels = unique(h$database_display))
if (!is.finite(reactome_min_q) || reactome_min_q < 0.05) {
  stop("Fig.2H renderer label contract expected no Reactome term with BH q < 0.05")
}
colour_breaks <- pretty(range(h$minus_log10_q), n = 3)
colour_breaks <- colour_breaks[colour_breaks >= min(h$minus_log10_q) & colour_breaks <= max(h$minus_log10_q)]
p_h <- ggplot(h, aes(y = term_axis, x = enrichment_ratio, colour = minus_log10_q)) +
  geom_point(alpha = 0.92, size = 2.5) +
  facet_wrap(~ database_display, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_x_continuous(expand = expansion(mult = c(0.20, 0.16))) +
  scale_colour_gradient(low = "#A8DADC", high = "#9B2226", breaks = colour_breaks,
                        name = "-log10(BH q)") +
  labs(
    title = "Yellow enrichment",
    x = "Enrichment\nratio", y = NULL,
    caption = "Reactome: all BH q >= 0.05"
  ) +
  guides(
    colour = guide_colourbar(order = 1, title.position = "top", title.hjust = 0.5,
                             direction = "vertical",
                             barwidth = grid::unit(1.8, "mm"),
                             barheight = grid::unit(19, "mm"))
  ) +
  theme_panel() +
  theme(
    axis.text.y = element_text(size = 7.2, lineheight = 0.90),
    legend.position = "right",
    legend.direction = "vertical",
    legend.box = "vertical",
    legend.box.just = "center",
    legend.spacing.y = grid::unit(0.2, "mm"),
    legend.title = element_text(size = 7.2, lineheight = 0.92),
    legend.text = element_text(size = 7.2),
    legend.margin = margin(0, 0, 0, 1),
    legend.box.spacing = grid::unit(0.5, "mm"),
    strip.text = element_text(size = 7.5, face = "bold", hjust = 0),
    plot.caption = element_text(size = 7.2, colour = "black", hjust = 0),
    plot.caption.position = "plot",
    plot.margin = margin(2, 2, 1.5, 1, unit = "mm")
  )

exports <- c(
  export_pair(p_e, "Fig2E_yellow_activity_major_lineages", 70, 50),
  export_pair(p_f, "Fig2F_yellow_activity_b_states", 60, 50),
  export_pair(p_g, "Fig2G_b_state_marker_enrichment_yellow", 60, 50),
  export_pair(p_h, "Fig2H_yellow_functional_enrichment", 85, 90)
)

legend_text <- c(
  "# Proposed legends for replacement Fig.2E-H",
  "",
  "**Fig. 2E.** Signed-kME Yellow-module activity across the 13 registered major cell lineages. Raw RNA counts were aggregated by patient and cell type; only pseudobulks containing at least 20 cells were eligible. Points are independent patients and the overlaid symbols and bars denote the median and interquartile range. This is a descriptive panel; no cell-level or lineage-level inferential P value was calculated.",
  "",
  "**Fig. 2F.** Signed-kME Yellow-module activity across five provisional B-lineage transcriptional states. Raw counts were aggregated by patient and state; only pseudobulks containing at least 10 cells were eligible. Points are independent patients and summaries are median and interquartile range. The historical marker-score-assigned state labels are provisional and are not independently reference-validated.",
  "",
  "**Fig. 2G.** Enrichment of donor-aware positive state markers within the Yellow module. State-versus-rest marker models used paired patient pseudobulks and edgeR quasi-likelihood modelling. Positive markers required log2 fold change greater than zero and significance after one global Benjamini-Hochberg correction across all gene-state tests. Two-sided Fisher tests used, for each state, the intersection between tested genes and the exact WGCNA input universe; q values shown in the panel are Benjamini-Hochberg-adjusted across the five enrichment tests. Labels report Yellow overlap/positive-marker count. State labels remain provisional and circularity from historical marker-score assignment limits interpretation.",
  "",
  "**Fig. 2H.** Over-representation analysis of all 539 Yellow-module genes against the exact WGCNA input-gene universe. GO Biological Process and Reactome were analysed separately, with Benjamini-Hochberg correction within each database. Up to five terms per database are displayed, ranked by adjusted P value and then enrichment ratio; complete unfiltered results and symbol-to-ENTREZ mapping ledgers are provided as Source Data. Counts are included in the term labels and colour indicates −log10(BH q). Term names are shortened solely for display and their exact database descriptions are retained in the Source Data and display-label ledger. No displayed Reactome term passed BH q < 0.05.",
  "",
  "These panels support cross-scale descriptive coherence of the historical Yellow module; they do not validate sparse4, establish definitive B-cell-state identities or provide causal evidence."
)
writeLines(legend_text, file.path(report_dir, "FIG2E_H_PROPOSED_LEGENDS.md"), useBytes = TRUE)

write.table(
  unique(h[c("database", "ID", "Description", "display_description")]),
  file.path(report_dir, "FIG2H_DISPLAY_LABEL_LEDGER.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

manifest <- data.frame(
  file = basename(exports),
  panel = rep(c("Fig.2E", "Fig.2F", "Fig.2G", "Fig.2H"), each = 2),
  format = rep(c("PDF", "PNG"), 4),
  width_mm = rep(c(70, 60, 60, 60), each = 2),
  height_mm = rep(c(50, 50, 50, 80), each = 2),
  stringsAsFactors = FALSE
)
write.table(manifest, file.path(report_dir, "RENDER_MANIFEST.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

writeLines(
  c(
    paste("Rendered at", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    "Backend: R only",
    "PDF device: grDevices::cairo_pdf; family=Arial; one page",
    "PNG device: ragg::agg_png; 600 dpi",
    "Statistics were not recomputed by this renderer.",
    paste("Output pairs:", length(exports) / 2)
  ),
  file.path(log_dir, "render.log")
)

message("Rendered four PDF/PNG panel pairs to ", figure_dir)
