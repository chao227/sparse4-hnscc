#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1d_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1D_source.R
# ==========================================================
# fig1D_final.R — 基于 fix_fig1D_echarts.R（柱状图版）修改
# 修改项：输出路径→final/fig1/、画布55×65mm、标题9pt居中无副标题、坐标轴7pt、legend7pt、600DPI cairo
# 其他代码严格保持原脚本不变
#
# Input:  results/01_scRNAseq/annotation/cell_metadata.csv
# Output: final/fig1/fig1D_cv_lollipop.{png,pdf}

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(ragg)
})

log_msg <- function(msg) cat(format(Sys.time(),"[%H:%M:%S]"), msg, "\n", flush=TRUE)
OUT <- "final/fig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

log_msg("Building Fig1D ECharts CV bar chart...")

meta <- read.csv("results/01_scRNAseq/annotation/cell_metadata.csv")

props <- meta %>%
  group_by(patient, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(patient) %>%
  mutate(prop = n / sum(n) * 100) %>%
  ungroup()

stats <- props %>%
  group_by(cell_type) %>%
  summarise(
    mean_p = mean(prop),
    sd_p   = sd(prop),
    cv     = sd_p / mean_p * 100,
    .groups = "drop"
  ) %>%
  arrange(desc(cv)) %>%
  mutate(
    cell_type = factor(cell_type, levels = rev(cell_type)),
    highlight = cell_type %in% c("Plasma_cell", "B_cell"),
    cv_norm   = (cv - min(cv)) / (max(cv) - min(cv))
  )

label_map <- c(
  Plasma_cell    = "Plasma cell",
  B_cell         = "B cell",
  CD8_T_Cytotoxic= "CD8+ T",
  CD4_T_Helper   = "CD4+ T",
  DC_pDC         = "DC/pDC",
  Tumor_cell     = "Tumour cell",
  Mast_cell      = "Mast cell",
  Endothelial    = "Endothelial",
  Macrophage     = "Macrophage",
  Pericyte       = "Pericyte",
  CAF            = "CAF",
  Cycling        = "Cycling",
  Treg           = "Treg"
)
stats$label <- label_map[as.character(stats$cell_type)]
stats$cell_label <- factor(stats$label, levels = rev(stats$label))

pal_grad <- colorRampPalette(c("#A8DADC","#457B9D","#1D3557"))(100)
cv_max   <- max(stats$cv)

p <- ggplot(stats, aes(x = cv, y = cell_label)) +

  geom_col(aes(x = cv_max * 1.05), fill = "#F0F4F8", width = 0.68) +

  geom_col(data = filter(stats, !highlight),
           aes(x = cv, fill = cv_norm),
           width = 0.68, color = NA, alpha = 0.92) +

  geom_col(data = filter(stats, highlight),
           aes(x = cv),
           fill = "#E63946", width = 0.68, color = NA, alpha = 0.92) +

  geom_text(aes(x = cv - 3,
                label = sprintf("%.0f%%", cv),
                color = highlight),
            hjust = 1, size = 2.4, fontface = "bold", family = "Arial") +

  scale_fill_gradientn(colors = pal_grad, guide = "none") +
  scale_color_manual(values = c("FALSE" = "white", "TRUE" = "white"),
                     guide  = "none") +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.04)),
    labels = function(x) paste0(x, "%"),
    breaks = seq(0, 200, 50)
  ) +

  labs(
    title = "Inter-patient Variability of Cell-type Proportions",
    x     = "Coefficient of Variation (%)",
    y     = NULL
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    plot.title        = element_text(face = "bold", size = 9, hjust = 0.5, color = "#1D3557"),
    plot.subtitle     = element_blank(),
    axis.text.y       = element_text(size = 7,
                           face  = ifelse(levels(stats$cell_label) %in%
                                            c("Plasma cell","B cell"),
                                          "bold", "plain"),
                           color = ifelse(levels(stats$cell_label) %in%
                                            c("Plasma cell","B cell"),
                                          "#E63946", "grey25")),
    axis.text.x       = element_text(size = 7, color = "grey40"),
    axis.title.x      = element_text(size = 7),
    panel.grid.major.y = element_blank(),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_line(color = "#E8ECF0", linewidth = 0.5),
    axis.ticks         = element_blank(),
    plot.background    = element_rect(fill = "white", color = NA),
    panel.background   = element_rect(fill = "white", color = NA),
    plot.margin        = margin(10, 20, 10, 10)
  )

# ── Save: 55×65 mm ───────────────────────────────────────────────────────────
w_mm <- 65; h_mm <- 55
w <- w_mm / 25.4; h <- h_mm / 25.4
ragg::agg_png(file.path(OUT, "fig1D_cv_lollipop.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p); dev.off()
grDevices::cairo_pdf(file.path(OUT, "fig1D_cv_lollipop.pdf"),
                     width = w, height = h, family = "Arial")
print(p); dev.off()
log_msg("Saved: final/fig1/fig1D_cv_lollipop.{png,pdf} [55×65mm]")
log_msg("=== fig1D_final.R complete ===")
