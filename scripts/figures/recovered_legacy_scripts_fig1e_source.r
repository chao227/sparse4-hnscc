#!/usr/bin/env Rscript

# ==========================================================
# Script   : recovered_legacy_scripts_fig1e_source.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:final_audit_reconciliation_20261022/legacy_panel_recovery/recovered_legacy_scripts/Fig1E_source.R
# ==========================================================
# fig1E_final.R — 基于 fig1e_bcell_umap_v2.R 修改
# 修改项：输出路径→final/fig1/、画布55×65mm、字体统一（标题9pt居中无副标题、坐标轴7-8pt、legend7pt）
# 其他代码严格保持原脚本不变
#
# Input:  results/subcluster/b_subclusters.rds
# Output: final/fig1/fig1E_bcell_subtype_umap.{png,pdf}

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(ggrepel); library(ragg); library(ggrastr)
})

log_msg <- function(x) cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), x))
OUT <- "final/fig1"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

B_COLORS <- c(
  "Naive B"        = "#1560A8",
  "Transitional B" = "#4D9DD4",
  "GC-like B"      = "#D97E12",
  "Memory B"       = "#38A038",
  "Plasma cell"    = "#C8351A"
)

log_msg("Loading B-cell subclusters...")
b_seu <- readRDS("results/subcluster/b_subclusters.rds")
b_seu$b_subtype <- factor(b_seu$b_subtype, levels = names(B_COLORS))
n <- ncol(b_seu)
log_msg(sprintf("  %d B cells", n))

pt <- if (n < 30000) 0.5 else if (n < 100000) 0.35 else 0.22
al <- if (n < 30000) 0.90 else if (n < 100000) 0.82 else 0.75

b_umap <- as.data.frame(Embeddings(b_seu, "umap"))
b_umap$b_subtype <- b_seu$b_subtype
b_centroids <- b_umap %>%
  group_by(b_subtype) %>%
  summarise(U1 = median(umap_1), U2 = median(umap_2), .groups = "drop")

p1E <- ggplot(b_umap, aes(umap_1, umap_2, colour = b_subtype)) +
  geom_point(size = pt, alpha = al, stroke = 0) +
  ggrepel::geom_text_repel(
    data = b_centroids,
    aes(U1, U2, label = b_subtype),
    colour = "black",
    size = 2.2, fontface = "bold", family = "Arial",
    segment.size = 0.25, segment.colour = "grey50",
    min.segment.length = 0.2, box.padding = 0.3,
    force = 1.5, inherit.aes = FALSE
  ) +
  scale_colour_manual(values = B_COLORS, name = "B-cell subtype",
                      guide = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  labs(
    title = "B-cell lineage subclustering (5 maturation states)",
    x = "UMAP 1", y = "UMAP 2"
  ) +
  theme_classic(base_size = 7, base_family = "Arial") +
  theme(
    axis.line      = element_line(linewidth = 0.35),
    axis.ticks     = element_blank(),
    axis.text      = element_blank(),
    axis.title     = element_text(size = 7),
    legend.position = "right",
    legend.title   = element_text(size = 7),
    legend.text    = element_text(size = 7),
    legend.key.size = unit(3, "mm"),
    plot.title     = element_text(size = 9, face = "bold", hjust = 0.5),
    plot.subtitle  = element_blank(),
    panel.grid     = element_blank()
  )
p1E$layers[[1]]$aes_params$alpha <- al

# 光栅化点云层
p1E_rast <- ggrastr::rasterize(p1E, layers = "Point", dpi = 600)

# ── Save: 55×65 mm ───────────────────────────────────────────────────────────
w_mm <- 65; h_mm <- 55
w <- w_mm / 25.4; h <- h_mm / 25.4
ragg::agg_png(file.path(OUT, "fig1E_bcell_subtype_umap.png"),
              width = w, height = h, units = "in", res = 600, bg = "white")
print(p1E_rast); dev.off()
grDevices::cairo_pdf(file.path(OUT, "fig1E_bcell_subtype_umap.pdf"), width = w, height = h, family = "Arial")
print(p1E_rast); dev.off()
log_msg("Saved: final/fig1/fig1E_bcell_subtype_umap.{png,pdf} [55×65mm]")
log_msg("=== fig1E_final.R complete ===")
