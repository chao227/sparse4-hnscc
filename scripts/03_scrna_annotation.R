#!/usr/bin/env Rscript

# ==========================================================
# Script   : 03_scrna_annotation.R
# Purpose  : Generate broad single-cell annotations
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: cell annotation; Fig. 1
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/01_single_cell/01d_annotation.R
# ==========================================================
# =============================================================================
# 脚本: 01d_annotation.R
# 功能: 细胞类型注释（标记基因 + FindAllMarkers）
# 输入: results/01_scRNAseq/clustering/seurat_clustered.rds
# 输出: results/01_scRNAseq/annotation/seurat_annotated.rds
#        figures/module01/dotplot_markers.pdf
#        figures/module01/umap_celltype.pdf
# 预计内存: ~60 GB
# 预计时间: 1-2小时（FindAllMarkers）
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(future)
})

plan("multicore", workers=16)  # 并行加速FindAllMarkers
options(future.globals.maxSize = 80 * 1024^3)

# ── 路径设置 ──────────────────────────────────────────────────────────────────
base_path  <- "."
clust_out  <- file.path(base_path, "results/01_scRNAseq/clustering")
annot_out  <- file.path(base_path, "results/01_scRNAseq/annotation")
fig_out    <- file.path(base_path, "figures/module01")
log_file   <- file.path(base_path, "logs/01d_annotation.log")

dir.create(annot_out, recursive=TRUE, showWarnings=FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file=log_file, append=TRUE)
}

log_msg("=== Module 1d: Cell Type Annotation 开始 ===")

# ── 检查上游输出 ──────────────────────────────────────────────────────────────
input_file <- file.path(clust_out, "seurat_clustered.rds")
if (!file.exists(input_file)) stop("缺少输入文件: ", input_file)

log_msg("加载 seurat_clustered.rds ...")
seu <- readRDS(input_file)
log_msg(paste("细胞数:", ncol(seu), " | 聚类数:", length(unique(Idents(seu)))))

# ── HNSCC TME 标记基因定义 ────────────────────────────────────────────────────
cell_markers <- list(
  # ── 免疫细胞 (CD45p) ────────────────────────────────────────────────────
  "CD8_T_Cytotoxic"  = c("CD8A","CD8B","GZMB","GZMK","PRF1","NKG7"),
  "CD8_T_Exhausted"  = c("CD8A","PDCD1","LAG3","HAVCR2","TIGIT","TOX"),
  "CD4_T_Helper"     = c("CD4","IL7R","CCR7","TCF7","LEF1"),
  "Treg"             = c("FOXP3","IL2RA","CTLA4","IKZF2","TNFRSF9"),
  "NK"               = c("GNLY","NKG7","KLRD1","NCAM1","FCGR3A"),
  "B_cell"           = c("CD79A","CD79B","MS4A1","CD19","PAX5"),
  "Plasma_cell"      = c("IGHG1","JCHAIN","MZB1","SDC1","PRDM1"),
  "Macrophage"       = c("CD68","CD163","MRC1","FCGR3A","CSF1R"),
  "Macrophage_M2"    = c("CD163","MRC1","TGFB1","IL10","MARCO"),
  "DC_cDC"           = c("ITGAX","HLA-DRA","CD1C","CLEC10A","FLT3"),
  "DC_pDC"           = c("LILRA4","CLEC4C","IL3RA","TCF4"),
  "Mast_cell"        = c("TPSAB1","CPA3","KIT","MS4A2"),
  # ── 非免疫细胞 (CD45n) ──────────────────────────────────────────────────
  "Tumor_cell"       = c("EPCAM","KRT5","KRT14","KRT17","TP63","DSP"),
  "CAF"              = c("FAP","ACTA2","COL1A1","COL1A2","PDPN","S100A4"),
  "Endothelial"      = c("PECAM1","VWF","CDH5","CLDN5","ESAM"),
  "Pericyte"         = c("PDGFRB","RGS5","ACTA2","MCAM","NOTCH3"),
  # ── 通用 ─────────────────────────────────────────────────────────────────
  "Cycling"          = c("MKI67","TOP2A","PCNA","CDK1","STMN1")
)

# 平铺所有标记基因（用于 DotPlot）
all_markers_flat <- unique(unlist(cell_markers))

# ── Step 1: 计算模块分数 ──────────────────────────────────────────────────────
log_msg("Step 1: 计算各细胞类型模块分数 ...")

for (ct in names(cell_markers)) {
  genes_use <- intersect(cell_markers[[ct]], rownames(seu))
  if (length(genes_use) >= 3) {
    seu <- AddModuleScore(seu,
      features = list(genes_use),
      name     = paste0("score_", ct),
      ctrl     = 50,
      seed     = 42
    )
    # AddModuleScore 会在列名后加1，重命名
    old_col <- paste0("score_", ct, "1")
    new_col <- paste0("score_", ct)
    colnames(seu@meta.data)[colnames(seu@meta.data) == old_col] <- new_col
  }
}
log_msg("模块分数计算完成")

# ── Step 2: FindAllMarkers（差异基因）────────────────────────────────────────
log_msg("Step 2: FindAllMarkers (Wilcoxon test) ...")
all_deg <- FindAllMarkers(
  seu,
  only.pos    = TRUE,
  min.pct     = 0.25,
  logfc.threshold = 0.5,
  test.use    = "wilcox",
  verbose     = FALSE
)
log_msg(paste("差异基因总数:", nrow(all_deg)))

write.csv(all_deg, file.path(annot_out, "all_cluster_markers.csv"), row.names=FALSE)

# 每个cluster top markers
top_markers <- all_deg %>%
  group_by(cluster) %>%
  slice_max(order_by=avg_log2FC, n=10) %>%
  ungroup()
write.csv(top_markers, file.path(annot_out, "top10_markers_per_cluster.csv"), row.names=FALSE)

# ── Step 3: 自动注释（基于标记基因重叠 + 模块分数）──────────────────────────
log_msg("Step 3: 自动细胞类型注释 ...")

# 方法：对每个cluster，找与标记基因重叠最多的细胞类型
clusters <- levels(Idents(seu))
cluster_annot <- data.frame(
  cluster   = clusters,
  cell_type = NA_character_,
  score     = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(clusters)) {
  cl <- clusters[i]
  cl_markers <- all_deg %>%
    filter(cluster == cl, p_val_adj < 0.05) %>%
    pull(gene)

  if (length(cl_markers) == 0) {
    cluster_annot$cell_type[i] <- "Unknown"
    cluster_annot$score[i]     <- 0
    next
  }

  # 与每个细胞类型标记基因的重叠度 (Jaccard-like)
  overlaps <- sapply(cell_markers, function(markers) {
    length(intersect(cl_markers, markers)) / length(markers)
  })

  best_ct    <- names(which.max(overlaps))
  best_score <- max(overlaps)

  cluster_annot$cell_type[i] <- best_ct
  cluster_annot$score[i]     <- best_score
}

# 额外规则：用 fraction 信息辅助注释
# CD45n fraction 中不应有纯免疫细胞类型
immune_types <- c("CD8_T_Cytotoxic","CD8_T_Exhausted","CD4_T_Helper","Treg",
                  "NK","B_cell","Plasma_cell","Macrophage","Macrophage_M2",
                  "DC_cDC","DC_pDC","Mast_cell")
nonimmune_types <- c("Tumor_cell","CAF","Endothelial","Pericyte")

# 按 fraction 调整：计算每个 cluster 中 CD45p/CD45n 细胞比例
cluster_frac <- table(Idents(seu), seu$fraction)
cluster_frac_pct <- prop.table(cluster_frac, margin=1)

for (i in seq_along(clusters)) {
  cl <- clusters[i]
  ct <- cluster_annot$cell_type[i]

  if (cl %in% rownames(cluster_frac_pct)) {
    pct_cd45n <- cluster_frac_pct[cl, "CD45n"]
    # 如果>60%为CD45n细胞但注释为免疫类型，重新评估
    if (!is.na(pct_cd45n) && pct_cd45n > 0.6 && ct %in% immune_types) {
      # 找最高的非免疫类型分数
      cl_markers_v <- all_deg %>%
        filter(cluster == cl, p_val_adj < 0.05) %>%
        pull(gene)
      ni_overlaps <- sapply(cell_markers[nonimmune_types], function(m)
        length(intersect(cl_markers_v, m)) / max(length(m), 1))
      if (max(ni_overlaps) > 0) {
        cluster_annot$cell_type[i] <- names(which.max(ni_overlaps))
      }
    }
  }
}

write.csv(cluster_annot, file.path(annot_out, "cluster_annotation.csv"), row.names=FALSE)
log_msg("聚类注释完成")
log_msg(paste("注释细胞类型:", paste(unique(cluster_annot$cell_type), collapse=", ")))

# ── Step 4: 添加注释到 Seurat 对象 ───────────────────────────────────────────
annot_map <- setNames(cluster_annot$cell_type, cluster_annot$cluster)
# Seurat v5: 赋值时去除向量名字，避免 "No cell overlap" 错误
seu$cell_type <- unname(annot_map[as.character(Idents(seu))])

# 定义粗细两级注释
coarse_map <- c(
  "CD8_T_Cytotoxic"  = "T_cell",
  "CD8_T_Exhausted"  = "T_cell",
  "CD4_T_Helper"     = "T_cell",
  "Treg"             = "T_cell",
  "NK"               = "NK_cell",
  "B_cell"           = "B_cell",
  "Plasma_cell"      = "Plasma_cell",
  "Macrophage"       = "Myeloid",
  "Macrophage_M2"    = "Myeloid",
  "DC_cDC"           = "Myeloid",
  "DC_pDC"           = "Myeloid",
  "Mast_cell"        = "Mast_cell",
  "Tumor_cell"       = "Tumor_cell",
  "CAF"              = "CAF",
  "Endothelial"      = "Endothelial",
  "Pericyte"         = "Pericyte",
  "Cycling"          = "Cycling",
  "Unknown"          = "Unknown"
)
seu$cell_type_coarse <- unname(coarse_map[seu$cell_type])

# ── Step 5: 可视化 ────────────────────────────────────────────────────────────
log_msg("Step 5: 绘制注释图 ...")

# 细胞类型调色板
ct_colors <- c(
  "CD8_T_Cytotoxic" = "#D62728","CD8_T_Exhausted" = "#FF7F0E",
  "CD4_T_Helper"    = "#AEC7E8","Treg"             = "#FFBB78",
  "NK"              = "#98DF8A","B_cell"           = "#1F77B4",
  "Plasma_cell"     = "#9467BD","Macrophage"       = "#8C564B",
  "Macrophage_M2"   = "#C49C94","DC_cDC"           = "#E377C2",
  "DC_pDC"          = "#F7B6D2","Mast_cell"        = "#7F7F7F",
  "Tumor_cell"      = "#BCBD22","CAF"              = "#17BECF",
  "Endothelial"     = "#9EDAE5","Pericyte"         = "#393B79",
  "Cycling"         = "#FFDD89","Unknown"          = "#CCCCCC"
)

# UMAP 细胞类型
p_annot <- DimPlot(seu, reduction="umap", group.by="cell_type",
                    label=TRUE, label.size=3.5, pt.size=0.1,
                    cols=ct_colors[unique(seu$cell_type)]) +
  labs(title="UMAP - Cell Type Annotation") +
  theme_bw() + theme(legend.position="right")
ggsave(file.path(fig_out, "umap_celltype.pdf"), p_annot, width=14, height=10)
ggsave(file.path(fig_out, "umap_celltype.png"), p_annot, width=14, height=10, dpi=150)
log_msg("  图: umap_celltype — UMAP细胞类型注释图，各细胞类型按颜色区分")

# 粗分类 UMAP
coarse_colors <- c(
  "T_cell"="#D62728","NK_cell"="#98DF8A","B_cell"="#1F77B4",
  "Plasma_cell"="#9467BD","Myeloid"="#8C564B","Mast_cell"="#7F7F7F",
  "Tumor_cell"="#BCBD22","CAF"="#17BECF","Endothelial"="#9EDAE5",
  "Pericyte"="#393B79","Cycling"="#FFDD89","Unknown"="#CCCCCC"
)
p_coarse <- DimPlot(seu, reduction="umap", group.by="cell_type_coarse",
                     label=TRUE, label.size=4, pt.size=0.1,
                     cols=coarse_colors[unique(seu$cell_type_coarse)]) +
  labs(title="UMAP - Broad Cell Type Annotation") +
  theme_bw()
ggsave(file.path(fig_out, "umap_celltype_coarse.pdf"), p_coarse, width=12, height=10)
ggsave(file.path(fig_out, "umap_celltype_coarse.png"), p_coarse, width=12, height=10, dpi=150)
log_msg("  图: umap_celltype_coarse — UMAP粗分类注释图（T/NK/B/Myeloid/Tumor/CAF/Endo等）")

# DotPlot 标记基因
markers_to_plot <- c(
  "CD8A","GZMB","PRF1","PDCD1","LAG3","CD4","FOXP3","GNLY","NKG7",
  "CD79A","MS4A1","JCHAIN","MZB1","CD68","CD163","CLEC9A","LILRA4",
  "TPSAB1","EPCAM","KRT14","FAP","COL1A1","PECAM1","VWF","MKI67"
)
markers_in_data <- intersect(markers_to_plot, rownames(seu))

p_dot <- DotPlot(seu, features=markers_in_data, group.by="cell_type",
                  dot.scale=6, col.min=-2, col.max=2) +
  scale_color_gradient2(low="#2166AC", mid="white", high="#D6604D") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=9),
        axis.text.y=element_text(size=9)) +
  labs(title="Canonical Marker Gene Expression by Cell Type",
       x="Gene", y="Cell Type")
ggsave(file.path(fig_out, "dotplot_markers.pdf"), p_dot, width=18, height=10)
ggsave(file.path(fig_out, "dotplot_markers.png"), p_dot, width=18, height=10, dpi=150)
log_msg("  图: dotplot_markers — 经典标记基因在各细胞类型中的表达DotPlot")

# FeaturePlot 关键标记
key_genes <- c("CD8A","FOXP3","CD68","EPCAM","FAP","PECAM1","MKI67","PDCD1")
key_genes <- intersect(key_genes, rownames(seu))
for (g in key_genes) {
  p_feat <- FeaturePlot(seu, features=g, reduction="umap", pt.size=0.1,
                         order=TRUE) +
    scale_color_gradient(low="lightgrey", high="#D62728") +
    theme_bw() +
    labs(title=paste("Expression:", g))
  ggsave(file.path(fig_out, paste0("featureplot_", g, ".pdf")), p_feat,
         width=8, height=7)
  ggsave(file.path(fig_out, paste0("featureplot_", g, ".png")), p_feat,
         width=8, height=7, dpi=150)
}
log_msg(paste("  图: featureplot_*.pdf/png — 关键基因在UMAP上的表达分布（",
              paste(key_genes, collapse="/"), "）"))

# 细胞类型比例图（按患者）
celltype_prop <- as.data.frame(table(seu$patient, seu$cell_type_coarse))
colnames(celltype_prop) <- c("Patient","CellType","Count")
celltype_prop <- celltype_prop %>%
  group_by(Patient) %>%
  mutate(Pct = Count / sum(Count) * 100)

p_prop <- ggplot(celltype_prop, aes(x=Patient, y=Pct, fill=CellType)) +
  geom_col(position="stack") +
  scale_fill_manual(values=coarse_colors) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1)) +
  labs(title="Cell Type Composition per Patient",
       x="Patient", y="Proportion (%)")
ggsave(file.path(fig_out, "celltype_composition.pdf"), p_prop, width=14, height=7)
ggsave(file.path(fig_out, "celltype_composition.png"), p_prop, width=14, height=7, dpi=150)
log_msg("  图: celltype_composition — 各患者细胞类型比例堆积柱状图")

# ── Step 6: 保存 ──────────────────────────────────────────────────────────────
log_msg("Step 6: 保存注释后 Seurat 对象 ...")
saveRDS(seu, file.path(annot_out, "seurat_annotated.rds"))

# 导出 meta.data 便于后续分析
write.csv(seu@meta.data, file.path(annot_out, "cell_metadata.csv"), row.names=TRUE)

# 细胞类型数量统计
ct_counts <- as.data.frame(table(seu$cell_type))
colnames(ct_counts) <- c("cell_type","n_cells")
ct_counts$pct <- round(ct_counts$n_cells / sum(ct_counts$n_cells) * 100, 2)
write.csv(ct_counts, file.path(annot_out, "celltype_counts.csv"), row.names=FALSE)

# ── 自检 ─────────────────────────────────────────────────────────────────────
log_msg("=== 自检 Module 1d ===")
n_types <- length(unique(seu$cell_type[seu$cell_type != "Unknown"]))
log_msg(paste("识别的细胞类型数:", n_types, "(预期: ≥8)"))
if (n_types < 8) log_msg("⚠ 识别细胞类型<8，注释可能不完整")

required_types <- c("Tumor_cell","CAF","Endothelial")
cd45n_cells <- seu@meta.data[seu$fraction == "CD45n", ]
has_nonimmune <- any(required_types %in% unique(cd45n_cells$cell_type))
if (!has_nonimmune) {
  log_msg("⚠ CD45n中未找到Tumor/CAF/Endothelial，需要检查注释")
} else {
  log_msg("✓ CD45n中存在非免疫细胞类型")
}

log_msg("=== Module 1d 完成 ===")
log_msg(paste("输出: results/01_scRNAseq/annotation/seurat_annotated.rds"))
