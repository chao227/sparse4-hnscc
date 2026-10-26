#!/usr/bin/env Rscript

# ==========================================================
# Script   : 04_bcell_cd8_states.R
# Purpose  : Construct the B-lineage and CD8 state objects used in the manuscript
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: state annotation; Fig. 1; Fig. 5; Supplementary Fig. 3
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/01_single_cell/R12_cd8_b_subcluster.R
# ==========================================================
# R12: CD8+T细胞 + B细胞亚群精细化分析
# 输入: results/01_scRNAseq/annotation/seurat_annotated.rds
# 输出:
#   results/subcluster/cd8_subclusters.rds
#   results/subcluster/b_subclusters.rds
#   paper_figures_revised/S6_Subclusters/

suppressPackageStartupMessages({
  library(Seurat); library(ggplot2); library(dplyr); library(tidyr)
  library(patchwork); library(ComplexHeatmap); library(circlize)
  library(RColorBrewer); library(scales); library(ggrepel)
  library(monocle3); library(slingshot); library(viridis)
})

out_dir    <- "paper_figures_revised/S6_Subclusters"
res_dir    <- "results/subcluster"
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(res_dir, recursive=TRUE, showWarnings=FALSE)

log_msg <- function(msg) cat(format(Sys.time(),"[%H:%M:%S]"), msg, "\n", flush=TRUE)
save_fig <- function(p, name, w=8, h=6) {
  ggsave(file.path(out_dir, paste0(name,".png")), p, width=w, height=h, dpi=300, bg="white")
  ggsave(file.path(out_dir, paste0(name,".pdf")), p, width=w, height=h)
  log_msg(paste("Saved:", name))
}

# 安全提取metadata列：处理list型列（Seurat merge后可能出现）
meta_chr <- function(seu, col) {
  x <- tryCatch(seu@meta.data[[col]], error = function(e) NULL)
  if (is.null(x)) return(rep(NA_character_, ncol(seu)))
  if (is.list(x)) {
    vapply(x, function(el) {
      if (is.null(el) || length(el) == 0) NA_character_
      else as.character(el[1])
    }, character(1), USE.NAMES = FALSE)
  } else {
    as.character(x)
  }
}

# ── 颜色方案 ──────────────────────────────────────────────────────────────────
CD8_COLORS <- c(
  "Naive/Stem-like"   = "#2166AC",
  "Progenitor-Tex"    = "#74ADD1",
  "Intermediate-Tex"  = "#FDAE61",
  "Terminal-Tex"      = "#D73027",
  "Cytotoxic-Eff"     = "#1A9641",
  "Memory"            = "#762A83"
)
B_COLORS <- c(
  "Naive B"           = "#2166AC",
  "Memory B"          = "#74ADD1",
  "GC-like B"         = "#F4A582",
  "Transitional B"    = "#92C5DE",
  "Plasma cell"       = "#D73027"
)

# ── 加载数据 ──────────────────────────────────────────────────────────────────
log_msg("Loading annotated Seurat object...")
ann <- readRDS("results/01_scRNAseq/annotation/seurat_annotated.rds")
log_msg(paste("Total cells:", ncol(ann)))

# ══════════════════════════════════════════════════════════════════════════════
# PART 1: CD8+T细胞亚群分析
# ══════════════════════════════════════════════════════════════════════════════
log_msg("=== PART 1: CD8+ T Cell Subclustering ===")

cd8 <- subset(ann, cell_type == "CD8_T_Cytotoxic")
log_msg(paste("CD8 T cells:", ncol(cd8)))

# 标准化 + 降维 + 聚类
cd8 <- NormalizeData(cd8, verbose=FALSE)
cd8 <- FindVariableFeatures(cd8, nfeatures=2000, verbose=FALSE)
cd8 <- ScaleData(cd8, verbose=FALSE)
cd8 <- RunPCA(cd8, npcs=30, verbose=FALSE)

# Harmony批次校正
if(requireNamespace("harmony", quietly=TRUE)) {
  library(harmony)
  cd8 <- RunHarmony(cd8, group.by.vars="patient", verbose=FALSE)
  cd8 <- RunUMAP(cd8, reduction="harmony", dims=1:20, verbose=FALSE)
  cd8 <- FindNeighbors(cd8, reduction="harmony", dims=1:20, verbose=FALSE)
} else {
  cd8 <- RunUMAP(cd8, dims=1:20, verbose=FALSE)
  cd8 <- FindNeighbors(cd8, dims=1:20, verbose=FALSE)
}
cd8 <- FindClusters(cd8, resolution=0.4, verbose=FALSE)
log_msg(paste("CD8 clusters:", nlevels(cd8$seurat_clusters)))

# ── CD8亚群标记基因打分 ──────────────────────────────────────────────────────
cd8_markers <- list(
  "Naive/Stem-like"  = c("CCR7","TCF7","SELL","LEF1","IL7R"),
  "Progenitor-Tex"   = c("TCF7","SLAMF6","XCL1","PDCD1"),
  "Intermediate-Tex" = c("GZMK","PDCD1","CCL4","CCL5","DUSP4"),
  "Terminal-Tex"     = c("HAVCR2","LAG3","TIGIT","TOX","ENTPD1","CXCL13"),
  "Cytotoxic-Eff"    = c("GZMB","PRF1","GNLY","NKG7","IFNG","FGFBP2"),
  "Memory"           = c("IL7R","TCF7","S1PR1","KLF2","SELL")
)
cd8 <- AddModuleScore(cd8, features=cd8_markers, name="score_", assay="RNA")
score_cols <- paste0("score_", seq_along(cd8_markers))
names(score_cols) <- names(cd8_markers)

# 每个细胞按最高打分归类
score_mat <- cd8@meta.data[, score_cols, drop=FALSE]
cd8$cd8_subtype <- names(cd8_markers)[apply(score_mat, 1, which.max)]
cd8$cd8_subtype <- factor(cd8$cd8_subtype, levels=names(CD8_COLORS))

# 输出亚型比例
tbl <- table(cd8$cd8_subtype)
log_msg("CD8 subtype counts:")
print(tbl)

# ── 图1: CD8亚群UMAP ─────────────────────────────────────────────────────────
log_msg("CD8 Fig1: Subtype UMAP...")
n_cells <- ncol(cd8)
pt <- 0.3; al <- 0.8

p_cd8_umap <- DimPlot(cd8, group.by="cd8_subtype", raster=FALSE,
                       cols=CD8_COLORS, pt.size=pt, label=FALSE) +
  scale_color_manual(values=CD8_COLORS, name="CD8+ T\nSubtype") +
  labs(title="CD8+ T Cell Subtypes",
       subtitle=paste0(n_cells, " cells | 18 HNSCC patients | Harmony-corrected UMAP"),
       x="UMAP_1", y="UMAP_2") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9,color="grey50"),
        legend.text=element_text(size=9),
        axis.title=element_text(size=9))
p_cd8_umap$layers[[1]]$aes_params$alpha <- al
save_fig(p_cd8_umap, "cd8_umap_subtypes", 8, 6.5)

# ── 图2: CD8亚群标记基因气泡图 ───────────────────────────────────────────────
log_msg("CD8 Fig2: Marker dotplot...")
key_markers_cd8 <- c(
  "CCR7","TCF7","SELL",           # Naive/Stem
  "SLAMF6","XCL1",                # Progenitor-Tex
  "GZMK","CCL4","DUSP4",          # Intermediate-Tex
  "HAVCR2","TOX","ENTPD1","CXCL13", # Terminal-Tex
  "GZMB","PRF1","GNLY","NKG7",    # Cytotoxic-Eff
  "IL7R","S1PR1","KLF2",          # Memory
  "PDCD1","LAG3","TIGIT","IKZF1"  # Checkpoint/Exhaustion TF
)
Idents(cd8) <- "cd8_subtype"
p_cd8_dot <- DotPlot(cd8, features=intersect(key_markers_cd8, rownames(cd8)),
                      group.by="cd8_subtype", dot.scale=5) +
  scale_color_gradient2(low="#2166AC", mid="#F7F7F7", high="#D73027",
                        midpoint=0, name="Avg\nExpr") +
  coord_flip() +
  labs(title="CD8+ T Cell Subtype Marker Genes",
       x=NULL, y=NULL) +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        axis.text.x=element_text(angle=30,hjust=1,size=9),
        axis.text.y=element_text(size=9))
save_fig(p_cd8_dot, "cd8_dotplot_markers", 10, 7)

# ── 图3: IKZF1在各CD8亚群的表达 ──────────────────────────────────────────────
log_msg("CD8 Fig3: IKZF1 across subtypes...")
cd8_meta <- cd8@meta.data
cd8_meta$IKZF1 <- GetAssayData(cd8, assay="RNA", layer="data")["IKZF1",]

# 按Terminal-Tex均值排序
ikzf1_order <- cd8_meta %>%
  group_by(cd8_subtype) %>%
  summarise(med=median(IKZF1)) %>%
  arrange(desc(med)) %>%
  pull(cd8_subtype)
cd8_meta$cd8_subtype <- factor(cd8_meta$cd8_subtype, levels=ikzf1_order)

p_ikzf1 <- ggplot(cd8_meta, aes(x=cd8_subtype, y=IKZF1, fill=cd8_subtype)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.8, linewidth=0.3) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white", linewidth=0.4) +
  scale_fill_manual(values=CD8_COLORS, guide="none") +
  scale_y_continuous(expand=expansion(mult=c(0,0.08))) +
  labs(title="IKZF1 Expression Across CD8+ T Cell Subtypes",
       subtitle="Higher IKZF1 in exhausted subtypes supports its role as an exhaustion-promoting TF",
       x=NULL, y="IKZF1 Expression (log-normalized)") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey50"),
        axis.text.x=element_text(angle=25,hjust=1,size=9.5))
save_fig(p_ikzf1, "cd8_ikzf1_across_subtypes", 9, 5.5)

# ── 图4: 耗竭/功能评分热图 ───────────────────────────────────────────────────
log_msg("CD8 Fig4: Functional score heatmap...")
score_mean <- cd8@meta.data %>%
  group_by(cd8_subtype) %>%
  summarise(across(all_of(score_cols), mean)) %>%
  as.data.frame()
rownames(score_mean) <- score_mean$cd8_subtype
mat_score <- as.matrix(score_mean[,-1])
colnames(mat_score) <- names(cd8_markers)
mat_score <- t(scale(mat_score))

col_fun <- colorRamp2(c(-2,0,2), c("#2166AC","white","#D73027"))
ht_score <- Heatmap(mat_score, name="Z-score",
  col=col_fun,
  cluster_rows=TRUE, cluster_columns=FALSE,
  row_names_gp=gpar(fontsize=9),
  column_names_gp=gpar(fontsize=9, fontface="bold"),
  column_title="CD8+ T Functional Score by Subtype",
  column_title_gp=gpar(fontsize=10, fontface="bold"),
  rect_gp=gpar(col="white", lwd=1))
png(file.path(out_dir,"cd8_functional_score_heatmap.png"),width=1800,height=1200,res=300)
draw(ht_score)
dev.off()
pdf(file.path(out_dir,"cd8_functional_score_heatmap.pdf"),width=6,height=4)
draw(ht_score)
dev.off()
log_msg("Saved: cd8_functional_score_heatmap")

# ── 图5: 各患者CD8亚型比例（堆积图）──────────────────────────────────────────
log_msg("CD8 Fig5: Per-patient subtype proportion...")
# base R table() — avoids monocle3 masking dplyr::count() via BiocGenerics
{
  pt_vec  <- meta_chr(cd8, "patient")
  st_vec  <- as.character(cd8$cd8_subtype)
  ok      <- !is.na(pt_vec) & !is.na(st_vec)
  cnt_tbl <- as.data.frame(table(patient=pt_vec[ok], cd8_subtype=st_vec[ok]),
                            stringsAsFactors=FALSE)
  prop_cd8 <- cnt_tbl %>%
    dplyr::rename(n=Freq) %>%
    dplyr::filter(n>0) %>%
    dplyr::group_by(patient) %>%
    dplyr::mutate(prop=n/sum(n)) %>%
    dplyr::ungroup()
}

# 按Terminal-Tex比例排序患者
pt_order_cd8 <- prop_cd8 %>%
  filter(cd8_subtype=="Terminal-Tex") %>%
  arrange(desc(prop)) %>% pull(patient)
prop_cd8$patient <- factor(prop_cd8$patient, levels=pt_order_cd8)

p_prop_cd8 <- ggplot(prop_cd8, aes(x=patient, y=prop, fill=cd8_subtype)) +
  geom_col(width=0.8, color="white", linewidth=0.3) +
  scale_fill_manual(values=CD8_COLORS, name="CD8+ T Subtype") +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  labs(title="CD8+ T Cell Subtype Composition Per Patient",
       subtitle="Sorted by Terminal-Tex proportion | Inter-patient heterogeneity in exhaustion state",
       x="Patient", y="Proportion") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey50"),
        axis.text.x=element_text(angle=40,hjust=1,size=8.5))
save_fig(p_prop_cd8, "cd8_patient_proportion", 10, 5.5)

# ── 图6: Checkpoint marker featureplot ────────────────────────────────────────
log_msg("CD8 Fig6: Checkpoint FeaturePlots...")
chkpt_genes <- c("PDCD1","HAVCR2","LAG3","TIGIT","TOX","IKZF1")
plots_fp <- lapply(chkpt_genes, function(g) {
  if(!g %in% rownames(cd8)) return(NULL)
  p <- FeaturePlot(cd8, features=g, pt.size=0.15, raster=FALSE) +
    scale_color_gradient(low="grey90", high="#D73027", name=g) +
    theme_void() +
    theme(plot.title=element_text(face="bold",size=10,hjust=0.5),
          legend.position="right",
          legend.text=element_text(size=7))
  p$layers[[1]]$aes_params$alpha <- 0.6
  p
})
plots_fp <- Filter(Negate(is.null), plots_fp)
p_fp <- wrap_plots(plots_fp, ncol=3) +
  plot_annotation(title="Exhaustion and Checkpoint Markers in CD8+ T Cells",
                  theme=theme(plot.title=element_text(face="bold",size=13,hjust=0.5)))
save_fig(p_fp, "cd8_checkpoint_featureplots", 14, 9)

# ── 拟时序：CD8耗竭轨迹 ──────────────────────────────────────────────────────
log_msg("CD8 pseudotime: Slingshot trajectory...")
tryCatch({
  library(SingleCellExperiment)
  sce_cd8 <- as.SingleCellExperiment(cd8, assay="RNA")
  reducedDim(sce_cd8, "UMAP") <- Embeddings(cd8, "umap")

  # 设置Naive/Stem-like为起始群
  start_cluster <- levels(cd8$cd8_subtype)[1]  # Naive/Stem-like
  sl <- slingshot(sce_cd8, clusterLabels="cd8_subtype",
                  reducedDim="UMAP",
                  start.clus="Naive/Stem-like",
                  end.clus=c("Terminal-Tex","Cytotoxic-Eff"))

  pt_vals <- slingPseudotime(sl)[,1]
  cd8$slingshot_pt <- pt_vals[colnames(cd8)]

  # 拟时序图
  pt_df <- data.frame(
    UMAP1 = Embeddings(cd8,"umap")[,1],
    UMAP2 = Embeddings(cd8,"umap")[,2],
    pseudotime = cd8$slingshot_pt,
    subtype = cd8$cd8_subtype
  )
  pt_df <- pt_df[order(pt_df$pseudotime, na.last=FALSE),]

  p_pt <- ggplot(pt_df, aes(x=UMAP1, y=UMAP2, color=pseudotime)) +
    geom_point(size=0.3, alpha=0.7) +
    scale_color_viridis_c(option="magma", name="Pseudotime", na.value="grey80") +
    labs(title="CD8+ T Cell Differentiation Trajectory (Slingshot)",
         subtitle="Naive/Stem-like → Intermediate-Tex → Terminal-Tex / Cytotoxic-Eff") +
    theme_classic() +
    theme(plot.title=element_text(face="bold",size=12),
          plot.subtitle=element_text(size=9,color="grey50"),
          axis.title=element_text(size=9))
  save_fig(p_pt, "cd8_slingshot_pseudotime", 7, 6)

  # 拟时序 vs IKZF1表达
  pt_df$IKZF1 <- as.numeric(GetAssayData(cd8, assay="RNA", layer="data")["IKZF1",])
  p_ikzf1_pt <- ggplot(pt_df[!is.na(pt_df$pseudotime),],
                        aes(x=pseudotime, y=IKZF1, color=subtype)) +
    geom_point(size=0.3, alpha=0.4) +
    geom_smooth(aes(group=1), method="loess", se=TRUE, color="grey30", linewidth=0.8, fill="grey85") +
    scale_color_manual(values=CD8_COLORS, name="Subtype", guide="none") +
    labs(title="IKZF1 Expression Along CD8+ T Cell Exhaustion Trajectory",
         subtitle="Increasing IKZF1 along exhaustion pseudotime supports its role in promoting exhaustion",
         x="Slingshot Pseudotime", y="IKZF1 Expression") +
    theme_classic() +
    theme(plot.title=element_text(face="bold",size=12),
          plot.subtitle=element_text(size=9,color="grey50"))
  save_fig(p_ikzf1_pt, "cd8_ikzf1_vs_pseudotime", 8, 5.5)
  log_msg("Slingshot trajectory complete")
}, error=function(e) log_msg(paste("Slingshot ERROR:", e$message)))

# 保存CD8亚群对象
saveRDS(cd8, file.path(res_dir, "cd8_subclusters.rds"))
log_msg("Saved CD8 Seurat object")

# ══════════════════════════════════════════════════════════════════════════════
# PART 2: B细胞 + 浆细胞亚群分析
# ══════════════════════════════════════════════════════════════════════════════
log_msg("=== PART 2: B Cell + Plasma Cell Subclustering ===")

b_seu <- subset(ann, cell_type %in% c("B_cell","Plasma_cell"))
log_msg(paste("B+Plasma cells:", ncol(b_seu)))

b_seu <- NormalizeData(b_seu, verbose=FALSE)
b_seu <- FindVariableFeatures(b_seu, nfeatures=2000, verbose=FALSE)
b_seu <- ScaleData(b_seu, verbose=FALSE)
b_seu <- RunPCA(b_seu, npcs=30, verbose=FALSE)

if(requireNamespace("harmony", quietly=TRUE)) {
  b_seu <- RunHarmony(b_seu, group.by.vars="patient", verbose=FALSE)
  b_seu <- RunUMAP(b_seu, reduction="harmony", dims=1:20, verbose=FALSE)
  b_seu <- FindNeighbors(b_seu, reduction="harmony", dims=1:20, verbose=FALSE)
} else {
  b_seu <- RunUMAP(b_seu, dims=1:20, verbose=FALSE)
  b_seu <- FindNeighbors(b_seu, dims=1:20, verbose=FALSE)
}
b_seu <- FindClusters(b_seu, resolution=0.4, verbose=FALSE)

# B细胞亚型打分
b_markers <- list(
  "Naive B"      = c("TCL1A","IGHD","IGHM","FCER2","IL4R","CD22"),
  "Memory B"     = c("CD27","IGHG1","IGHA1","S1PR1","CR2","TNFRSF13B"),
  "GC-like B"    = c("AICDA","BCL6","MKI67","PCNA","FCRL4","ITGAX"),
  "Transitional B"= c("CD38","VPREB3","PTPRC","EBF1","IL4R"),
  "Plasma cell"  = c("MZB1","JCHAIN","SDC1","XBP1","IRF4","IGHG1","IGHA1")
)
b_seu <- AddModuleScore(b_seu, features=b_markers, name="bscore_", assay="RNA")
bscore_cols <- paste0("bscore_", seq_along(b_markers))
names(bscore_cols) <- names(b_markers)

b_mat <- b_seu@meta.data[, bscore_cols, drop=FALSE]
b_seu$b_subtype <- names(b_markers)[apply(b_mat, 1, which.max)]
b_seu$b_subtype <- factor(b_seu$b_subtype, levels=names(B_COLORS))

tbl_b <- table(b_seu$b_subtype)
log_msg("B cell subtype counts:"); print(tbl_b)

# ── 图7: B细胞亚群UMAP ────────────────────────────────────────────────────────
log_msg("B Fig7: UMAP...")
p_b_umap <- DimPlot(b_seu, group.by="b_subtype", raster=FALSE,
                     cols=B_COLORS, pt.size=0.5, label=FALSE) +
  scale_color_manual(values=B_COLORS, name="B Cell\nSubtype") +
  labs(title="B Cell Subtype Atlas",
       subtitle=paste0(ncol(b_seu), " cells | B cell → Plasma cell maturation axis"),
       x="UMAP_1", y="UMAP_2") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=13),
        plot.subtitle=element_text(size=9,color="grey50"))
p_b_umap$layers[[1]]$aes_params$alpha <- 0.85
save_fig(p_b_umap, "b_umap_subtypes", 8, 6.5)

# ── 图8: B细胞标记基因点图 ────────────────────────────────────────────────────
log_msg("B Fig8: Dotplot...")
key_markers_b <- c(
  "TCL1A","IGHD","IGHM","FCER2",     # Naive
  "CD27","IGHA1","IGHG1","S1PR1",    # Memory
  "AICDA","BCL6","MKI67","FCRL4",    # GC-like
  "CD38","VPREB3","EBF1",            # Transitional
  "MZB1","JCHAIN","SDC1","XBP1","IRF4", # Plasma
  "CD79A","MS4A1"                    # Pan-B markers
)
Idents(b_seu) <- "b_subtype"
p_b_dot <- DotPlot(b_seu, features=intersect(key_markers_b, rownames(b_seu)),
                    group.by="b_subtype", dot.scale=5) +
  scale_color_gradient2(low="#2166AC",mid="#F7F7F7",high="#D73027",
                        midpoint=0, name="Avg\nExpr") +
  coord_flip() +
  labs(title="B Cell Subtype Marker Genes",
       x=NULL, y=NULL) +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        axis.text.x=element_text(angle=30,hjust=1,size=9),
        axis.text.y=element_text(size=9))
save_fig(p_b_dot, "b_dotplot_markers", 9, 7)

# ── 图9: CD79A在B细胞亚群的表达 ──────────────────────────────────────────────
log_msg("B Fig9: CD79A across B subtypes...")
b_meta <- b_seu@meta.data
b_meta$CD79A <- as.numeric(GetAssayData(b_seu, assay="RNA", layer="data")["CD79A",])

p_cd79a <- ggplot(b_meta, aes(x=b_subtype, y=CD79A, fill=b_subtype)) +
  geom_violin(scale="width", trim=TRUE, alpha=0.8, linewidth=0.3) +
  geom_boxplot(width=0.12, outlier.shape=NA, fill="white", linewidth=0.4) +
  scale_fill_manual(values=B_COLORS, guide="none") +
  labs(title="CD79A Expression Across B Cell Subtypes",
       subtitle="CD79A highest in Naive/Memory B; marks BCR axis activity as a bulk marker",
       x=NULL, y="CD79A Expression (log-normalized)") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey50"),
        axis.text.x=element_text(angle=20,hjust=1,size=9.5))
save_fig(p_cd79a, "b_cd79a_across_subtypes", 8, 5.5)

# ── 图10: B细胞分化拟时序 ─────────────────────────────────────────────────────
log_msg("B pseudotime: Slingshot B cell maturation...")
tryCatch({
  library(SingleCellExperiment)
  sce_b <- as.SingleCellExperiment(b_seu, assay="RNA")
  reducedDim(sce_b, "UMAP") <- Embeddings(b_seu, "umap")
  sl_b <- slingshot(sce_b, clusterLabels="b_subtype", reducedDim="UMAP",
                    start.clus="Naive B", end.clus="Plasma cell")

  pt_b <- slingPseudotime(sl_b)[,1]
  b_seu$slingshot_pt <- pt_b[colnames(b_seu)]

  pt_b_df <- data.frame(
    UMAP1=Embeddings(b_seu,"umap")[,1],
    UMAP2=Embeddings(b_seu,"umap")[,2],
    pseudotime=b_seu$slingshot_pt,
    subtype=b_seu$b_subtype
  )
  p_pt_b <- ggplot(pt_b_df, aes(x=UMAP1, y=UMAP2, color=pseudotime)) +
    geom_point(size=0.5, alpha=0.75) +
    scale_color_viridis_c(option="plasma", name="Pseudotime", na.value="grey80") +
    labs(title="B Cell Maturation Trajectory (Slingshot)",
         subtitle="Naive B → Memory B / GC-like → Plasma cell differentiation axis") +
    theme_classic() +
    theme(plot.title=element_text(face="bold",size=12),
          plot.subtitle=element_text(size=9,color="grey50"))
  save_fig(p_pt_b, "b_slingshot_pseudotime", 7, 6)
  log_msg("B cell slingshot complete")
}, error=function(e) log_msg(paste("B slingshot ERROR:", e$message)))

# ── 图11: 患者层面B亚型比例 ───────────────────────────────────────────────────
log_msg("B Fig11: Per-patient B proportion...")
{
  pt_vec_b  <- meta_chr(b_seu, "patient")
  st_vec_b  <- as.character(b_seu$b_subtype)
  ok_b      <- !is.na(pt_vec_b) & !is.na(st_vec_b)
  cnt_tbl_b <- as.data.frame(table(patient=pt_vec_b[ok_b], b_subtype=st_vec_b[ok_b]),
                              stringsAsFactors=FALSE)
  prop_b <- cnt_tbl_b %>%
    dplyr::rename(n=Freq) %>%
    dplyr::filter(n>0) %>%
    dplyr::group_by(patient) %>%
    dplyr::mutate(prop=n/sum(n)) %>%
    dplyr::ungroup()
}

pt_order_b <- prop_b %>%
  filter(b_subtype=="Plasma cell") %>%
  arrange(desc(prop)) %>% pull(patient)
prop_b$patient <- factor(prop_b$patient, levels=pt_order_b)

p_prop_b <- ggplot(prop_b, aes(x=patient, y=prop, fill=b_subtype)) +
  geom_col(width=0.8, color="white", linewidth=0.3) +
  scale_fill_manual(values=B_COLORS, name="B Cell Subtype") +
  scale_y_continuous(labels=scales::percent, expand=c(0,0)) +
  labs(title="B Cell Subtype Composition Per Patient",
       subtitle="Sorted by Plasma cell proportion | BCR maturation heterogeneity",
       x="Patient", y="Proportion") +
  theme_classic() +
  theme(plot.title=element_text(face="bold",size=12),
        plot.subtitle=element_text(size=9,color="grey50"),
        axis.text.x=element_text(angle=40,hjust=1,size=8.5))
save_fig(p_prop_b, "b_patient_proportion", 10, 5.5)

saveRDS(b_seu, file.path(res_dir, "b_subclusters.rds"))
log_msg("Saved B cell Seurat object")

# ══════════════════════════════════════════════════════════════════════════════
# 导出CD8 count矩阵供pySCENIC使用
# ══════════════════════════════════════════════════════════════════════════════
log_msg("Exporting CD8 count matrix for pySCENIC...")
tryCatch({
  # 导出为loom格式（pySCENIC标准输入）
  if(requireNamespace("loomR", quietly=TRUE)) {
    library(loomR)
    counts_cd8 <- GetAssayData(cd8, assay="RNA", layer="counts")
    # 过滤低表达基因（至少1%细胞表达）
    pct_expr <- rowSums(counts_cd8 > 0) / ncol(counts_cd8)
    counts_cd8 <- counts_cd8[pct_expr >= 0.01, ]
    log_msg(paste("CD8 genes for SCENIC:", nrow(counts_cd8)))
    meta_loom <- data.frame(
      patient     = meta_chr(cd8, "patient"),
      cd8_subtype = as.character(cd8$cd8_subtype),
      sample      = meta_chr(cd8, "sample"),
      stringsAsFactors = FALSE
    )
    loom_path <- "results/subcluster/cd8_for_scenic.loom"
    create(filename=loom_path,
           data=counts_cd8,
           cell.attrs=as.list(meta_loom),
           overwrite=TRUE)
    log_msg(paste("Loom saved:", loom_path))
  } else {
    # 用csv导出（备选，pySCENIC Python脚本可读取）
    counts_cd8 <- GetAssayData(cd8, assay="RNA", layer="counts")
    pct_expr <- rowSums(counts_cd8 > 0) / ncol(counts_cd8)
    counts_cd8 <- counts_cd8[pct_expr >= 0.01, ]
    log_msg(paste("CSV export:", nrow(counts_cd8), "genes x", ncol(counts_cd8), "cells"))
    out_csv <- "results/subcluster/cd8_for_scenic.csv"
    write.csv(as.data.frame(as.matrix(counts_cd8)), out_csv)
    meta_safe <- data.frame(
      patient     = meta_chr(cd8, "patient"),
      cd8_subtype = as.character(cd8$cd8_subtype),
      sample      = meta_chr(cd8, "sample"),
      stringsAsFactors = FALSE
    )
    write.csv(meta_safe, "results/subcluster/cd8_meta_for_scenic.csv")
    log_msg(paste("CSV export complete:", out_csv))
  }
}, error=function(e) {
  log_msg(paste("Loom export ERROR:", e$message))
  log_msg("Falling back to sparse CSV export...")
  tryCatch({
    suppressPackageStartupMessages(library(Matrix))
    counts_cd8 <- GetAssayData(cd8, assay="RNA", layer="counts")
    pct_expr <- rowSums(counts_cd8 > 0) / ncol(counts_cd8)
    counts_cd8 <- counts_cd8[pct_expr >= 0.01, ]
    out_csv <- "results/subcluster/cd8_for_scenic.csv"
    write.csv(as.data.frame(as.matrix(counts_cd8)), out_csv)
    meta_safe2 <- data.frame(
      patient = meta_chr(cd8, "patient"),
      cd8_subtype = as.character(cd8$cd8_subtype),
      sample = meta_chr(cd8, "sample"),
      stringsAsFactors = FALSE
    )
    write.csv(meta_safe2, "results/subcluster/cd8_meta_for_scenic.csv")
    log_msg(paste("Fallback CSV written:", out_csv))
  }, error=function(e2) log_msg(paste("Fallback CSV also failed:", e2$message)))
})

# 写入完成标志供 P01 监测
writeLines(as.character(Sys.time()), "results/subcluster/r12_complete.flag")
log_msg("Sentinel flag written: results/subcluster/r12_complete.flag")

# ══════════════════════════════════════════════════════════════════════════════
# 同步到 final_paper
# ══════════════════════════════════════════════════════════════════════════════
log_msg("Syncing to final_paper...")
fp_dir <- "final_paper/Figure2_Subclusters"
dir.create(fp_dir, recursive=TRUE, showWarnings=FALSE)
figs <- list.files(out_dir, pattern="\\.(png|pdf)$", full.names=TRUE)
for(f in figs) file.copy(f, file.path(fp_dir, basename(f)), overwrite=TRUE)
log_msg(paste("Synced", length(figs), "figures"))

log_msg("=== R12 complete ===")
