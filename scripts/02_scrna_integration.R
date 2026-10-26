#!/usr/bin/env Rscript

# ==========================================================
# Script   : 02_scrna_integration.R
# Purpose  : Merge single-cell objects and perform donor-aware integration
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: single-cell data processing; Fig. 1
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/01_single_cell/01b_merge_harmony.R
# ==========================================================
# =============================================================================
# 脚本: 01b_merge_harmony.R
# 功能: 合并所有样本 Seurat 对象，Harmony 批次校正整合
# 输入: results/01_scRNAseq/QC/seurat_per_sample_postQC.rds
# 输出: results/01_scRNAseq/clustering/seurat_harmony.rds
# 预计内存: ~80 GB
# 预计时间: 1-2小时
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(future)
})

# ── 路径设置 ──────────────────────────────────────────────────────────────────
base_path    <- "."
qc_out       <- file.path(base_path, "results/01_scRNAseq/QC")
clust_out    <- file.path(base_path, "results/01_scRNAseq/clustering")
fig_out      <- file.path(base_path, "figures/module01")
log_file     <- file.path(base_path, "logs/01b_harmony.log")

dir.create(clust_out, recursive=TRUE, showWarnings=FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file=log_file, append=TRUE)
}

log_msg("=== Module 1b: Merge + Harmony 开始 ===")

# ── 检查上游输出 ──────────────────────────────────────────────────────────────
input_file <- file.path(qc_out, "seurat_per_sample_postQC.rds")
if (!file.exists(input_file)) stop("缺少输入文件: ", input_file, "\n先运行 01a_qc_per_sample.R")
log_msg(paste("加载输入文件:", input_file))

seurat_list <- readRDS(input_file)
log_msg(paste("样本数:", length(seurat_list)))
log_msg(paste("样本列表:", paste(names(seurat_list), collapse=", ")))

# ── Step 1: 合并所有样本 ──────────────────────────────────────────────────────
log_msg("Step 1: 合并所有 Seurat 对象 ...")

# Seurat v5 merge
merged <- merge(
  x    = seurat_list[[1]],
  y    = seurat_list[-1],
  add.cell.ids = names(seurat_list),
  merge.data   = TRUE
)

# Seurat v5: 合并后需要 JoinLayers
merged <- JoinLayers(merged)

log_msg(paste("合并后总细胞数:", ncol(merged)))
log_msg(paste("总基因数:", nrow(merged)))

# ── Step 2: 标准化 + 高变基因 ────────────────────────────────────────────────
log_msg("Step 2: 标准化 + 高变基因筛选 ...")

merged <- NormalizeData(merged, normalization.method="LogNormalize",
                        scale.factor=10000, verbose=FALSE)
merged <- FindVariableFeatures(merged, selection.method="vst",
                               nfeatures=3000, verbose=FALSE)

log_msg(paste("高变基因数:", length(VariableFeatures(merged))))

# ── Step 3: 标准化 + PCA ─────────────────────────────────────────────────────
log_msg("Step 3: ScaleData + PCA ...")

# 只对高变基因 scale（节省内存）
merged <- ScaleData(merged, features=VariableFeatures(merged),
                    vars.to.regress=c("pct_mt", "nCount_RNA"), verbose=FALSE)
merged <- RunPCA(merged, npcs=50, verbose=FALSE)

# PCA 解释方差
pca_var <- (merged@reductions$pca@stdev)^2
pca_var_pct <- pca_var / sum(pca_var) * 100
cum_var <- cumsum(pca_var_pct)
n_pcs <- which(cum_var >= 90)[1]
log_msg(paste("PC1-50 累积方差:", round(cum_var[50], 1), "%"))
log_msg(paste("达到90%累积方差所需PC数:", n_pcs))

# Elbow plot
p_elbow <- ElbowPlot(merged, ndims=50) +
  geom_vline(xintercept=30, linetype="dashed", color="red") +
  labs(title="PCA Elbow Plot - Variance Explained by PC",
       x="Principal Component", y="Standard Deviation") +
  theme_bw()
ggsave(file.path(fig_out, "pca_elbow.pdf"), p_elbow, width=8, height=5)
ggsave(file.path(fig_out, "pca_elbow.png"), p_elbow, width=8, height=5, dpi=150)
log_msg("  图: pca_elbow — PCA各主成分解释方差，红虚线=PC30")

# ── Step 4: UMAP（未校正，用于对比）─────────────────────────────────────────
log_msg("Step 4: 计算未校正 UMAP（对比用）...")
merged <- RunUMAP(merged, dims=1:30, reduction="pca",
                  reduction.name="umap.uncorrected", verbose=FALSE)

p_pre_harmony <- DimPlot(merged, reduction="umap.uncorrected",
                          group.by="sample", pt.size=0.1) +
  labs(title="UMAP before Harmony (colored by sample)") +
  theme_bw() + theme(legend.position="right")
ggsave(file.path(fig_out, "umap_pre_harmony_sample.pdf"), p_pre_harmony, width=12, height=8)
ggsave(file.path(fig_out, "umap_pre_harmony_sample.png"), p_pre_harmony, width=12, height=8, dpi=150)
log_msg("  图: umap_pre_harmony_sample — Harmony校正前UMAP，按样本着色，显示批次效应")

p_pre_frac <- DimPlot(merged, reduction="umap.uncorrected",
                       group.by="fraction", pt.size=0.1) +
  scale_color_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8")) +
  labs(title="UMAP before Harmony (colored by fraction)") +
  theme_bw()
ggsave(file.path(fig_out, "umap_pre_harmony_fraction.pdf"), p_pre_frac, width=10, height=8)
ggsave(file.path(fig_out, "umap_pre_harmony_fraction.png"), p_pre_frac, width=10, height=8, dpi=150)
log_msg("  图: umap_pre_harmony_fraction — Harmony校正前UMAP，CD45p（红）vs CD45n（蓝）分布")

# ── Step 5: Harmony 批次校正 ──────────────────────────────────────────────────
log_msg("Step 5: Harmony 批次校正（batch=sample）...")

merged <- RunHarmony(
  merged,
  group.by.vars   = "sample",
  reduction       = "pca",
  dims.use        = 1:30,
  reduction.save  = "harmony",
  theta           = 2,
  max.iter.harmony = 50,
  verbose         = TRUE
)

log_msg("Harmony 完成")

# ── Step 6: UMAP（Harmony校正后）─────────────────────────────────────────────
log_msg("Step 6: 计算 Harmony 校正后 UMAP ...")

merged <- RunUMAP(merged, dims=1:30, reduction="harmony",
                  reduction.name="umap", verbose=FALSE)

# UMAP by sample
p_post_sample <- DimPlot(merged, reduction="umap", group.by="sample",
                           pt.size=0.1) +
  labs(title="UMAP after Harmony (colored by sample)") +
  theme_bw()
ggsave(file.path(fig_out, "umap_post_harmony_sample.pdf"), p_post_sample, width=12, height=8)
ggsave(file.path(fig_out, "umap_post_harmony_sample.png"), p_post_sample, width=12, height=8, dpi=150)
log_msg("  图: umap_post_harmony_sample — Harmony校正后UMAP，各样本充分混合说明批次校正成功")

# UMAP by fraction
p_post_frac <- DimPlot(merged, reduction="umap", group.by="fraction",
                        pt.size=0.1) +
  scale_color_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8")) +
  labs(title="UMAP after Harmony (colored by fraction)") +
  theme_bw()
ggsave(file.path(fig_out, "umap_post_harmony_fraction.pdf"), p_post_frac, width=10, height=8)
ggsave(file.path(fig_out, "umap_post_harmony_fraction.png"), p_post_frac, width=10, height=8, dpi=150)
log_msg("  图: umap_post_harmony_fraction — Harmony校正后CD45p/CD45n分布，两者既有分离也有共享区域")

# UMAP by patient
p_post_patient <- DimPlot(merged, reduction="umap", group.by="patient",
                            pt.size=0.1) +
  labs(title="UMAP after Harmony (colored by patient)") +
  theme_bw()
ggsave(file.path(fig_out, "umap_post_harmony_patient.pdf"), p_post_patient, width=12, height=8)
ggsave(file.path(fig_out, "umap_post_harmony_patient.png"), p_post_patient, width=12, height=8, dpi=150)
log_msg("  图: umap_post_harmony_patient — Harmony校正后UMAP，按患者着色")

# ── Step 7: 保存 ──────────────────────────────────────────────────────────────
log_msg("Step 7: 保存整合后 Seurat 对象 ...")
saveRDS(merged, file.path(clust_out, "seurat_harmony.rds"))

# ── 自检 ─────────────────────────────────────────────────────────────────────
log_msg("=== 自检 Module 1b ===")
log_msg(paste("总细胞数:", ncol(merged)))
log_msg(paste("Harmony reduction dims:", ncol(merged@reductions$harmony@cell.embeddings)))
if (!"harmony" %in% names(merged@reductions)) {
  log_msg("⚠ 错误: Harmony reduction 缺失")
} else {
  log_msg("✓ Harmony reduction 存在")
}

log_msg("=== Module 1b 完成 ===")
log_msg(paste("输出: results/01_scRNAseq/clustering/seurat_harmony.rds"))
