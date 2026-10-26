#!/usr/bin/env Rscript

# ==========================================================
# Script   : 07_tcga_wgcna_yellow.R
# Purpose  : Reproduce TCGA WGCNA and the historical Yellow module
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: WGCNA; Fig. 2A-D
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/02_wgcna/02c_wgcna.R
# ==========================================================
# =============================================================================
# 脚本: 02c_wgcna.R
# 功能: WGCNA 网络构建 + 模块-性状关联
# 输入: results/02_WGCNA/vst_top25mad.rds
#        results/02_WGCNA/trait_matrix/trait_matrix.rds
# 输出: results/02_WGCNA/network/wgcna_results.rds
#        figures/module02/wgcna_*.pdf/png
# 预计内存: ~60 GB
# 预计时间: 2-4小时
# =============================================================================

suppressPackageStartupMessages({
  library(WGCNA)
  library(ggplot2)
  library(dplyr)
  library(reshape2)
})

# WGCNA 全局设置
options(stringsAsFactors=FALSE)
enableWGCNAThreads(nThreads=32)  # 32核并行（保留余量）

# ── 路径设置 ──────────────────────────────────────────────────────────────────
base_path   <- "."
wgcna_out   <- file.path(base_path, "results/02_WGCNA")
network_out <- file.path(wgcna_out, "network")
trait_dir   <- file.path(wgcna_out, "trait_matrix")
fig_out     <- file.path(base_path, "figures/module02")
log_file    <- file.path(base_path, "logs/02c_wgcna.log")

dir.create(network_out, recursive=TRUE, showWarnings=FALSE)

log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file=log_file, append=TRUE)
}

log_msg("=== Module 2c: WGCNA 网络构建 开始 ===")
log_msg(paste("WGCNA version:", packageVersion("WGCNA")))

# ── 检查上游输出 ──────────────────────────────────────────────────────────────
vst_file   <- file.path(wgcna_out, "vst_top25mad.rds")
trait_file <- file.path(trait_dir, "trait_matrix.rds")
if (!file.exists(vst_file))   stop("缺少输入: ", vst_file)
if (!file.exists(trait_file)) stop("缺少输入: ", trait_file)

log_msg("加载数据 ...")
vst_mat   <- readRDS(vst_file)
trait_df  <- readRDS(trait_file)

# 转置（WGCNA 需要 samples x genes）
datExpr <- t(vst_mat)
log_msg(paste("表达矩阵:", nrow(datExpr), "样本 x", ncol(datExpr), "基因"))

# 对齐样本
common_samples <- intersect(rownames(datExpr), trait_df$sample_id)
datExpr  <- datExpr[common_samples, ]
trait_df <- trait_df[trait_df$sample_id %in% common_samples, ]
rownames(trait_df) <- trait_df$sample_id
trait_df <- trait_df[common_samples, ]
log_msg(paste("对齐后样本数:", length(common_samples)))

# ── Step 1: 检查基因/样本质量 ────────────────────────────────────────────────
log_msg("Step 1: 检查基因/样本质量 ...")

gsg <- goodSamplesGenes(datExpr, verbose=0)
if (!gsg$allOK) {
  n_bad_genes   <- sum(!gsg$goodGenes)
  n_bad_samples <- sum(!gsg$goodSamples)
  log_msg(paste("移除:", n_bad_genes, "个坏基因,", n_bad_samples, "个坏样本"))
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
}
log_msg(paste("质量检查后:", nrow(datExpr), "样本 x", ncol(datExpr), "基因"))

# 样本树（检测离群样本）
sampleTree <- hclust(dist(datExpr), method="average")
pdf(file.path(fig_out, "wgcna_sample_tree.pdf"), width=16, height=8)
plot(sampleTree, main="Sample Clustering (outlier detection)",
     sub="", xlab="", ylab="Height", cex.axis=1, cex.main=1, cex=0.6)
abline(h=quantile(sampleTree$height, 0.95), col="red", lty=2)
dev.off()
png(file.path(fig_out, "wgcna_sample_tree.png"), width=1600, height=800, res=100)
plot(sampleTree, main="Sample Clustering (outlier detection)",
     sub="", xlab="", ylab="Height", cex.axis=1, cex.main=1, cex=0.6)
abline(h=quantile(sampleTree$height, 0.95), col="red", lty=2)
dev.off()
log_msg("  图: wgcna_sample_tree — 样本层次聚类树（红虚线=95百分位高度，用于识别离群样本）")

# 切除离群样本（高度>95百分位）
height_cutoff <- quantile(sampleTree$height, 0.99)  # 保守移除（99%）
clust <- cutreeStatic(sampleTree, cutHeight=height_cutoff, minSize=10)
keep_samples <- (clust == 1)
n_removed <- sum(!keep_samples)
if (n_removed > 0) {
  log_msg(paste("移除离群样本:", n_removed, "个"))
  datExpr  <- datExpr[keep_samples, ]
  trait_df <- trait_df[keep_samples, ]
}
log_msg(paste("最终分析样本数:", nrow(datExpr)))

# ── Step 2: 选择软阈值 ───────────────────────────────────────────────────────
log_msg("Step 2: 选择软阈值（soft-thresholding power）...")

powers <- c(1:10, seq(12, 20, by=2))
sft <- pickSoftThreshold(datExpr, powerVector=powers, verbose=0,
                          networkType="unsigned", corFnc="bicor")

# 找到 R² > 0.85 的最小 power
sft_df      <- sft$fitIndices
best_power  <- sft_df$Power[sft_df$SFT.R.sq >= 0.85][1]
if (is.na(best_power)) {
  best_power <- sft_df$Power[which.max(sft_df$SFT.R.sq)]
  log_msg(paste("⚠ R²<0.85，使用最高R²对应power:", best_power,
                "(R²=", round(max(sft_df$SFT.R.sq),3), ")"))
} else {
  log_msg(paste("最优 soft power:", best_power,
                "(R²=", round(sft_df$SFT.R.sq[sft_df$Power==best_power],3), ")"))
}

# Soft threshold 图
p1 <- ggplot(sft_df, aes(x=Power, y=SFT.R.sq)) +
  geom_point(color="#377EB8", size=3) +
  geom_line(color="#377EB8") +
  geom_label(aes(label=Power), size=3) +
  geom_hline(yintercept=0.85, linetype="dashed", color="red") +
  geom_vline(xintercept=best_power, linetype="dashed", color="darkgreen") +
  theme_bw() +
  labs(title="Scale-free Topology Fit",
       x="Soft Threshold (power)", y=expression(R^2))

p2 <- ggplot(sft_df, aes(x=Power, y=mean.k.)) +
  geom_point(color="#E41A1C", size=3) +
  geom_line(color="#E41A1C") +
  geom_vline(xintercept=best_power, linetype="dashed", color="darkgreen") +
  theme_bw() +
  labs(title="Mean Connectivity",
       x="Soft Threshold (power)", y="Mean Connectivity")

library(patchwork)
p_sft <- p1 + p2
ggsave(file.path(fig_out, "wgcna_soft_threshold.pdf"), p_sft, width=14, height=6)
ggsave(file.path(fig_out, "wgcna_soft_threshold.png"), p_sft, width=14, height=6, dpi=150)
log_msg("  图: wgcna_soft_threshold — 软阈值选择图，左：无尺度拓扑R²，右：平均连接度")

# ── Step 3: WGCNA 网络构建 ───────────────────────────────────────────────────
log_msg(paste("Step 3: blockwiseModules（power=", best_power, "，32核并行）..."))

net <- blockwiseModules(
  datExpr,
  power             = best_power,
  networkType       = "unsigned",
  corType           = "bicor",
  TOMType           = "unsigned",
  minModuleSize     = 30,
  reassignThreshold = 0,
  mergeCutHeight    = 0.25,
  numericLabels     = FALSE,
  pamRespectsDendro = FALSE,
  maxBlockSize      = 20000,
  nThreads          = 32,
  saveTOMs          = TRUE,
  saveTOMFileBase   = file.path(network_out, "TOM"),
  verbose           = 3
)

module_colors <- net$colors
n_modules     <- length(unique(module_colors)) - 1  # 不计 grey
log_msg(paste("识别模块数（不含grey）:", n_modules))
log_msg(paste("模块颜色:", paste(unique(module_colors)[unique(module_colors)!="grey"],
                                collapse=", ")))

# 模块大小统计
module_sizes <- table(module_colors)
log_msg(paste("模块大小分布（前10大）:",
              paste(sort(module_sizes, decreasing=TRUE)[1:min(10,length(module_sizes))],
                    collapse=", ")))

# ── Step 4: 聚类树 + 模块颜色图 ──────────────────────────────────────────────
log_msg("Step 4: 绘制模块聚类树 ...")

mergedMEs <- net$MEs
sampleTree2 <- hclust(dist(datExpr), method="average")

pdf(file.path(fig_out, "wgcna_dendrogram.pdf"), width=20, height=8)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels=FALSE, hang=0.03,
  addGuide=TRUE, guideHang=0.05,
  main="Gene Clustering Dendrogram and Module Colors"
)
dev.off()
png(file.path(fig_out, "wgcna_dendrogram.png"), width=2000, height=800, res=100)
plotDendroAndColors(
  net$dendrograms[[1]],
  module_colors[net$blockGenes[[1]]],
  "Module colors",
  dendroLabels=FALSE, hang=0.03,
  addGuide=TRUE, guideHang=0.05,
  main="Gene Clustering Dendrogram and Module Colors"
)
dev.off()
log_msg("  图: wgcna_dendrogram — 基因聚类树与模块颜色条，展示各基因所属模块")

# ── Step 5: 模块特征向量 (ME) 相关性 ─────────────────────────────────────────
log_msg("Step 5: 计算模块特征向量相关性 ...")

ME_order <- orderMEs(mergedMEs)
ME_cor   <- cor(ME_order, use="p")
ME_p     <- corPvalueStudent(ME_cor, nrow(datExpr))

pdf(file.path(fig_out, "wgcna_me_correlation.pdf"), width=12, height=10)
plotEigengeneNetworks(ME_order, "",
  marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),
  plotDendrograms=TRUE, xLabelsAngle=90)
dev.off()
png(file.path(fig_out, "wgcna_me_correlation.png"), width=1200, height=1000, res=100)
plotEigengeneNetworks(ME_order, "",
  marDendro=c(0,4,1,2), marHeatmap=c(3,4,1,2),
  plotDendrograms=TRUE, xLabelsAngle=90)
dev.off()
log_msg("  图: wgcna_me_correlation — 模块特征向量相关性网络（ME间的相似性）")

# ── Step 6: 模块-性状相关性 ──────────────────────────────────────────────────
log_msg("Step 6: 计算模块-性状相关性 ...")

# 提取数值型 trait 列
trait_num_cols <- colnames(trait_df)[sapply(trait_df, is.numeric)]
datTraits      <- trait_df[common_samples[rownames(datExpr) %in% common_samples],
                            trait_num_cols]
datTraits      <- datTraits[rownames(datExpr), , drop=FALSE]

# 计算相关性
moduleTraitCor <- cor(mergedMEs, datTraits, use="p")
moduleTraitP   <- corPvalueStudent(moduleTraitCor, nrow(datExpr))

# 显著性标注
sig_text <- paste(signif(moduleTraitCor, 2), "\n(",
                  signif(moduleTraitP, 1), ")", sep="")
dim(sig_text) <- dim(moduleTraitCor)

# 热图
pdf(file.path(fig_out, "wgcna_module_trait_heatmap.pdf"), width=16, height=12)
par(mar=c(6,8.5,3,3))
labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(datTraits),
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = sig_text,
  setStdMargins = FALSE,
  cex.text    = 0.5,
  zlim        = c(-1, 1),
  main        = "Module-Trait Relationships"
)
dev.off()
png(file.path(fig_out, "wgcna_module_trait_heatmap.png"), width=1600, height=1200, res=100)
par(mar=c(6,8.5,3,3))
labeledHeatmap(
  Matrix      = moduleTraitCor,
  xLabels     = colnames(datTraits),
  yLabels     = rownames(moduleTraitCor),
  ySymbols    = rownames(moduleTraitCor),
  colorLabels = FALSE,
  colors      = blueWhiteRed(50),
  textMatrix  = sig_text,
  setStdMargins = FALSE,
  cex.text    = 0.5,
  zlim        = c(-1, 1),
  main        = "Module-Trait Relationships"
)
dev.off()
log_msg("  图: wgcna_module_trait_heatmap — 模块-性状相关性热图（颜色=相关系数，括号内=p值）")

# 识别关键模块（与 Immune_Score 相关 |r|>0.4, p<0.05）
key_trait    <- "Immune_Score"
if (key_trait %in% colnames(moduleTraitCor)) {
  immune_cor   <- moduleTraitCor[, key_trait]
  immune_p     <- moduleTraitP[, key_trait]
  key_modules  <- rownames(moduleTraitCor)[abs(immune_cor) > 0.4 & immune_p < 0.05]
  key_modules  <- gsub("^ME", "", key_modules)

  log_msg(paste("与Immune_Score显著相关的模块:", paste(key_modules, collapse=", ")))
  log_msg(paste("关联r值:", paste(round(immune_cor[abs(immune_cor)>0.4 & immune_p<0.05], 3),
                                   collapse=", ")))
  writeLines(key_modules, file.path(network_out, "key_modules.txt"))
} else {
  key_modules <- gsub("^ME", "", rownames(moduleTraitCor)[1:3])  # 默认取前3个
  log_msg("⚠ Immune_Score列不存在，使用前3个模块作为关键模块")
}

# ── Step 7: 保存结果 ──────────────────────────────────────────────────────────
log_msg("Step 7: 保存 WGCNA 结果 ...")

wgcna_results <- list(
  net              = net,
  module_colors    = module_colors,
  MEs              = mergedMEs,
  moduleTraitCor   = moduleTraitCor,
  moduleTraitP     = moduleTraitP,
  key_modules      = key_modules,
  best_power       = best_power,
  n_samples        = nrow(datExpr),
  n_genes          = ncol(datExpr),
  n_modules        = n_modules,
  sample_order     = rownames(datExpr),
  gene_names       = colnames(datExpr)
)

saveRDS(wgcna_results, file.path(network_out, "wgcna_results.rds"))

# 基因-模块映射表
gene_module_df <- data.frame(
  gene   = colnames(datExpr),
  module = module_colors,
  stringsAsFactors = FALSE
)
write.csv(gene_module_df, file.path(network_out, "gene_module_mapping.csv"), row.names=FALSE)

# ── 自检 ─────────────────────────────────────────────────────────────────────
log_msg("=== 自检 Module 2c ===")
log_msg(paste("软阈值:", best_power, "(预期: 4-20)"))
if (best_power > 20) log_msg("⚠ 软阈值>20，数据质量可能存在问题")

log_msg(paste("模块数:", n_modules, "(预期: 10-30)"))
if (n_modules < 10) log_msg("⚠ 模块数<10，minModuleSize可能过大")
if (n_modules > 30) log_msg("⚠ 模块数>30，可考虑提高mergeCutHeight")

log_msg(paste("关键模块数:", length(key_modules), "(预期: ≥1)"))
if (length(key_modules) == 0) log_msg("⚠ 未找到与Immune_Score显著相关的模块")

log_msg("=== Module 2c 完成 ===")
log_msg("输出: results/02_WGCNA/network/wgcna_results.rds")
