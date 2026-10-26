#!/usr/bin/env Rscript

# ==========================================================
# Script   : 01_scrna_qc.R
# Purpose  : Single-cell quality control
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: single-cell data processing; Fig. 1; Supplementary Fig. 1
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/01_single_cell/01a_qc_per_sample.R
# ==========================================================
# =============================================================================
# 脚本: 01a_qc_per_sample.R
# 功能: GSE164690数据解压、组织、每样本QC过滤 + DoubletFinder去双细胞
# 输入: GSE164690_RAW.tar
# 输出: results/01_scRNAseq/QC/seurat_per_sample_postQC.rds
#        results/01_scRNAseq/doublet/doublet_stats.csv
#        figures/module01/qc_violin_*.pdf/png
# 预计内存: ~60 GB（所有样本同时加载）
# 预计时间: 2-4小时
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(dplyr)
  library(future)
})

# ── 路径设置（相对路径，从项目根目录运行）────────────────────────────────────
base_path      <- "."
tar_file       <- file.path(base_path, "GSE164690_RAW.tar")
raw_dir        <- file.path(base_path, "data/scRNAseq/GSE164690_raw")
sample_dir     <- file.path(base_path, "data/scRNAseq/GSE164690")
qc_out         <- file.path(base_path, "results/01_scRNAseq/QC")
doublet_out    <- file.path(base_path, "results/01_scRNAseq/doublet")
fig_out        <- file.path(base_path, "figures/module01")
log_file       <- file.path(base_path, "logs/01a_qc.log")

dir.create(qc_out,      recursive=TRUE, showWarnings=FALSE)
dir.create(doublet_out, recursive=TRUE, showWarnings=FALSE)
dir.create(fig_out,     recursive=TRUE, showWarnings=FALSE)
dir.create(raw_dir,     recursive=TRUE, showWarnings=FALSE)
dir.create(sample_dir,  recursive=TRUE, showWarnings=FALSE)

# ── 日志函数 ─────────────────────────────────────────────────────────────────
log_msg <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  line <- paste0("[", ts, "] ", msg)
  cat(line, "\n")
  cat(line, "\n", file=log_file, append=TRUE)
}

log_msg("=== Module 1a: QC + DoubletFinder 开始 ===")
log_msg(paste("R version:", R.version$major, R.version$minor))
log_msg(paste("Seurat version:", packageVersion("Seurat")))

# ── Step 1: 解压 tar 文件 ─────────────────────────────────────────────────────
log_msg("Step 1: 解压 GSE164690_RAW.tar ...")

# 检查是否已解压
existing_files <- list.files(raw_dir, pattern="\\.gz$")
if (length(existing_files) < 100) {
  log_msg(paste("解压到:", raw_dir))
  system(paste("tar -xf", tar_file, "-C", raw_dir))
  log_msg(paste("解压完成，文件数:", length(list.files(raw_dir))))
} else {
  log_msg(paste("已存在", length(existing_files), "个文件，跳过解压"))
}

# ── Step 2: 组织成 per-sample 目录结构 ──────────────────────────────────────
log_msg("Step 2: 组织 per-sample 目录结构 ...")

all_files <- normalizePath(list.files(raw_dir, full.names=TRUE), mustWork=FALSE)

# 解析文件名获取样本信息
# 格式: GSM*_HN##_CD45p/CD45n/PBL_barcodes/features/matrix.tsv/mtx.gz
file_info <- data.frame(
  full_path = all_files,
  filename  = basename(all_files),
  stringsAsFactors = FALSE
)

# 提取样本名 (e.g. HN01_CD45p) 和文件类型
file_info$sample <- gsub("^GSM\\d+_(HN\\d+_(?:CD45p|CD45n|PBL))_.*$", "\\1",
                          file_info$filename, perl=TRUE)
file_info$filetype <- gsub("^GSM\\d+_HN\\d+_(?:CD45p|CD45n|PBL)_(.*\\.(?:tsv|mtx)\\.gz)$", "\\1",
                            file_info$filename, perl=TRUE)

# 标准化文件名映射
# barcodes.tsv.gz → barcodes.tsv.gz
# features.tsv.gz → features.tsv.gz
# matrix.mtx.gz   → matrix.mtx.gz
get_std_name <- function(fn) {
  if (grepl("barcodes", fn)) return("barcodes.tsv.gz")
  if (grepl("features", fn)) return("features.tsv.gz")
  if (grepl("matrix", fn))   return("matrix.mtx.gz")
  return(fn)
}

# 创建符号链接或复制文件到 per-sample 目录
samples <- unique(file_info$sample)
log_msg(paste("发现样本数:", length(samples)))

for (smp in samples) {
  smp_dir <- file.path(sample_dir, smp)
  dir.create(smp_dir, recursive=TRUE, showWarnings=FALSE)

  smp_files <- file_info[file_info$sample == smp, ]

  for (i in seq_len(nrow(smp_files))) {
    src <- smp_files$full_path[i]
    std_name <- get_std_name(smp_files$filename[i])
    dst <- file.path(smp_dir, std_name)
    if (!file.exists(dst)) {
      file.symlink(src, dst)  # 使用符号链接节省磁盘空间
    }
  }
}
log_msg("Per-sample 目录组织完成")

# ── Step 3: 确定分析样本（仅CD45p + CD45n，不含PBL）────────────────────────
# CD45p: 18例 (HN01-HN18) — 免疫细胞
# CD45n: 15例 (HN05-HN18) — 肿瘤/CAF/内皮细胞（HN02/03/04缺失）
analysis_samples <- samples[!grepl("_PBL$", samples)]
cd45p_samples    <- sort(analysis_samples[grepl("_CD45p$", analysis_samples)])
cd45n_samples    <- sort(analysis_samples[grepl("_CD45n$", analysis_samples)])
log_msg(paste("CD45p样本:", length(cd45p_samples)))
log_msg(paste("CD45n样本:", length(cd45n_samples)))

# ── Step 4: 加载数据 + QC 过滤 + DoubletFinder ──────────────────────────────
log_msg("Step 4: 逐样本加载、QC、DoubletFinder ...")

# QC 参数
MIN_FEATURES  <- 200
MAX_FEATURES  <- 6000
MIN_COUNTS    <- 500
MAX_COUNTS    <- 50000
MAX_MT_PCT    <- 25
MIN_RIBO_PCT  <- 0   # 不过滤核糖体

# DoubletFinder 参数
DF_pN <- 0.25

# 存储结果
seurat_list    <- list()
doublet_stats  <- data.frame()
qc_stats_all   <- data.frame()

process_sample <- function(smp) {
  log_msg(paste("  处理样本:", smp))
  smp_dir <- file.path(sample_dir, smp)

  # 检查文件存在
  req_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")
  missing <- req_files[!file.exists(file.path(smp_dir, req_files))]
  if (length(missing) > 0) {
    log_msg(paste("  ⚠ 缺少文件:", paste(missing, collapse=", "), "，跳过", smp))
    return(NULL)
  }

  # 读取数据 (gene.column=2 使用基因Symbol)
  counts <- Read10X(smp_dir, gene.column=2)

  # 提取患者ID和分组信息
  patient <- gsub("_(CD45p|CD45n)$", "", smp)
  fraction <- gsub("^HN\\d+_", "", smp)

  # 创建 Seurat 对象
  seu <- CreateSeuratObject(
    counts   = counts,
    project  = smp,
    min.cells = 3,
    min.features = 100
  )

  # 添加元数据
  seu$sample   <- smp
  seu$patient  <- patient
  seu$fraction <- fraction

  # 计算线粒体基因百分比
  seu[["pct_mt"]] <- PercentageFeatureSet(seu, pattern="^MT-")

  # 计算核糖体基因百分比
  seu[["pct_ribo"]] <- PercentageFeatureSet(seu, pattern="^RP[SL]")

  # QC stats 过滤前
  n_pre_qc <- ncol(seu)

  # QC 过滤
  seu <- subset(seu,
    subset = nFeature_RNA >= MIN_FEATURES &
             nFeature_RNA <= MAX_FEATURES &
             nCount_RNA   >= MIN_COUNTS   &
             nCount_RNA   <= MAX_COUNTS   &
             pct_mt       <= MAX_MT_PCT
  )

  n_post_qc <- ncol(seu)

  if (n_post_qc < 100) {
    log_msg(paste("  ⚠ 样本", smp, "QC后仅剩", n_post_qc, "细胞，跳过"))
    return(NULL)
  }

  log_msg(paste("  ", smp, ": QC前", n_pre_qc, "→ QC后", n_post_qc, "细胞"))

  # ── DoubletFinder ────────────────────────────────────────────────────────
  # 预处理（DoubletFinder 需要 PCA）
  seu <- NormalizeData(seu, verbose=FALSE)
  seu <- FindVariableFeatures(seu, nfeatures=2000, verbose=FALSE)
  seu <- ScaleData(seu, verbose=FALSE)
  seu <- RunPCA(seu, npcs=20, verbose=FALSE)

  # 估计 pK (使用 paramSweep，计算代价较高但更准确)
  tryCatch({
    # 对于大样本 (>5000 cells) 采样加速
    if (ncol(seu) > 5000) {
      sweep_seu <- seu[, sample(colnames(seu), 5000)]
    } else {
      sweep_seu <- seu
    }

    sweep_res <- paramSweep(sweep_seu, PCs=1:20, sct=FALSE)
    sweep_stats <- summarizeSweep(sweep_res, GT=FALSE)
    bcmvn <- find.pK(sweep_stats)
    optimal_pK <- as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)]))
  }, error=function(e) {
    log_msg(paste("  paramSweep失败，使用默认pK=0.09:", e$message))
    optimal_pK <<- 0.09
  })

  # 估计双细胞数量（基于 10x 双细胞率表：每1000 cells约0.8%）
  n_cells <- ncol(seu)
  doublet_rate <- min(0.1, n_cells / 100000)  # 上限10%
  nExp_poi <- round(doublet_rate * n_cells)

  log_msg(paste("  ", smp, ": pK=", optimal_pK, ", nExp=", nExp_poi))

  # 运行 DoubletFinder
  seu <- doubletFinder(
    seu,
    PCs    = 1:20,
    pN     = DF_pN,
    pK     = optimal_pK,
    nExp   = nExp_poi,
    reuse.pANN = FALSE,
    sct = FALSE
  )

  # 获取 doublet 分类列名
  df_col <- grep("DF.classifications", colnames(seu@meta.data), value=TRUE)
  pann_col <- grep("pANN", colnames(seu@meta.data), value=TRUE)

  n_doublets   <- sum(seu@meta.data[[df_col]] == "Doublet")
  n_singlets   <- sum(seu@meta.data[[df_col]] == "Singlet")
  doublet_pct  <- round(n_doublets / n_cells * 100, 2)

  log_msg(paste("  ", smp, ": 双细胞=", n_doublets, "(", doublet_pct, "%)"))

  # 保存 doublet 统计
  stats_row <- data.frame(
    sample       = smp,
    patient      = patient,
    fraction     = fraction,
    n_pre_qc     = n_pre_qc,
    n_post_qc    = n_post_qc,
    n_doublets   = n_doublets,
    n_singlets   = n_singlets,
    doublet_pct  = doublet_pct,
    pK_used      = optimal_pK,
    nExp_used    = nExp_poi
  )

  # 过滤双细胞，保留单细胞
  seu <- subset(seu, cells = colnames(seu)[seu@meta.data[[df_col]] == "Singlet"])

  # 清理 DoubletFinder 列（减少内存）
  seu@meta.data <- seu@meta.data[, !colnames(seu@meta.data) %in% c(df_col, pann_col)]

  log_msg(paste("  ", smp, ": 最终保留", ncol(seu), "单细胞"))

  return(list(seurat=seu, stats=stats_row))
}

# 逐样本处理（内存限制下串行）
for (smp in analysis_samples) {
  result <- process_sample(smp)
  if (!is.null(result)) {
    seurat_list[[smp]] <- result$seurat
    doublet_stats <- rbind(doublet_stats, result$stats)
  }
}

log_msg(paste("成功处理样本数:", length(seurat_list)))

# ── Step 5: QC 汇总图 ────────────────────────────────────────────────────────
log_msg("Step 5: 绘制 QC 汇总图 ...")

# 合并所有样本的 meta.data 用于出图
all_meta <- do.call(rbind, lapply(names(seurat_list), function(smp) {
  m <- seurat_list[[smp]]@meta.data[, c("sample","patient","fraction","nFeature_RNA","nCount_RNA","pct_mt","pct_ribo")]
  m
}))

# QC violin 图
p_qc <- ggplot(all_meta, aes(x=sample, y=nFeature_RNA, fill=fraction)) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.3, fill="white") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(title="Genes per Cell after QC", x="Sample", y="Number of Genes") +
  scale_fill_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8"))

ggsave(file.path(fig_out, "qc_violin_genes.pdf"), p_qc, width=16, height=5)
ggsave(file.path(fig_out, "qc_violin_genes.png"), p_qc, width=16, height=5, dpi=150)
log_msg("  图: qc_violin_genes — 各样本QC后每细胞基因数分布，CD45p（红）与CD45n（蓝）分类展示")

p_mt <- ggplot(all_meta, aes(x=sample, y=pct_mt, fill=fraction)) +
  geom_violin(scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.3, fill="white") +
  geom_hline(yintercept=MAX_MT_PCT, linetype="dashed", color="red") +
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1, size=6)) +
  labs(title="Mitochondrial Gene Percentage after QC", x="Sample", y="MT%") +
  scale_fill_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8"))

ggsave(file.path(fig_out, "qc_violin_mt.pdf"), p_mt, width=16, height=5)
ggsave(file.path(fig_out, "qc_violin_mt.png"), p_mt, width=16, height=5, dpi=150)
log_msg("  图: qc_violin_mt — 各样本线粒体基因百分比（红虚线=25%阈值）")

# 细胞数柱状图
cell_counts <- doublet_stats
p_cells <- ggplot(cell_counts, aes(x=sample, fill=fraction)) +
  geom_col(aes(y=n_post_qc), position="dodge") +
  geom_col(aes(y=-n_doublets), fill="grey50", alpha=0.6) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7)) +
  labs(title="Cell Counts: Singlets (positive) vs Doublets (negative)",
       x="Sample", y="Cell Count") +
  scale_fill_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8"))

ggsave(file.path(fig_out, "qc_cell_counts.pdf"), p_cells, width=16, height=5)
ggsave(file.path(fig_out, "qc_cell_counts.png"), p_cells, width=16, height=5, dpi=150)
log_msg("  图: qc_cell_counts — 各样本单细胞（正值）与双细胞（负值）数量")

# Doublet 率图
p_dbr <- ggplot(doublet_stats, aes(x=sample, y=doublet_pct, fill=fraction)) +
  geom_col() +
  geom_hline(yintercept=c(1,5,10), linetype="dashed", color=c("grey","orange","red")) +
  theme_bw() +
  theme(axis.text.x=element_text(angle=45, hjust=1, size=7)) +
  labs(title="Doublet Rate per Sample", x="Sample", y="Doublet Rate (%)") +
  scale_fill_manual(values=c("CD45p"="#E41A1C","CD45n"="#377EB8"))

ggsave(file.path(fig_out, "qc_doublet_rate.pdf"), p_dbr, width=16, height=5)
ggsave(file.path(fig_out, "qc_doublet_rate.png"), p_dbr, width=16, height=5, dpi=150)
log_msg("  图: qc_doublet_rate — 各样本双细胞率（参考线：1%/5%/10%）")

# ── Step 6: 保存结果 ─────────────────────────────────────────────────────────
log_msg("Step 6: 保存结果 ...")

# 保存每个样本的 Seurat 对象列表
saveRDS(seurat_list, file.path(qc_out, "seurat_per_sample_postQC.rds"))
write.csv(doublet_stats, file.path(doublet_out, "doublet_stats.csv"), row.names=FALSE)
write.csv(all_meta, file.path(qc_out, "qc_meta_all_cells.csv"), row.names=TRUE)

# ── 自检 ─────────────────────────────────────────────────────────────────────
log_msg("=== 自检 Module 1a ===")
total_cells <- sum(doublet_stats$n_singlets)
log_msg(paste("总单细胞数:", total_cells, "(预期: 80,000-130,000)"))
if (total_cells < 80000) log_msg("⚠ 警告: 总细胞数<80,000，QC可能过严")
if (total_cells > 130000) log_msg("⚠ 警告: 总细胞数>130,000，需检查双细胞剔除")

min_cells_per_sample <- min(doublet_stats$n_singlets)
log_msg(paste("最小单样本细胞数:", min_cells_per_sample, "(预期: ≥500)"))
if (min_cells_per_sample < 500) {
  low_samp <- doublet_stats$sample[doublet_stats$n_singlets < 500]
  log_msg(paste("⚠ 以下样本细胞数<500:", paste(low_samp, collapse=", ")))
}

max_dbr <- max(doublet_stats$doublet_pct)
log_msg(paste("最高双细胞率:", max_dbr, "%(预期: <15%)"))
if (max_dbr > 15) log_msg("⚠ 警告: 部分样本双细胞率>15%")

log_msg("=== Module 1a 完成 ===")
log_msg(paste("输出: results/01_scRNAseq/QC/seurat_per_sample_postQC.rds"))
log_msg(paste("图片: figures/module01/qc_*.pdf/png"))
