#!/usr/bin/env Rscript

# ==========================================================
# Script   : scan_1019_remediation_20261019_12_render_scan1019_figures.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:scan_1019_remediation_20261019/scripts/12_render_scan1019_figures.R
# ==========================================================
suppressPackageStartupMessages({ library(ggplot2); library(scales) })
script_arg <- commandArgs(FALSE)[grep("^--file=", commandArgs(FALSE))]
root <- normalizePath(file.path(dirname(sub("^--file=", "", script_arg)), ".."))
out <- file.path(root, "output", "figures"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
read_tsv <- function(p) read.delim(p, check.names=FALSE, stringsAsFactors=FALSE)

theme_set(theme_classic(base_size=8, base_family="Arial") + theme(
  axis.line=element_line(linewidth=.35), axis.ticks=element_line(linewidth=.35),
  axis.text=element_text(size=7.5), axis.title=element_text(size=8),
  plot.title=element_text(size=8, face="bold", hjust=0),
  plot.subtitle=element_text(size=7.5), legend.position="none",
  plot.margin=margin(4,5,4,4)))
save_pair <- function(p, stem, w=90, h=70) {
  cairo_pdf(file.path(out,paste0(stem,".pdf")), width=w/25.4, height=h/25.4, family="Arial")
  print(p); dev.off()
  png(file.path(out,paste0(stem,".png")), width=w/25.4, height=h/25.4, units="in", res=600, type="cairo")
  print(p); dev.off()
}

# Fig. 4H: all 260 mapped patients, four prespecified molecular subtypes.
sub <- read_tsv(file.path(root,"results/P2_03_molecular_subtypes/P2_03_SUBTYPE_PATIENT_MAPPING.tsv"))
omni <- read_tsv(file.path(root,"results/P2_03_molecular_subtypes/P2_03_SUBTYPE_OMNIBUS.tsv"))
sub$Subtype_mRNA <- factor(sub$Subtype_mRNA, levels=c("Atypical","Basal","Classical","Mesenchymal"))
p1 <- ggplot(sub,aes(Subtype_mRNA,score_per_global_sd,fill=Subtype_mRNA))+
  geom_violin(scale="width",trim=FALSE,linewidth=.25,alpha=.75)+
  geom_boxplot(width=.18,outlier.shape=NA,fill="white",linewidth=.3)+
  scale_fill_manual(values=c("#86BBD8","#F4B183","#9BC995","#C6A0C8"))+
  labs(title="Sparse4 score across molecular subtypes",subtitle=sprintf("Kruskal–Wallis P = %.2g; n = 260",omni$p_value),x=NULL,y="Sparse4 score (global SD)")+
  theme(plot.title=element_text(hjust=0.5),axis.text.x=element_text(size=7.5,angle=25,hjust=1))
save_pair(p1,"Fig4H_sparse4_molecular_subtypes",100,75)

# Supplement: time-specific sensitivities only for PH-triggered strata.
tv <- read_tsv(file.path(root,"results/P1_06_P2_04_rosters_subgroups/P2_04_TRIGGERED_TIME_VARYING_SENSITIVITY.tsv"))
stopifnot(all(is.finite(tv$ci_lower)), all(tv$ci_lower > 0))
tv$label <- ifelse(tv$group_type=="treatment","Surgery + adjuvant","Oral cavity")
tv$time_label <- factor(paste0(tv$time_months," months"),levels=c("1 months","12 months","60 months"))
p2 <- ggplot(tv,aes(hazard_ratio,time_label,color=label))+
  geom_vline(xintercept=1,lty=2,color="grey55",linewidth=.35)+
  geom_errorbar(aes(xmin=ci_lower,xmax=ci_upper),orientation="y",width=.12,linewidth=.4,position=position_dodge(width=.35))+
  geom_point(size=1.8,position=position_dodge(width=.35))+
  scale_x_log10()+scale_color_manual(values=c("Oral cavity"="#4C78A8","Surgery + adjuvant"="#D95F5F"))+
  labs(title="PH-triggered time-varying score associations",subtitle="Patient-level Cox sensitivity; constant-HR interpretation prohibited",x="Time-specific HR per global SD (log scale)",y=NULL,color=NULL)+
  theme(legend.position="bottom",legend.text=element_text(size=6.8))
save_pair(p2,"SFig_score_time_varying_sensitivity",100,72)

# Supplement: score and TLS21/TLS20 sensitivity; points are definitions, not replicates.
tls <- read_tsv(file.path(root,"results/P2_02_P2_07_tls_marker/P2_07_SPARSE4_TLS21_ASSOCIATION.tsv"))
tls$label <- ifelse(tls$target == "TLS21","TLS21","TLS20 (without MS4A1)")
rho_col <- names(tls)[grepl("rho|spearman",names(tls),ignore.case=TRUE)][1]
tls$plot_label <- sprintf("%s  (rho = %.3f)", tls$label, tls[[rho_col]])
p3 <- ggplot(tls,aes(x=.data[[rho_col]],y=reorder(plot_label,.data[[rho_col]])))+
  geom_vline(xintercept=0,color="grey70",linewidth=.35)+geom_segment(aes(x=0,xend=.data[[rho_col]],yend=reorder(plot_label,.data[[rho_col]])),linewidth=1,color="#6A8EAE")+
  geom_point(size=2.4,color="#234F70")+
  coord_cartesian(xlim=c(-.72,.08),clip="off")+
  labs(title="Sparse4–TLS association",subtitle="Patient-level Spearman correlations (n = 476)",x="Spearman rho",y=NULL)+
  theme(axis.text=element_text(size=7.5),plot.subtitle=element_text(size=7.5))
save_pair(p3,"SFig_sparse4_TLS21_TLS20_sensitivity",115,65)

# F3: patient is the inferential unit; raw patient-state summaries are displayed.
f3 <- read_tsv(file.path(root,"inputs/IKZF1_patient_aware_AUCell.tsv"))
lrt <- read_tsv(file.path(root,"results/P3_07_F3_patient_state/P3_07_F3_GLOBAL_LRT.tsv"))
lev <- c("Naive/Stem-like","Progenitor-Tex","Intermediate-Tex","Terminal-Tex","Cytotoxic-Eff","Memory")
f3$CD8_state <- factor(f3$CD8_state,levels=lev)
p4 <- ggplot(f3,aes(CD8_state,IKZF1_AUCell_median,group=patient))+
  geom_line(color="grey75",linewidth=.25,alpha=.65)+geom_point(color="#3B6F8F",size=.75,alpha=.8)+
  stat_summary(aes(group=1),fun=median,geom="line",color="#B33A3A",linewidth=.8)+
  stat_summary(aes(group=1),fun=median,geom="point",color="#B33A3A",size=1.7)+
  labs(title="Patient-aware IKZF1 activity across CD8 states",subtitle=sprintf("Random-intercept global LRT P = %.3f; 18 patients × 6 states",lrt$p_value),x=NULL,y="Patient-state median IKZF1 AUCell")+
  theme(axis.text.x=element_text(angle=32,hjust=1))
save_pair(p4,"SFig_IKZF1_patient_aware_CD8_states",100,72)

writeLines(capture.output(sessionInfo()),file.path(root,"logs","figure_sessionInfo.txt"))
cat("Rendered four scan-1019 panels as PDF and PNG\n")
