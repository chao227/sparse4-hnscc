#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_fig4d_sfig2e_like_submitted_fig4g.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_fig4d_sfig2e_like_submitted_fig4g.R
# ==========================================================

# Style-only adaptation of the submitted Fig.4G Kaplan-Meier renderer.
# The two current panels retain their own frozen patient-level data, groups,
# endpoints and statistics; only the submitted Fig.4G layout is inherited:
# in-panel legend, in-panel log-rank/n annotation, median hv guides, a 30%
# risk table, and the same red/high versus blue/low colour semantics.

suppressPackageStartupMessages({
  library(ggplot2)
  library(survival)
  library(survminer)
  library(ragg)
})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
qa_dir <- file.path(root, "panel_layout_revision_20261020/qa")
report_dir <- file.path(root, "panel_layout_revision_20261020/reports")
dir.create(qa_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(report_dir, recursive=TRUE, showWarnings=FALSE)

pt_mm <- function(pt) pt / ggplot2::.pt
theme_submitted_fig4g <- function() {
  theme_classic(base_family="Arial", base_size=8) +
    theme(
      plot.title=element_text(size=9,face="bold",hjust=.5,lineheight=.96),
      axis.title=element_text(size=8,colour="black"),
      axis.text=element_text(size=7.5,colour="black"),
      legend.title=element_blank(),
      legend.text=element_text(size=7.5,colour="black"),
      panel.grid=element_blank(),
      plot.margin=margin(4,5,2,5)
    )
}

save_surv_pair <- function(g, stem, width_mm, height_mm) {
  dir.create(dirname(stem),recursive=TRUE,showWarnings=FALSE)
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=width_mm/25.4,height=height_mm/25.4,family="Arial",onefile=TRUE)
  print(g,newpage=FALSE); dev.off()
  ragg::agg_png(paste0(stem,".png"),width=width_mm/25.4,height=height_mm/25.4,units="in",res=600,background="white")
  print(g,newpage=FALSE); dev.off()
}

make_submitted_fig4g_style <- function(data, group_col, levels_high_low,
                                       title, x_breaks, x_limit, p_value,
                                       annotation_x, output_stem,
                                       width_mm, height_mm) {
  d <- data
  d$display_group <- factor(d[[group_col]],levels=levels_high_low)
  stopifnot(!anyNA(d$display_group),nrow(d)==476L,sum(d$OS_status)==220L)
  fit <- survfit(Surv(OS_time_months,OS_status)~display_group,data=d)
  calc_p <- pchisq(survdiff(Surv(OS_time_months,OS_status)~display_group,data=d)$chisq,df=1,lower.tail=FALSE)
  stopifnot(isTRUE(all.equal(as.numeric(calc_p),as.numeric(p_value),tolerance=1e-12)))
  p_lab <- if(p_value<0.001) sprintf("P = %.2g",p_value) else sprintf("P = %.3f",p_value)

  g <- ggsurvplot(
    fit,data=d,
    palette=c("#E64B35","#4DBBD5"),
    risk.table=TRUE,risk.table.height=.30,
    risk.table.fontsize=pt_mm(7.5),
    tables.theme=theme_submitted_fig4g()+theme(axis.title.y=element_blank()),
    pval=FALSE,
    censor.size=2.3,
    conf.int=FALSE,surv.median.line="hv",
    break.time.by=if(length(x_breaks)>1) diff(x_breaks)[1] else NULL,
    xlim=c(0,x_limit),
    xlab="Time (months)",ylab="Overall survival probability",
    title=title,
    legend.title="",
    legend.labs=levels_high_low,
    legend=c(.80,.86),
    ggtheme=theme_submitted_fig4g()
  )
  g$plot <- g$plot +
    scale_x_continuous(limits=c(0,x_limit),breaks=x_breaks,expand=expansion(mult=c(0,.01))) +
    annotate("text",x=annotation_x,y=.94,
             label=sprintf("%s (log-rank)\nn = %d",p_lab,nrow(d)),
             hjust=0,vjust=1,size=pt_mm(7.5),family="Arial",colour="grey20") +
    theme(plot.margin=margin(4,5,1,5))
  # survminer renders risk rows in reverse vertical order; provide explicit
  # labels while retaining the submitted Fig.4G colored-row convention.
  g$table <- g$table +
    scale_x_continuous(limits=c(0,x_limit),breaks=x_breaks,expand=expansion(mult=c(.07,.02))) +
    labs(title="Number at risk",x="Time (months)",y=NULL) +
    theme(
      plot.title=element_text(size=8,face="bold",hjust=0),
      axis.title.x=element_text(size=7.5,margin=margin(t=2)),
      axis.text.x=element_text(size=7.5,colour="black"),
      axis.line.x=element_line(linewidth=.3,colour="black"),
      axis.ticks.x=element_line(linewidth=.3,colour="black"),
      axis.ticks.length.x=grid::unit(1,"mm"),
      plot.margin=margin(0,5,5,5)
    )
  # survminer can install element_markdown for these labels; direct slot
  # replacement avoids an incompatible theme-element merge in recent ggplot2.
  g$table$theme$axis.text.y <- element_text(
    family="Arial",size=7.5,face="bold",colour="black"
  )
  save_surv_pair(g,output_stem,width_mm,height_mm)
  list(fit=fit,p=calc_p)
}

# Fig.4D: current sparse4 median split; High is the adverse red group.
fig4d <- read.delim(file.path(root,"replacement_figures_20261018/source_data/Fig2H_Bcell_score_KM_display_data.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
fig4d_effect <- read.delim(file.path(root,"agent_work/sparse4_tcga_clinical_20261016/MEDIAN_KM_DESCRIPTIVE_EFFECT.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(as.integer(table(factor(fig4d$score_median_group,levels=c("High","Low")))),c(238L,238L)))
fit4d <- make_submitted_fig4g_style(
  fig4d,"score_median_group",c("High","Low"),
  "B-cell score and overall survival",seq(0,120,24),120,
  fig4d_effect$logrank_p,8,
  file.path(root,"replacement_figures_20261018/figures/Fig2H_sparse4_KM"),
  70,82
)

# SFig.2E: current traceable TLS21 median split; High is protective but keeps
# the submitted Fig.4G red/high and blue/low colour mapping.
sfig2e <- read.delim(file.path(root,"panel_layout_revision_20261020/source_data/SFIG2E_TLS21_SURVIVAL_SOURCE.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
sfig2e_summary <- read.delim(file.path(root,"panel_layout_revision_20261020/source_data/SFIG2E_TLS21_SURVIVAL_SUMMARY.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(identical(as.integer(table(factor(sfig2e$TLS21_group,levels=c("TLS-High","TLS-Low")))),c(238L,238L)))
fit2e <- make_submitted_fig4g_style(
  sfig2e,"TLS21_group",c("TLS-High","TLS-Low"),
  "TLS21 score and overall survival",seq(0,200,50),215,
  sfig2e_summary$logrank_p,8,
  file.path(root,"panel_layout_revision_20261020/figures/rebuilt_raster_supplementary/SFig2E_TLS21_survival_vector"),
  82,92
)

risk_check <- function(fit,times,panel,expected) {
  s <- summary(fit,times=times,extend=TRUE)
  observed <- as.integer(s$n.risk)
  stopifnot(identical(observed,as.integer(expected)))
  data.frame(panel=panel,time=s$time,stratum=sub("display_group=","",s$strata),n_risk=observed)
}
qa <- rbind(
  risk_check(fit4d$fit,seq(0,120,24),"Fig4D",c(238,93,39,16,6,3,238,115,45,16,9,7)),
  risk_check(fit2e$fit,seq(0,200,50),"SFig2E",c(238,40,8,5,1,238,35,7,2,0))
)
write.table(qa,file.path(qa_dir,"FIG4D_SFIG2E_SUBMITTED_FIG4G_STYLE_RISK_QA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(c(
  "Fig4D and SFig2E submitted-Fig4G style adaptation",
  "Reuse level: style-only inheritance.",
  "Retained: each panel's frozen patient roster, endpoint, groups, log-rank result and risk counts.",
  "Inherited: in-panel legend, log-rank/n annotation, median hv guides, 30% risk table and high-red/low-blue palette.",
  "No model fitting, cutpoint optimization, bootstrap or patient exclusion was performed."
),file.path(report_dir,"FIG4D_SFIG2E_SUBMITTED_FIG4G_STYLE_REPORT.txt"))
cat("Fig4D and SFig2E submitted-Fig4G-style rerender complete\n")
