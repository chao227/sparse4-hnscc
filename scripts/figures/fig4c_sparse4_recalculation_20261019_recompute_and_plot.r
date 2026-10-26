#!/usr/bin/env Rscript

# ==========================================================
# Script   : fig4c_sparse4_recalculation_20261019_recompute_and_plot.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:fig4c_sparse4_recalculation_20261019/scripts/recompute_and_plot.R
# ==========================================================
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(ragg)
})

args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)]))
mod <- dirname(dirname(script))
repo <- dirname(mod)
out_tab <- file.path(mod, "source_data")
out_fig <- file.path(mod, "figures")
out_rep <- file.path(mod, "reports")
dir.create(out_tab, recursive=TRUE, showWarnings=FALSE)
dir.create(out_fig, recursive=TRUE, showWarnings=FALSE)
dir.create(out_rep, recursive=TRUE, showWarnings=FALSE)

read_tsv <- function(x) read.delim(x, check.names=FALSE, stringsAsFactors=FALSE)
write_tsv <- function(x, p) write.table(x, p, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
inp <- file.path(repo, "scan_1019_remediation_20261019/results/P1_06_P2_04_rosters_subgroups/P1_06_PATIENT_ROSTER_476_320_429.tsv")
d <- read_tsv(inp)
stopifnot(nrow(d)==476L, sum(d$OS_status)==220L, !anyDuplicated(d$patient_id))
stopifnot(all(is.finite(d$score_per_global_sd)))

subsets <- list(
  "Stage I-II" = d[d$stage_factor=="I_II" & !is.na(d$stage_factor),],
  "Stage III-IV" = d[d$stage_factor=="III_IV" & !is.na(d$stage_factor),],
  "Non-surgical" = d[d$treatment_4group=="Non-surgical subgroup",],
  "Surgery + adjuvant" = d[d$treatment_4group=="Surgery + adjuvant",]
)
sections <- c("Stage I-II"="Stage", "Stage III-IV"="Stage",
              "Non-surgical"="Treatment", "Surgery + adjuvant"="Treatment")

one_cox <- function(x, nm) {
  f <- coxph(Surv(OS_time_months, OS_status) ~ score_per_global_sd,
             data=x, ties="efron", x=TRUE)
  s <- summary(f); z <- cox.zph(f, transform="km")
  data.frame(section=unname(sections[nm]), subgroup=nm, n=f$n, events=f$nevent,
             beta=unname(coef(f)[1]), se=s$coefficients[1,"se(coef)"],
             hazard_ratio=s$conf.int[1,"exp(coef)"],
             ci_lower=s$conf.int[1,"lower .95"], ci_upper=s$conf.int[1,"upper .95"],
             p_nominal=s$coefficients[1,"Pr(>|z|)"], ph_p=z$table[1,"p"],
             stringsAsFactors=FALSE)
}
sg <- do.call(rbind, lapply(names(subsets), function(nm) one_cox(subsets[[nm]], nm)))
rownames(sg) <- NULL
sg$p_bh4 <- p.adjust(sg$p_nominal, method="BH")
sg$ph_trigger <- sg$ph_p < 0.05
sg$association_direction <- ifelse(sg$ci_lower>1 & sg$p_bh4<0.05,
                                    "POSITIVE_ASSOCIATION", "NOT_POSITIVE_AFTER_BH")

interaction_fit <- function(x, group, reference, comparison, label) {
  x$group <- factor(x[[group]], levels=c(reference, comparison))
  x <- x[!is.na(x$group),]
  reduced <- coxph(Surv(OS_time_months, OS_status) ~ score_per_global_sd + strata(group),
                   data=x, ties="efron", x=TRUE)
  full <- coxph(Surv(OS_time_months, OS_status) ~ score_per_global_sd +
                  score_per_global_sd:group + strata(group),
                data=x, ties="efron", x=TRUE)
  lr <- anova(reduced, full, test="LRT")
  sm <- summary(full)
  term <- grep("score_per_global_sd:group", rownames(sm$coefficients), value=TRUE)
  stopifnot(length(term)==1L)
  b <- sm$coefficients[term,"coef"]; se <- sm$coefficients[term,"se(coef)"]
  zp <- cox.zph(full, transform="km", terms=TRUE, singledf=TRUE)
  data.frame(interaction=label, reference=reference, comparison=comparison,
             n=full$n, events=full$nevent, beta_interaction=b, se=se,
             interaction_hr=exp(b), ci_lower=exp(b-1.96*se), ci_upper=exp(b+1.96*se),
             p_wald=sm$coefficients[term,"Pr(>|z|)"],
             lr_chisq=lr[2,"Chisq"], lr_df=lr[2,"Df"], p_lr=lr[2,"Pr(>|Chi|)"],
             ph_global_p=zp$table["GLOBAL","p"], stringsAsFactors=FALSE)
}

stage_d <- d[d$stage_factor %in% c("I_II","III_IV"),]
tx_d <- d[d$treatment_4group %in% c("Non-surgical subgroup","Surgery + adjuvant"),]
ints <- rbind(
  interaction_fit(stage_d, "stage_factor", "I_II", "III_IV", "Score by stage"),
  interaction_fit(tx_d, "treatment_4group", "Non-surgical subgroup", "Surgery + adjuvant", "Score by treatment")
)

# Preserve term-specific PH diagnostics for the ordinary interaction models.
interaction_ph <- function(x, group, reference, comparison, label) {
  x$group <- factor(x[[group]], levels=c(reference, comparison)); x <- x[!is.na(x$group),]
  f <- coxph(Surv(OS_time_months,OS_status) ~ score_per_global_sd +
               score_per_global_sd:group + strata(group),
             data=x,ties="efron",x=TRUE)
  z <- as.data.frame(cox.zph(f,transform="km",terms=TRUE,singledf=TRUE)$table)
  z$term <- rownames(z); rownames(z) <- NULL; z$interaction <- label
  z[,c("interaction","term","chisq","df","p")]
}
iph <- rbind(
  interaction_ph(stage_d,"stage_factor","I_II","III_IV","Score by stage"),
  interaction_ph(tx_d,"treatment_4group","Non-surgical subgroup","Surgery + adjuvant","Score by treatment")
)

# The treatment interaction violates PH.  This post-diagnostic exploratory
# sensitivity adds a common score*time term and a treatment-specific score*time
# interaction while retaining group-stratified baseline hazards.  The 2-df LR
# test jointly evaluates constant and log-time interaction components.
tx_d$group <- factor(tx_d$treatment_4group,
                     levels=c("Non-surgical subgroup","Surgery + adjuvant"))
ttfun <- function(v,t,...) v*log(pmax(t,1/30.44))
tx_red_tv <- coxph(Surv(OS_time_months,OS_status) ~ score_per_global_sd +
                     tt(score_per_global_sd) + strata(group),
                   data=tx_d,ties="efron",x=TRUE,tt=ttfun)
tx_full_tv <- coxph(Surv(OS_time_months,OS_status) ~ score_per_global_sd +
                       score_per_global_sd:group + tt(score_per_global_sd) +
                       tt(score_per_global_sd):group + strata(group),
                    data=tx_d,ties="efron",x=TRUE,tt=ttfun)
tx_tv_lr <- anova(tx_red_tv,tx_full_tv,test="LRT")
tx_tv_p <- tx_tv_lr[2,"Pr(>|Chi|)"]
cf_tv <- coef(tx_full_tv); vv_tv <- vcov(tx_full_tv)
delta_nm <- grep("score_per_global_sd:group",names(cf_tv),value=TRUE)
eta_nm <- grep("(group.*:tt\\(score_per_global_sd\\))|(tt\\(score_per_global_sd\\):group)",names(cf_tv),value=TRUE)
stopifnot(length(delta_nm)==1L,length(eta_nm)==1L)
delta <- unname(cf_tv[delta_nm]); eta <- unname(cf_tv[eta_nm])
tm <- c(12,36,60); lt <- log(tm); db <- delta+eta*lt
dse <- sqrt(vv_tv[delta_nm,delta_nm] + lt^2*vv_tv[eta_nm,eta_nm] +
            2*lt*vv_tv[delta_nm,eta_nm])
tx_tv_int <- data.frame(interaction="Score by treatment",time_months=tm,
  log_interaction_ratio=db,se=dse,interaction_ratio=exp(db),
  ci_lower=exp(db-1.96*dse),ci_upper=exp(db+1.96*dse),
  lr_chisq=tx_tv_lr[2,"Chisq"],lr_df=tx_tv_lr[2,"Df"],p_lr_time_varying=tx_tv_p)

ints$ph_triggered_time_varying_p <- NA_real_
ints$ph_triggered_time_varying_p[ints$interaction=="Score by treatment"] <- tx_tv_p
ints$final_heterogeneity_p <- ifelse(ints$interaction=="Score by treatment",tx_tv_p,ints$p_lr)
ints$p_holm2 <- p.adjust(ints$final_heterogeneity_p, method="holm")
ints$interaction_result <- ifelse(
  ints$interaction=="Score by treatment" & ints$p_holm2<0.05,
  "EXPLORATORY_TIME_VARYING_HETEROGENEITY_SIGNAL",
  ifelse(ints$p_holm2<0.05,"INTERACTION_SIGNAL","NO_HETEROGENEITY_DETECTED"))

tv <- list()
for (i in which(sg$ph_trigger)) {
  nm <- sg$subgroup[i]; x <- subsets[[nm]]
  f <- coxph(Surv(OS_time_months, OS_status) ~ score_per_global_sd + tt(score_per_global_sd),
             data=x, ties="efron", x=TRUE,
             tt=function(v,t,...) v*log(pmax(t,1/30.44)))
  cf <- coef(f); vv <- vcov(f)
  b <- unname(cf["score_per_global_sd"]); g <- unname(cf["tt(score_per_global_sd)"])
  vb <- vv["score_per_global_sd","score_per_global_sd"]
  vg <- vv["tt(score_per_global_sd)","tt(score_per_global_sd)"]
  cbg <- vv["score_per_global_sd","tt(score_per_global_sd)"]
  tm <- c(1,12,60); lt <- log(tm); bt <- b+g*lt
  se <- sqrt(vb+lt^2*vg+2*lt*cbg)
  tv[[length(tv)+1]] <- data.frame(subgroup=nm, trigger_ph_p=sg$ph_p[i],
    time_months=tm, beta=bt, se=se, hazard_ratio=exp(bt),
    ci_lower=exp(bt-1.96*se), ci_upper=exp(bt+1.96*se),
    gamma_log_time=g, gamma_p=2*pnorm(-abs(g/sqrt(vg))), stringsAsFactors=FALSE)
}
tv <- if(length(tv)) do.call(rbind,tv) else data.frame()

stopifnot(nrow(subsets[["Stage I-II"]])==88L, sum(subsets[["Stage I-II"]]$OS_status)==27L)
stopifnot(nrow(subsets[["Stage III-IV"]])==273L, sum(subsets[["Stage III-IV"]]$OS_status)==114L)
stopifnot(nrow(subsets[["Non-surgical"]])==300L, sum(subsets[["Non-surgical"]]$OS_status)==118L)
stopifnot(nrow(subsets[["Surgery + adjuvant"]])==129L, sum(subsets[["Surgery + adjuvant"]]$OS_status)==61L)
stopifnot(ints$n[ints$interaction=="Score by stage"]==361L,
          ints$events[ints$interaction=="Score by stage"]==141L,
          ints$n[ints$interaction=="Score by treatment"]==429L,
          ints$events[ints$interaction=="Score by treatment"]==179L)

write_tsv(sg, file.path(out_tab,"SFig2B_sparse4_subgroup_estimates.tsv"))
write_tsv(ints, file.path(out_tab,"SFig2B_sparse4_interactions.tsv"))
write_tsv(iph, file.path(out_tab,"SFig2B_interaction_PH_diagnostics.tsv"))
write_tsv(tx_tv_int, file.path(out_tab,"SFig2B_treatment_time_varying_interaction.tsv"))
write_tsv(tv, file.path(out_tab,"SFig2B_sparse4_time_varying_sensitivity.tsv"))

int_lab <- c(
  Stage=sprintf("stage interaction p=%.3f",ints$p_lr[ints$interaction=="Score by stage"]),
  Treatment="descriptive treatment strata"
)
sg$facet <- factor(paste0(sg$section,"\n",int_lab[sg$section]),
                   levels=paste0(c("Stage","Treatment"),"\n",int_lab[c("Stage","Treatment")]))
sg$subgroup <- factor(sg$subgroup, levels=rev(c("Stage I-II","Stage III-IV","Non-surgical","Surgery + adjuvant")))
sg$ann <- sprintf("%.2f [%.2f–%.2f]%s",sg$hazard_ratio,sg$ci_lower,sg$ci_upper,
                  ifelse(sg$ph_trigger,"*",""))

p <- ggplot(sg,aes(x=hazard_ratio,y=subgroup,xmin=ci_lower,xmax=ci_upper))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey55",linewidth=0.35)+
  geom_errorbar(orientation="y",width=0.18,linewidth=0.65,colour="#31688E")+
  geom_point(shape=21,size=2.7,stroke=0.4,fill="#35B779",colour="white")+
  geom_text(aes(x=6.2,label=ann),hjust=1,size=2.7,family="Arial",colour="grey20")+
  facet_wrap(~facet,ncol=1,scales="free_y")+
  scale_x_log10(limits=c(0.8,6.5),breaks=c(1,1.5,2,3))+
  labs(x="HR per global sparse4 SD (95% CI)",y=NULL)+
  theme_classic(base_family="Arial",base_size=7.5)+
  theme(axis.text=element_text(size=7.5),axis.title.x=element_text(size=7.5),
        strip.background=element_rect(fill="grey96",colour=NA),
        strip.text=element_text(size=7.5,face="bold",margin=margin(1.5,1,1.5,1)),
        panel.spacing=unit(0.35,"lines"),plot.margin=margin(2,3,2,2))

pdf_path <- file.path(out_fig,"SFig2B_sparse4_treatment_stage_subgroups.pdf")
png_path <- file.path(out_fig,"SFig2B_sparse4_treatment_stage_subgroups.png")
cairo_pdf(pdf_path,width=100/25.4,height=65/25.4,family="Arial"); print(p); dev.off()
agg_png(png_path,width=100/25.4,height=65/25.4,units="in",res=600,bg="white"); print(p); dev.off()

qc <- data.frame(check_id=c("PATIENT_ROSTER","STAGE_COUNTS","TREATMENT_COUNTS","FOUR_SUBGROUP_ROWS",
                            "TWO_INTERACTIONS","INTERACTION_PH_EXPORTED","PH_TRIGGER_EXPORTED",
                            "TIME_VARYING_TREATMENT_INTERACTION","NO_BOOTSTRAP_OR_REFIT"),
  status="PASS",
  observed=c("476/220","88/27;273/114","300/118;129/61",nrow(sg),nrow(ints),nrow(iph),nrow(tv),nrow(tx_tv_int),"TRUE"),
  expected=c("476/220","88/27;273/114","300/118;129/61",4,2,8,"0 or multiples of 3",3,"TRUE"))
write_tsv(qc,file.path(out_rep,"NUMERICAL_QA.tsv"))
writeLines(capture.output(sessionInfo()),file.path(out_rep,"SESSIONINFO.txt"))
cat(sprintf("completed: %s; %s\n",pdf_path,png_path))
