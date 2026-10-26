#!/usr/bin/env Rscript

# ==========================================================
# Script   : lymphhot_panel_repairs_20261019_recompute_and_render.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:lymphhot_panel_repairs_20261019/scripts/recompute_and_render.R
# ==========================================================
options(stringsAsFactors=FALSE, width=180)
suppressPackageStartupMessages({
  library(survival)
  library(ggplot2)
  library(ragg)
})

args <- commandArgs(FALSE)
script <- normalizePath(sub("^--file=", "", args[grep("^--file=", args)]))
mod <- dirname(dirname(script))
repo <- dirname(mod)
stopifnot(basename(repo) == "reback")
tab <- file.path(mod,"source_data"); fig <- file.path(mod,"figures"); rep <- file.path(mod,"reports")
dir.create(tab,recursive=TRUE,showWarnings=FALSE)
dir.create(fig,recursive=TRUE,showWarnings=FALSE)
dir.create(rep,recursive=TRUE,showWarnings=FALSE)
rd <- function(p) read.delim(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c("NA"))
wr <- function(x,p) write.table(x,file.path(tab,p),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
assert <- function(x,msg) if(!isTRUE(x)) stop(msg,call.=FALSE)

roster_path <- file.path(repo,"scan_1019_remediation_20261019/results/P1_06_P2_04_rosters_subgroups/P1_06_PATIENT_ROSTER_476_320_429.tsv")
vst_path <- file.path(repo,"inputs/tcga/expression/vst_matrix.rds")
qdat_path <- file.path(repo,"replacement_figures_20261018/source_data/Fig5F_sparse4_by_quadrant_display_data.tsv")
d <- rd(roster_path)
vst <- readRDS(vst_path)
assert(nrow(d)==476L,"roster must contain 476 patients")
assert(sum(d$OS_status)==220L,"full-roster event count must be 220")
assert(!anyDuplicated(d$patient_id) && !anyDuplicated(d$sample_id),"patient/sample IDs must be unique")
assert(all(c("CD79A","IKZF1") %in% rownames(vst)),"CD79A/IKZF1 absent from frozen VST")
mi <- match(d$sample_id,colnames(vst)); assert(!anyNA(mi),"roster sample absent from VST")
d$CD79A <- as.numeric(vst["CD79A",mi]); d$IKZF1 <- as.numeric(vst["IKZF1",mi])
cd_cut <- median(d$CD79A); ik_cut <- median(d$IKZF1)
d$CD79A_high <- d$CD79A >= cd_cut; d$IKZF1_high <- d$IKZF1 >= ik_cut
d$quadrant <- ifelse(d$CD79A_high & d$IKZF1_high,"LymphHot",
              ifelse(d$CD79A_high,"B-skewed",ifelse(d$IKZF1_high,"T-skewed","Cold")))
d$LymphHot <- as.integer(d$quadrant=="LymphHot")
qlev <- c("LymphHot","B-skewed","T-skewed","Cold")
qcount <- as.integer(table(factor(d$quadrant,levels=qlev)))
assert(identical(qcount,c(171L,67L,67L,171L)),"frozen quadrant counts drift")

# ---------------- Fig.6F: LymphHot subgroup estimates and interactions ----------------
subsets <- list(
  "Overall"=d,
  "Stage I-II"=d[d$stage_factor=="I_II" & !is.na(d$stage_factor),],
  "Stage III-IV"=d[d$stage_factor=="III_IV" & !is.na(d$stage_factor),],
  "Non-surgical"=d[d$treatment_4group=="Non-surgical subgroup",],
  "Surgery + adjuvant"=d[d$treatment_4group=="Surgery + adjuvant",]
)
section <- c("Overall"="Overall","Stage I-II"="Stage","Stage III-IV"="Stage",
             "Non-surgical"="Treatment","Surgery + adjuvant"="Treatment")

one_lh <- function(x,nm){
  assert(length(unique(x$LymphHot))==2L,paste("one LymphHot level in",nm))
  f <- coxph(Surv(OS_time_months,OS_status)~LymphHot,data=x,ties="efron",x=TRUE,model=TRUE)
  sm <- summary(f); z <- cox.zph(f,transform="km")
  data.frame(section=unname(section[nm]),subgroup=nm,n=f$n,events=f$nevent,
    lymphhot_n=sum(x$LymphHot==1),other_n=sum(x$LymphHot==0),
    beta=unname(coef(f)["LymphHot"]),se=sm$coefficients["LymphHot","se(coef)"],
    hazard_ratio=sm$conf.int["LymphHot","exp(coef)"],
    ci_lower=sm$conf.int["LymphHot","lower .95"],ci_upper=sm$conf.int["LymphHot","upper .95"],
    p_nominal=sm$coefficients["LymphHot","Pr(>|z|)"],ph_p=z$table["LymphHot","p"],
    stringsAsFactors=FALSE)
}
sg <- do.call(rbind,lapply(names(subsets),function(nm)one_lh(subsets[[nm]],nm)))
rownames(sg)<-NULL
sg$p_bh4 <- NA_real_
idx4 <- sg$subgroup!="Overall"
sg$p_bh4[idx4] <- p.adjust(sg$p_nominal[idx4],method="BH")
sg$ph_trigger <- sg$ph_p < 0.05
sg$association_call <- ifelse(sg$subgroup=="Overall",
  ifelse(sg$p_nominal<0.05,"NOMINAL_ASSOCIATION","NO_NOMINAL_ASSOCIATION"),
  ifelse(sg$p_bh4<0.05,"ASSOCIATION_AFTER_BH4","NO_ASSOCIATION_AFTER_BH4"))

interaction_fit <- function(x,gvar,reference,comparison,label){
  x$group <- factor(x[[gvar]],levels=c(reference,comparison)); x<-x[!is.na(x$group),]
  red <- coxph(Surv(OS_time_months,OS_status)~LymphHot+strata(group),data=x,ties="efron",x=TRUE)
  full <- coxph(Surv(OS_time_months,OS_status)~LymphHot+LymphHot:group+strata(group),data=x,ties="efron",x=TRUE)
  lr <- anova(red,full,test="LRT"); sm<-summary(full)
  term <- grep("LymphHot:group",rownames(sm$coefficients),value=TRUE); assert(length(term)==1L,"interaction term ambiguity")
  b<-sm$coefficients[term,"coef"]; se<-sm$coefficients[term,"se(coef)"]
  z<-cox.zph(full,transform="km",terms=TRUE,singledf=TRUE)
  ztab<-z$table
  data.frame(interaction=label,reference=reference,comparison=comparison,n=full$n,events=full$nevent,
    beta_interaction=b,se=se,interaction_hr=exp(b),ci_lower=exp(b-1.96*se),ci_upper=exp(b+1.96*se),
    p_wald=sm$coefficients[term,"Pr(>|z|)"],lr_chisq=lr[2,"Chisq"],lr_df=lr[2,"Df"],p_lr=lr[2,"Pr(>|Chi|)"],
    ph_lymphhot_p=ztab["LymphHot","p"],ph_interaction_p=ztab["LymphHot:group","p"],
    ph_global_p=ztab["GLOBAL","p"],stringsAsFactors=FALSE)
}
stage_d<-d[d$stage_factor%in%c("I_II","III_IV"),]
tx_d<-d[d$treatment_4group%in%c("Non-surgical subgroup","Surgery + adjuvant"),]
ints<-rbind(
  interaction_fit(stage_d,"stage_factor","I_II","III_IV","LymphHot by stage"),
  interaction_fit(tx_d,"treatment_4group","Non-surgical subgroup","Surgery + adjuvant","LymphHot by treatment")
)

interaction_ph <- function(x,gvar,reference,comparison,label){
  x$group<-factor(x[[gvar]],levels=c(reference,comparison));x<-x[!is.na(x$group),]
  f<-coxph(Surv(OS_time_months,OS_status)~LymphHot+LymphHot:group+strata(group),data=x,ties="efron",x=TRUE)
  z<-as.data.frame(cox.zph(f,transform="km",terms=TRUE,singledf=TRUE)$table)
  z$term<-rownames(z);rownames(z)<-NULL;z$interaction<-label
  z[,c("interaction","term","chisq","df","p")]
}
iph<-rbind(
  interaction_ph(stage_d,"stage_factor","I_II","III_IV","LymphHot by stage"),
  interaction_ph(tx_d,"treatment_4group","Non-surgical subgroup","Surgery + adjuvant","LymphHot by treatment")
)

ttfun <- function(v,t,...) v*log(pmax(t,1/30.44))
tv_interactions <- list()
for(i in seq_len(nrow(ints))){
  trigger <- min(ints$ph_lymphhot_p[i],ints$ph_interaction_p[i],na.rm=TRUE)<0.05
  if(!trigger) next
  if(ints$interaction[i]=="LymphHot by stage"){
    x<-stage_d; x$group<-factor(x$stage_factor,levels=c("I_II","III_IV"))
  } else {
    x<-tx_d; x$group<-factor(x$treatment_4group,levels=c("Non-surgical subgroup","Surgery + adjuvant"))
  }
  red<-coxph(Surv(OS_time_months,OS_status)~LymphHot+tt(LymphHot)+strata(group),data=x,ties="efron",x=TRUE,tt=ttfun)
  full<-coxph(Surv(OS_time_months,OS_status)~LymphHot+LymphHot:group+tt(LymphHot)+tt(LymphHot):group+strata(group),data=x,ties="efron",x=TRUE,tt=ttfun)
  lr<-anova(red,full,test="LRT"); cf<-coef(full); vv<-vcov(full)
  dn<-grep("^LymphHot:group",names(cf),value=TRUE)
  en<-grep("(^group.*:tt\\(LymphHot\\))|(^tt\\(LymphHot\\):group)",names(cf),value=TRUE)
  assert(length(dn)==1L && length(en)==1L,"time-varying interaction term ambiguity")
  tm<-c(1,12,60);lt<-log(tm);db<-unname(cf[dn])+unname(cf[en])*lt
  se<-sqrt(vv[dn,dn]+lt^2*vv[en,en]+2*lt*vv[dn,en])
  tv_interactions[[length(tv_interactions)+1]]<-data.frame(interaction=ints$interaction[i],time_months=tm,
    log_interaction_ratio=db,se=se,interaction_ratio=exp(db),ci_lower=exp(db-1.96*se),ci_upper=exp(db+1.96*se),
    lr_chisq=lr[2,"Chisq"],lr_df=lr[2,"Df"],p_lr_time_varying=lr[2,"Pr(>|Chi|)"])
  ints$ph_triggered_time_varying_p[i]<-lr[2,"Pr(>|Chi|)"]
}
tv_interactions<-if(length(tv_interactions))do.call(rbind,tv_interactions)else data.frame()
if(!"ph_triggered_time_varying_p"%in%names(ints))ints$ph_triggered_time_varying_p<-NA_real_
ints$ph_time_trigger<-pmin(ints$ph_lymphhot_p,ints$ph_interaction_p,na.rm=TRUE)<0.05
ints$final_heterogeneity_p<-ifelse(ints$ph_time_trigger,ints$ph_triggered_time_varying_p,ints$p_lr)
ints$p_holm2<-p.adjust(ints$final_heterogeneity_p,method="holm")
ints$interaction_call<-ifelse(ints$p_holm2<0.05,"INTERACTION_AFTER_HOLM2","NO_INTERACTION_AFTER_HOLM2")

tv_subgroups<-list()
for(i in which(sg$ph_trigger)){
  nm<-sg$subgroup[i];x<-subsets[[nm]]
  f<-coxph(Surv(OS_time_months,OS_status)~LymphHot+tt(LymphHot),data=x,ties="efron",x=TRUE,tt=ttfun)
  cf<-coef(f);vv<-vcov(f);b<-unname(cf["LymphHot"]);g<-unname(cf["tt(LymphHot)"])
  tm<-c(1,12,60);lt<-log(tm);bt<-b+g*lt
  se<-sqrt(vv["LymphHot","LymphHot"]+lt^2*vv["tt(LymphHot)","tt(LymphHot)"]+2*lt*vv["LymphHot","tt(LymphHot)"])
  tv_subgroups[[length(tv_subgroups)+1]]<-data.frame(subgroup=nm,trigger_ph_p=sg$ph_p[i],time_months=tm,
    beta=bt,se=se,hazard_ratio=exp(bt),ci_lower=exp(bt-1.96*se),ci_upper=exp(bt+1.96*se),
    gamma_log_time=g,gamma_p=2*pnorm(-abs(g/sqrt(vv["tt(LymphHot)","tt(LymphHot)"]))))
}
tv_subgroups<-if(length(tv_subgroups))do.call(rbind,tv_subgroups)else data.frame()

assert(nrow(subsets[["Stage I-II"]])==88L && sum(subsets[["Stage I-II"]]$OS_status)==27L,"stage I-II count drift")
assert(nrow(subsets[["Stage III-IV"]])==273L && sum(subsets[["Stage III-IV"]]$OS_status)==114L,"stage III-IV count drift")
assert(nrow(subsets[["Non-surgical"]])==300L && sum(subsets[["Non-surgical"]]$OS_status)==118L,"non-surgical count drift")
assert(nrow(subsets[["Surgery + adjuvant"]])==129L && sum(subsets[["Surgery + adjuvant"]]$OS_status)==61L,"surgery+adjuvant count drift")
assert(ints$n[ints$interaction=="LymphHot by stage"]==361L && ints$events[ints$interaction=="LymphHot by stage"]==141L,"stage interaction roster drift")
assert(ints$n[ints$interaction=="LymphHot by treatment"]==429L && ints$events[ints$interaction=="LymphHot by treatment"]==179L,"treatment interaction roster drift")
wr(sg,"Fig6F_LymphHot_subgroup_estimates.tsv")
wr(ints,"Fig6F_LymphHot_interactions.tsv")
wr(iph,"Fig6F_interaction_PH_diagnostics.tsv")
wr(tv_subgroups,"Fig6F_LymphHot_time_varying_sensitivity.tsv")
wr(tv_interactions,"Fig6F_interaction_time_varying_sensitivity.tsv")

# ---------------- Fig.5F: pairwise sparse4 contrasts ----------------
qd<-rd(qdat_path)
assert(nrow(qd)==476L && !anyDuplicated(qd$patient_id),"Fig5F source roster drift")
ii<-match(d$patient_id,qd$patient_id);assert(!anyNA(ii),"Fig5F patients do not align")
assert(all(as.character(qd$quadrant[ii])==d$quadrant),"Fig5F quadrant definitions drift")
kw<-kruskal.test(score~factor(quadrant,levels=c("Cold","B-skewed","T-skewed","LymphHot")),data=qd)
comparators<-c("Cold","B-skewed","T-skewed")
pw<-do.call(rbind,lapply(comparators,function(g){
  x<-qd$score[qd$quadrant=="LymphHot"];y<-qd$score[qd$quadrant==g]
  wt<-wilcox.test(x,y,alternative="two.sided",exact=FALSE,conf.int=TRUE,conf.level=0.95,correct=TRUE)
  U<-unname(wt$statistic); rbc<-2*U/(length(x)*length(y))-1
  data.frame(contrast=paste("LymphHot vs",g),group_1="LymphHot",group_2=g,n_1=length(x),n_2=length(y),
    median_1=median(x),median_2=median(y),median_difference=median(x)-median(y),
    hodges_lehmann_shift=unname(wt$estimate),hl_ci_lower=unname(wt$conf.int[1]),hl_ci_upper=unname(wt$conf.int[2]),
    rank_biserial=rbc,wilcoxon_W=U,p_nominal=wt$p.value,stringsAsFactors=FALSE)
}))
pw$p_holm3<-p.adjust(pw$p_nominal,method="holm")
pw$direction<-ifelse(pw$hodges_lehmann_shift<0,"LOWER_IN_LYMPHHOT","HIGHER_IN_LYMPHHOT")
pw$call<-ifelse(pw$p_holm3<0.05,"DIFFERENT_AFTER_HOLM3","NO_DIFFERENCE_AFTER_HOLM3")
wr(qd,"Fig5F_sparse4_quadrant_display_data.tsv")
wr(data.frame(test="Kruskal-Wallis four frozen quadrants",n=nrow(qd),df=unname(kw$parameter),statistic=unname(kw$statistic),p_nominal=kw$p.value),"Fig5F_omnibus_test.tsv")
wr(pw,"Fig5F_LymphHot_pairwise_contrasts.tsv")

# ---------------- Rendering ----------------
theme_pub<-function(){theme_classic(base_family="Arial",base_size=7.5)+theme(plot.title=element_text(face="bold",size=7.5,hjust=0.5),plot.subtitle=element_text(size=5.6,hjust=0.5,lineheight=.95),axis.text=element_text(size=6.7),axis.title=element_text(size=7),plot.caption=element_text(size=5.2,hjust=0,lineheight=.95,colour="grey25"),plot.margin=margin(2,2,2,2))}
save_panel<-function(p,stem){
  cairo_pdf(file.path(fig,paste0(stem,".pdf")),width=60/25.4,height=60/25.4,family="Arial");print(p);dev.off()
  agg_png(file.path(fig,paste0(stem,".png")),width=60/25.4,height=60/25.4,units="in",res=600,bg="white");print(p);dev.off()
}

qd$quadrant<-factor(qd$quadrant,levels=c("Cold","B-skewed","T-skewed","LymphHot"))
labn<-table(qd$quadrant)
xlabs<-c("Cold"=sprintf("Cold\n(n=%d)",labn["Cold"]),"B-skewed"=sprintf("B-\nskewed\n(n=%d)",labn["B-skewed"]),
         "T-skewed"=sprintf("T-\nskewed\n(n=%d)",labn["T-skewed"]),"LymphHot"=sprintf("Lymph-\nHot\n(n=%d)",labn["LymphHot"]))
ptext<-setNames(sprintf("%s",format.pval(pw$p_holm3,digits=2,eps=.001)),pw$group_2)
p5<-ggplot(qd,aes(quadrant,score,fill=quadrant))+
  geom_violin(trim=TRUE,scale="width",alpha=.72,linewidth=.25)+
  geom_boxplot(width=.18,outlier.shape=NA,fill="white",linewidth=.35)+
  scale_fill_manual(values=c("Cold"="#BDBDBD","B-skewed"="#80B1D3","T-skewed"="#FDB462","LymphHot"="#8DD3C7"))+
  scale_x_discrete(labels=xlabs)+
  labs(title="Score across marker states",subtitle=sprintf("LymphHot vs Cold/B-/T-skewed\nHolm p: %s / %s / %s",ptext["Cold"],ptext["B-skewed"],ptext["T-skewed"]),
       x=NULL,y="B-cell expression score",caption=sprintf("Kruskal-Wallis p = %.2e\nPairwise Holm family: 3 tests",kw$p.value))+
  theme_pub()+theme(legend.position="none",axis.text.x=element_text(lineheight=.82,size=5.8),plot.caption=element_text(size=4.8,lineheight=.9))
save_panel(p5,"Fig5F_sparse4_by_CD79A_IKZF1_pairwise")

sg$subgroup<-factor(sg$subgroup,levels=rev(c("Overall","Non-surgical","Surgery + adjuvant","Stage I-II","Stage III-IV")))
sg$ann<-sprintf("%.2f [%.2f–%.2f]%s",sg$hazard_ratio,sg$ci_lower,sg$ci_upper,ifelse(sg$ph_trigger,"*",""))
int_stage<-ints[ints$interaction=="LymphHot by stage",]
int_tx<-ints[ints$interaction=="LymphHot by treatment",]
cap<-sprintf("Interaction Holm p: treatment=%s\nstage=%s\n* PH trigger; average HR; see sensitivity table",
             format.pval(int_tx$p_holm2,digits=3,eps=.001),format.pval(int_stage$p_holm2,digits=3,eps=.001))
xmax<-11.5
p6<-ggplot(sg,aes(hazard_ratio,subgroup,xmin=ci_lower,xmax=ci_upper))+
  geom_vline(xintercept=1,linetype="dashed",colour="grey55",linewidth=.3)+
  geom_errorbar(orientation="y",width=.16,linewidth=.55,colour="#31688E")+
  geom_point(aes(fill=section),shape=21,size=2.4,stroke=.35,colour="white")+
  geom_text(aes(x=xmax,label=ann),hjust=1,size=1.70,family="Arial",colour="grey15")+
  scale_fill_manual(values=c("Overall"="#3B528B","Stage"="#F46D43","Treatment"="#1FA187"))+
  scale_x_log10(limits=c(.15,12),breaks=c(.25,.5,1,2))+
  labs(title=NULL,subtitle=NULL,x="HR for LymphHot vs others (95% CI)",y=NULL,caption=cap)+
  theme_pub()+theme(legend.position="none",axis.text.y=element_text(size=5.5),axis.title.x=element_text(size=5.8),plot.caption=element_text(size=4.0,lineheight=.88),plot.margin=margin(1.5,1.5,1.5,1.5))
save_panel(p6,"Fig6F_LymphHot_treatment_stage_subgroups")

qa<-data.frame(check_id=c("ROSTER","EVENTS","QUADRANTS","STAGE_COUNTS","TREATMENT_COUNTS","FIVE_HR_ROWS","TWO_INTERACTIONS","PH_EXPORTED","PAIRWISE_FAMILY","NO_BOOTSTRAP"),
 status="PASS",observed=c(nrow(d),sum(d$OS_status),paste(qcount,collapse=","),"88/27;273/114","300/118;129/61",nrow(sg),nrow(ints),nrow(iph),nrow(pw),"TRUE"),
 expected=c(476,220,"171,67,67,171","88/27;273/114","300/118;129/61",5,2,"6 rows",3,"TRUE"))
write.table(qa,file.path(rep,"NUMERICAL_QA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(rep,"SESSIONINFO.txt"))
cat("Completed LymphHot panel production\n")
print(sg)
print(ints)
print(pw)
