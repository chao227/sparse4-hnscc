#!/usr/bin/env Rscript

# ==========================================================
# Script   : 15_external_cohorts.R
# Purpose  : Apply the locked scoring rules to external cohorts
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/
# Paper    : Methods: post hoc external evaluation; Fig. 4G; Supplementary Tables 11-12 and 18
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:public_code_release_20261023/HNSCC_sparse4_reproducibility/analysis/08_external/run_v31_yellow_enet_external.R
# ==========================================================
suppressPackageStartupMessages({
  library(AnnotationDbi); library(hgu133plus2.db); library(illuminaHumanv4.db)
  library(GEOquery); library(survival); library(pROC); library(jsonlite)
})
options(stringsAsFactors = FALSE)
need <- function(ok, ...) if (!isTRUE(ok)) stop(sprintf(...), call. = FALSE)
sha <- function(path) sub("[[:space:]].*$", "", system2("sha256sum", normalizePath(path, mustWork = TRUE), stdout = TRUE))
atomic_tsv <- function(x, path) { dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE); q<-paste0(path,".tmp-",Sys.getpid()); write.table(x,q,sep="\t",quote=FALSE,row.names=FALSE,na="NA"); need(file.rename(q,path),"cannot install %s",path) }
atomic_json <- function(x, path) { dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE); q<-paste0(path,".tmp-",Sys.getpid()); write_json(x,q,auto_unbox=TRUE,pretty=TRUE,null="null",digits=16); need(file.rename(q,path),"cannot install %s",path) }
atomic_rds <- function(x, path) { dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE); q<-paste0(path,".tmp-",Sys.getpid()); saveRDS(x,q); need(file.rename(q,path),"cannot install %s",path) }

args <- commandArgs(trailingOnly=TRUE)
need(length(args) %in% 1:2 && args[[1L]] %in% c("freeze_coverage","evaluate"),
     "usage: freeze_coverage|evaluate [sparse]")
mode <- args[[1L]]
variant <- if (length(args) == 2L) args[[2L]] else "sparse"
need(identical(variant, "sparse"), "this release entry point permits only the current sparse4 procedure")
repo <- normalizePath(getwd(), mustWork=TRUE); need(!grepl("(^|/)reback2(/|$)",repo),"reback2 prohibited")
model_dir <- file.path(repo,"results/v3_1/P1_02_YELLOW_SPARSE_ENET_20261016")
out <- file.path(repo,"agent_work/yellow_sparse_enet_external_20261016")
prefix <- "YELLOW_SPARSE"
contract_name <- "YELLOW_SPARSE_ENET_FROZEN_CONTRACT.json"
model_id <- "V31_HISTORICAL_YELLOW_SPARSE_ENET_20261016"
dir.create(out,recursive=TRUE,showWarnings=FALSE)
coef_path <- file.path(model_dir,paste0(prefix,"_FULL_COEFFICIENTS.tsv"))
scale_path <- file.path(model_dir,paste0(prefix,"_TRAINING_SCALE.tsv"))
score_path <- file.path(model_dir,paste0(prefix,"_FULL_SCORE.tsv"))
receipt_path <- file.path(model_dir,paste0(prefix,"_FULL_FIT_RECEIPT.json"))
contract_path <- file.path(model_dir,contract_name)
identity_path <- file.path(out,"PRE_EXTERNAL_YELLOW_MODEL_IDENTITY.json")
ref_path <- file.path(out,"PRE_EXTERNAL_YELLOW_TRAINING_REFERENCE.tsv")
coverage_path <- file.path(out,"EXTERNAL_COHORT_COVERAGE_CLASSIFICATION.tsv")
ledger_path <- file.path(out,"EXTERNAL_GENE_COVERAGE.tsv")
mapped_path <- file.path(out,"EXTERNAL_MAPPED_EXPRESSION_OUTCOME_BLIND.rds")
plan_path <- file.path(out,"EXTERNAL_ANALYSIS_PLAN_FROZEN.json")

if (mode == "freeze_coverage") {
  need(!file.exists(identity_path) && !file.exists(plan_path),"identity/plan already exists; refuse overwrite")
  co <- read.delim(coef_path,check.names=FALSE); sc <- read.delim(scale_path,check.names=FALSE); tcga <- read.delim(score_path,check.names=FALSE)
  fit_receipt <- fromJSON(receipt_path,simplifyVector=FALSE)
  need(identical(fit_receipt$status,"FROZEN_BEFORE_EXTERNAL_EVALUATION_PENDING_INDEPENDENT_AUDIT"),"full fit not frozen")
  need(nrow(co)>0L && !anyDuplicated(co$gene_symbol),"coefficient identity invalid")
  z <- sc[match(co$gene_symbol,sc$gene_symbol),]; need(!anyNA(z$gene_symbol) && all(z$selected),"training scale mismatch")
  ref <- data.frame(coefficient_order=seq_len(nrow(co)),gene_symbol=co$gene_symbol,beta_vst_original_scale=co$coefficient,
                    tcga_vst_mean=z$training_mean,tcga_vst_sd=z$training_sd,
                    gamma_external_z_weight=co$coefficient*z$training_sd)
  need(all(is.finite(as.matrix(ref[,3:6]))) && all(ref$tcga_vst_sd>0),"nonfinite training reference")
  vst_path <- Sys.getenv("HNSCC_TCGA_VST_TOP5000", unset = file.path(repo, "inputs/tcga/expression/vst_top5000.rds"))
  need(file.exists(vst_path), "TCGA top-5000 VST matrix not found; set HNSCC_TCGA_VST_TOP5000")
  vst <- readRDS(vst_path)
  X <- t(vst[ref$gene_symbol,tcga$sample_id,drop=FALSE]); Z <- sweep(sweep(X,2L,ref$tcga_vst_mean,"-"),2L,ref$tcga_vst_sd,"/")
  transported <- as.numeric(Z %*% ref$gamma_external_z_weight)
  centered_original <- tcga$linear_predictor - mean(tcga$linear_predictor)
  residual <- max(abs((transported-mean(transported))-centered_original))
  need(residual < 1e-10,"TCGA transport equivalence failure %.17g",residual)
  atomic_tsv(ref,ref_path)
  identity <- list(schema="V31_HISTORICAL_YELLOW_PRE_EXTERNAL_IDENTITY_V1",
                   status="IDENTITY_AND_TRANSPORT_FROZEN_BEFORE_EXTERNAL_OUTCOME_EVALUATION_PENDING_INDEPENDENT_AUDIT",
                   created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),model_id=model_id,model_variant=variant,
                   patients=fit_receipt$patients,events=fit_receipt$events,historical_yellow_n=fit_receipt$historical_yellow_n,
                   candidate_n=fit_receipt$candidate_n,selected_n=nrow(ref),score_direction="higher score denotes higher hazard",
                   transport_rule="within-cohort per-gene z score using gamma=beta*TCGA_training_SD; no cohort pooling",
                   tcga_centered_equivalence_max_abs_residual=residual,
                   evidence=list(contract_sha256=sha(contract_path),coefficients_sha256=sha(coef_path),training_scale_sha256=sha(scale_path),training_reference_sha256=sha(ref_path)))
  atomic_json(identity,identity_path)

  read_series <- function(path) {
    con<-gzfile(path,"rt"); on.exit(close(con)); n<-0L; begin<-end<-NA_integer_; platform<-character()
    repeat { line<-readLines(con,n=1L,warn=FALSE); if(!length(line)) break; n<-n+1L
      if(startsWith(line,"!Sample_platform_id\t")) platform<-unique(gsub('^"|"$',"",strsplit(line,"\t",fixed=TRUE)[[1L]][-1L]))
      if(identical(line,"!series_matrix_table_begin")) begin<-n
      if(identical(line,"!series_matrix_table_end")){end<-n;break}
    }
    need(is.finite(begin)&&is.finite(end)&&length(platform)==1L,"series matrix markers/platform invalid")
    ex<-as.matrix(read.delim(gzfile(path),skip=begin,nrows=end-begin-2L,row.names=1L,check.names=FALSE,quote='"',comment.char="")); storage.mode(ex)<-"double"
    list(expression=ex,platform=platform)
  }
  collapse_array <- function(path,db,cohort) {
    s<-read_series(path); ex<-s$expression; ann<-AnnotationDbi::select(db,keys=ref$gene_symbol,keytype="SYMBOL",columns="PROBEID")
    ann<-unique(ann[!is.na(ann$PROBEID)&ann$PROBEID%in%rownames(ex),c("SYMBOL","PROBEID")]); vv<-apply(ex,1L,var,na.rm=TRUE)
    led<-do.call(rbind,lapply(ref$gene_symbol,function(g){ p<-sort(unique(ann$PROBEID[ann$SYMBOL==g])); p<-p[p%in%rownames(ex)]
      if(!length(p)) return(data.frame(cohort=cohort,gene_symbol=g,present_probe_n=0L,selected_feature=NA_character_,selection_variance=NA_real_,feature_sd=NA_real_,finite_complete=FALSE))
      pick<-p[order(-vv[p],p)[1L]]; val<-as.numeric(ex[pick,]); data.frame(cohort=cohort,gene_symbol=g,present_probe_n=length(p),selected_feature=pick,selection_variance=vv[pick],feature_sd=sd(val),finite_complete=all(is.finite(val))) }))
    sel<-led$selected_feature[!is.na(led$selected_feature)]; need(!anyDuplicated(sel),"%s feature maps to multiple genes",cohort)
    m<-ex[sel,,drop=FALSE]; rownames(m)<-led$gene_symbol[!is.na(led$selected_feature)]
    list(matrix=m,ledger=led,sample_ids=colnames(ex),platform=s$platform,nrow=nrow(ex),ncol=ncol(ex))
  }
  g416<-collapse_array("inputs/external/GSE41613/GSE41613_series_matrix.txt.gz",hgu133plus2.db,"GSE41613")
  g658<-collapse_array("inputs/external/GSE65858/GSE65858_series_matrix.txt.gz",illuminaHumanv4.db,"GSE65858")
  m159<-as.matrix(read.delim(gzfile("inputs/external/GSE159067/GSE159067_IHN_log2cpm_data.txt.gz"),row.names=1L,check.names=FALSE)); storage.mode(m159)<-"double"
  l159<-do.call(rbind,lapply(ref$gene_symbol,function(g){
    if(!g%in%rownames(m159)) return(data.frame(cohort="GSE159067",gene_symbol=g,present_probe_n=0L,selected_feature=NA_character_,selection_variance=NA_real_,feature_sd=NA_real_,finite_complete=FALSE))
    v<-as.numeric(m159[g,]); data.frame(cohort="GSE159067",gene_symbol=g,present_probe_n=1L,selected_feature=g,selection_variance=var(v),feature_sd=sd(v),finite_complete=all(is.finite(v))) }))
  all_led<-rbind(g416$ledger,g658$ledger,l159); all_led$eligible_feature<-!is.na(all_led$selected_feature)&all_led$finite_complete&is.finite(all_led$feature_sd)&all_led$feature_sd>0
  all_led$gamma_external_z_weight<-ref$gamma_external_z_weight[match(all_led$gene_symbol,ref$gene_symbol)]; all_led$coefficient_order<-ref$coefficient_order[match(all_led$gene_symbol,ref$gene_symbol)]
  label<-function(n) if(n==nrow(ref)) "FULL" else if(n>1L) "PARTIAL" else if(n==1L) "PROXY" else "NOT_EVALUABLE"
  cov<-do.call(rbind,lapply(c("GSE41613","GSE65858","GSE159067"),function(cn){q<-all_led[all_led$cohort==cn,]; n<-sum(q$eligible_feature)
    info<-switch(cn,GSE41613=list(platform="GPL570 Affymetrix HG-U133 Plus 2.0",n=g416$ncol,rows=g416$nrow),GSE65858=list(platform="GPL10558 Illumina HumanHT-12 v4",n=g658$ncol,rows=g658$nrow),GSE159067=list(platform="IHN targeted log2 CPM panel",n=ncol(m159),rows=nrow(m159)))
    data.frame(cohort=cn,platform=info$platform,expression_feature_n=info$rows,expression_sample_n=info$n,frozen_gene_n=nrow(ref),eligible_gene_n=n,missing_or_ineligible_n=nrow(ref)-n,coverage_class=label(n),found_genes=paste(q$gene_symbol[q$eligible_feature],collapse=","),missing_or_ineligible_genes=paste(q$gene_symbol[!q$eligible_feature],collapse=","),outcome_blind_gate="FROZEN_AFTER_MODEL_BEFORE_EXTERNAL_OUTCOME_MODEL") }))
  atomic_tsv(all_led,ledger_path); atomic_tsv(cov,coverage_path)
  atomic_rds(list(identity_sha256=sha(identity_path),training_reference_sha256=sha(ref_path),mapping_rule="exact SYMBOL; highest variance probe; lexical tie break",
                  GSE41613=list(matrix=g416$matrix,sample_ids=g416$sample_ids),GSE65858=list(matrix=g658$matrix,sample_ids=g658$sample_ids),GSE159067=list(matrix=m159[intersect(ref$gene_symbol,rownames(m159)),,drop=FALSE],sample_ids=colnames(m159))),mapped_path)
  plan<-list(schema="V31_HISTORICAL_YELLOW_EXTERNAL_PLAN_V1",status="FROZEN_BEFORE_EXTERNAL_OUTCOME_MODEL_EXECUTION",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),identity_sha256=sha(identity_path),training_reference_sha256=sha(ref_path),coverage_sha256=sha(coverage_path),mapped_expression_sha256=sha(mapped_path),primary_endpoint="OS",primary_complete_model_rule="only FULL coverage is complete-model validation",partial_rule="exploratory only; no renormalized surrogate claim")
  atomic_json(plan,plan_path); print(cov); quit(save="no",status=0L)
}

identity<-fromJSON(identity_path,simplifyVector=FALSE); plan<-fromJSON(plan_path,simplifyVector=FALSE)
need(identical(plan$status,"FROZEN_BEFORE_EXTERNAL_OUTCOME_MODEL_EXECUTION")&&identical(plan$identity_sha256,sha(identity_path))&&identical(plan$training_reference_sha256,sha(ref_path))&&identical(plan$coverage_sha256,sha(coverage_path))&&identical(plan$mapped_expression_sha256,sha(mapped_path)),"external plan binding failure")
ref<-read.delim(ref_path,check.names=FALSE); cov<-read.delim(coverage_path,check.names=FALSE); mapped<-readRDS(mapped_path); gamma<-setNames(ref$gamma_external_z_weight,ref$gene_symbol)
score_cohort<-function(cn){m<-mapped[[cn]]$matrix; genes<-rownames(m); need(length(genes)>1L&&all(genes%in%names(gamma)),"%s not evaluable",cn); z<-sweep(sweep(m,1L,rowMeans(m),"-"),1L,apply(m,1L,sd),"/"); setNames(as.numeric(crossprod(gamma[genes],z)),colnames(m))}
g416<-getGEO(filename="inputs/external/GSE41613/GSE41613_series_matrix.txt.gz",GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=FALSE); p416<-pData(g416); s416<-score_cohort("GSE41613")
vraw<-as.character(p416[["vital:ch1"]]); vn<-tolower(trimws(vraw)); amap<-c("alive"=0L,"dead -unk cause"=1L,"dead-non oc"=1L,"dead-oral ca"=1L)
d416<-data.frame(cohort="GSE41613",sample_id=rownames(p416),patient_id=as.character(p416$title),time_months=suppressWarnings(as.numeric(as.character(p416[["fu time:ch1"]]))),event=as.integer(unname(amap[vn])),score=unname(s416),endpoint_raw=vraw)
g658<-getGEO(filename="inputs/external/GSE65858/GSE65858_series_matrix.txt.gz",GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=FALSE); p658<-pData(g658); s658<-score_cohort("GSE65858")
eraw<-as.character(p658[["os_event:ch1"]]); en<-tolower(trimws(eraw)); emap<-c("false"=0L,"0"=0L,"true"=1L,"1"=1L)
d658<-data.frame(cohort="GSE65858",sample_id=rownames(p658),patient_id=as.character(p658[["ID:ch1"]]),time_months=suppressWarnings(as.numeric(as.character(p658[["os:ch1"]])))/30.44,event=as.integer(unname(emap[en])),score=unname(s658),endpoint_raw=eraw)
m159<-read.csv("inputs/external/legacy_derivatives/gse159067_ies.csv",check.names=FALSE); s159<-score_cohort("GSE159067")
d159<-data.frame(cohort="GSE159067",sample_id=m159$sample_id,patient_id=m159$sample_id,time_months=suppressWarnings(as.numeric(m159$OS_time)),event=suppressWarnings(as.integer(m159$OS_status)),score=unname(s159),endpoint_raw=as.character(m159$OS_status))
all<-rbind(d416,d658,d159); all$eligible<-is.finite(all$time_months)&all$time_months>0&all$event%in%0:1&is.finite(all$score); all$coverage_class<-cov$coverage_class[match(all$cohort,cov$cohort)]
fit_one<-function(d){z<-d[d$eligible,]; z$score_z<-as.numeric(scale(z$score)); fit<-coxph(Surv(time_months,event)~score_z,data=z,ties="efron",x=TRUE); ci<-exp(confint(fit))[1L,]; cc<-concordance(Surv(time_months,event)~score,data=z,reverse=TRUE); se<-sqrt(as.numeric(cc$var)); p<-summary(fit)$coefficients[1L,"Pr(>|z|)"]
  data.frame(cohort=z$cohort[[1L]],coverage_class=z$coverage_class[[1L]],n=nrow(z),events=sum(z$event),score_gene_n=cov$eligible_gene_n[match(z$cohort[[1L]],cov$cohort)],HR_per_score_SD=unname(exp(coef(fit))),CI95_lower=unname(ci[1]),CI95_upper=unname(ci[2]),Wald_p=p,Harrell_C=as.numeric(cc$concordance),Harrell_C_SE_IJ=se,Harrell_C_CI95_lower=max(0,as.numeric(cc$concordance)-qnorm(.975)*se),Harrell_C_CI95_upper=min(1,as.numeric(cc$concordance)+qnorm(.975)*se),PH_score_p=cox.zph(fit)$table["score_z","p"]) }
surv_res<-do.call(rbind,lapply(split(all,all$cohort),fit_one)); surv_res$claim_ceiling<-ifelse(surv_res$coverage_class=="FULL","COMPLETE_FROZEN_MODEL_EXTERNAL_EVALUATION_PENDING_INDEPENDENT_AUDIT","EXPLORATORY_PARTIAL_SCORE_ASSOCIATION_NOT_MODEL_VALIDATION")
resp<-toupper(trimws(as.character(m159$response))); responder<-as.integer(resp%in%c("CR","PR")); rocfit<-pROC::roc(responder,-unname(s159),levels=c(0,1),direction="<",quiet=TRUE); rci<-pROC::ci.auc(rocfit,method="delong")
response<-data.frame(cohort="GSE159067",coverage_class=cov$coverage_class[match("GSE159067",cov$cohort)],n=length(resp),responders=sum(responder),AUC=as.numeric(pROC::auc(rocfit)),CI95_lower=as.numeric(rci[1]),CI95_upper=as.numeric(rci[3]),claim_ceiling="EXPLORATORY_PARTIAL_SCORE_RESPONSE_DISCRIMINATION_NOT_MODEL_VALIDATION")
runout<-file.path(out,"EXTERNAL_EVALUATION"); need(!dir.exists(runout),"external evaluation already exists; refuse overwrite"); dir.create(runout)
atomic_tsv(all,file.path(runout,"EXTERNAL_SAMPLE_ENDPOINT_SCORE_LEDGER.tsv")); atomic_tsv(surv_res,file.path(runout,"EXTERNAL_SURVIVAL_RESULTS.tsv")); atomic_tsv(response,file.path(runout,"EXTERNAL_RESPONSE_RESULTS.tsv")); writeLines(capture.output(sessionInfo()),file.path(runout,"SESSIONINFO.txt"))
receipt<-list(schema="V31_HISTORICAL_YELLOW_EXTERNAL_EVALUATION_V1",status="COMPUTED_PENDING_INDEPENDENT_AUDIT",created_utc=format(Sys.time(),tz="UTC",usetz=TRUE),identity_sha256=sha(identity_path),plan_sha256=sha(plan_path),full_cohort_n=sum(cov$coverage_class=="FULL"),results=list(survival_sha256=sha(file.path(runout,"EXTERNAL_SURVIVAL_RESULTS.tsv")),response_sha256=sha(file.path(runout,"EXTERNAL_RESPONSE_RESULTS.tsv")),sample_ledger_sha256=sha(file.path(runout,"EXTERNAL_SAMPLE_ENDPOINT_SCORE_LEDGER.tsv"))))
atomic_json(receipt,file.path(runout,"EXTERNAL_EVALUATION_RECEIPT.json")); print(surv_res); print(response)
