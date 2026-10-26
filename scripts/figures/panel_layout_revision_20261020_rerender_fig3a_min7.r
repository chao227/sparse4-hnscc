#!/usr/bin/env Rscript

# ==========================================================
# Script   : panel_layout_revision_20261020_rerender_fig3a_min7.r
# Purpose  : Generate or render one or more manuscript panels from frozen source tables
# Inputs   : see data/README.md and config/paths.yml
# Outputs  : results/figures/
# Paper    : See docs/figure_table_map.md
# Seed     : see config/params.yml
# Runtime  : not benchmarked in the public release
# Last run : see provenance/source_provenance.csv
# Source   : guarded-workspace:panel_layout_revision_20261020/scripts/rerender_fig3a_min7.R
# ==========================================================

# Presentation-only successor for Fig.3A. Statistical/data inputs are unchanged.
# H&E scale bars are redrawn as vector line/text at the final panel scale so
# every visible semantic label is at least 7 pt.

suppressPackageStartupMessages({library(png); library(jsonlite); library(ragg)})

root <- normalizePath(Sys.getenv("SPARSE4_REPO_ROOT", unset = getwd()), mustWork = TRUE)
input_root <- file.path(root, "inputs/spatial")
manifest_path <- file.path(root, "replacement_spatial_figures_20261019/source_data/GSE208253_sample_manifest.tsv")
coord_path <- file.path(input_root, "official/GSE208253_official_spatial_coordinates.tsv")
qc_spot_path <- file.path(root, "yellow_module_spatial_20261019/production/source_data/YELLOW_SPATIAL_SPOT_LEVEL.tsv.gz")
out_stems <- c(
  file.path(root, "replacement_spatial_figures_20261019/figures/Fig3/Fig3A_spatial_cohort_overview"),
  file.path(root, "yellow_module_spatial_20261019/render/figures/Fig3/Fig3A_spatial_cohort_overview"),
  file.path(root, "fig3_split_20261020/figures/Fig3/Fig3A_spatial_cohort_overview")
)
tile_dir <- file.path(root, "panel_layout_revision_20261020/qa/fig3a_clean_tiles")
dir.create(tile_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.delim(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
official <- read.delim(coord_path, check.names = FALSE, stringsAsFactors = FALSE)
qc_spot <- read.delim(gzfile(qc_spot_path), check.names = FALSE, stringsAsFactors = FALSE)
stopifnot(nrow(manifest) == 12L, !anyDuplicated(official[c("GSM", "barcode")]),
          nrow(qc_spot) == 23882L, !anyDuplicated(qc_spot[c("GSM", "barcode")]))

read_gz_png <- function(path) { con <- gzfile(path, "rb"); on.exit(close(con)); readPNG(readBin(con, "raw", n = 1e9)) }
read_gz_json <- function(path) { con <- gzfile(path, "rt"); on.exit(close(con)); fromJSON(paste(readLines(con, warn = FALSE), collapse = "\n")) }
make_crop <- function(img, co, sf, pad = .035) {
  x <- co$pixel_col_fullres * sf$tissue_hires_scalef
  y <- co$pixel_row_fullres * sf$tissue_hires_scalef
  xr <- range(x); yr <- range(y); px <- max(10, diff(xr) * pad); py <- max(10, diff(yr) * pad)
  list(x0 = max(0, xr[1]-px), x1 = min(dim(img)[2], xr[2]+px),
       y0 = max(0, yr[1]-py), y1 = min(dim(img)[1], yr[2]+py),
       w = dim(img)[2], h = dim(img)[1])
}
bar_info <- function(sf, crop, prefer = c(1000, 500, 250, 100)) {
  um_per_hires <- (55 / as.numeric(sf$spot_diameter_fullres)) / as.numeric(sf$tissue_hires_scalef)
  avail <- crop$x1 - crop$x0
  choice <- prefer[which(prefer / um_per_hires < avail * .28)[1]]
  if (is.na(choice)) choice <- tail(prefer, 1)
  list(um = choice, fraction = (choice / um_per_hires) / avail)
}

sections <- lapply(seq_len(nrow(manifest)), function(i) {
  gsm <- manifest$GSM[i]; sid <- manifest$section_id[i]
  stem <- paste0(gsm, "_", tolower(sid))
  img <- read_gz_png(file.path(input_root, "histology", paste0(stem, "_tissue_hires_image.png.gz")))
  sf <- read_gz_json(file.path(input_root, "histology", paste0(stem, "_scalefactors_json.json.gz")))
  co <- official[official$GSM == gsm, , drop = FALSE]
  qc <- qc_spot[qc_spot$GSM == gsm, , drop = FALSE]
  stopifnot(nrow(qc) > 0L, all(qc$barcode %in% co$barcode))
  crop <- make_crop(img, co, sf); bar <- bar_info(sf, crop)
  tile <- file.path(tile_dir, paste0(gsm, "_he_no_text.png"))
  ragg::agg_png(tile, width = 1800, height = 1400, units = "px", res = 600, background = "white")
  par(mar = c(0,0,0,0)); plot.new(); plot.window(xlim=c(crop$x0,crop$x1), ylim=c(crop$y1,crop$y0), asp=1, xaxs="i", yaxs="i")
  rasterImage(img, 0, crop$h, crop$w, 0, interpolate = FALSE); box(col="grey45", lwd=.45); dev.off()
  list(gsm=gsm, patient=manifest$patient_id[i], section=sid,
       n_qc=nrow(qc), n_official=nrow(co), tile=tile, bar=bar)
})

draw_panel <- function() {
  par(family="Arial", mar=c(0,0,0,0), xpd=NA, ps=12)
  plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1))
  text(.5,.97,"Spatial cohort overview",font=2,cex=.82)
  rect(.018,.07,.285,.90,col="#F5F7FA",border="#355C7D",lwd=.8)
  text(.151,.85,"GSE208253",font=2,cex=.68)
  text(.045,.73,"12 Visium sections",adj=0,cex=.62)
  text(.045,.66,"10 patients",adj=0,cex=.62)
  text(.045,.59,"23,882 QC-retained spots",adj=0,cex=.62)
  text(.045,.52,"0 eligible section exclusions",adj=0,cex=.62)
  text(.045,.41,"Repeat sections:",adj=0,font=2,cex=.62)
  text(.045,.34,"P7: S7 + S8",adj=0,cex=.62)
  text(.045,.28,"P8: S9 + S10",adj=0,cex=.62)
  text(.151,.14,"Scale bars are section-calibrated",cex=.62,col="grey25")
  xs <- seq(.32,.99,length.out=5); ys <- c(.87,.60,.33)
  for (i in seq_len(12)) {
    r <- ceiling(i/4); cc <- (i-1)%%4+1
    x0 <- xs[cc]; x1 <- xs[cc+1]-.008; y1 <- ys[r]; y0 <- y1-.18; z <- sections[[i]]
    rasterImage(readPNG(z$tile),x0,y0,x1,y1,interpolate=FALSE)
    rect(x0,y0,x1,y1,border=ifelse(z$patient %in% c("P7","P8"),"#C26D3A","grey55"),lwd=ifelse(z$patient %in% c("P7","P8"),1,.4))
    # Vector calibrated scale bar and label; text is 7.44 pt at cex=.62.
    bx1 <- x1 - .018*(x1-x0); bx0 <- bx1 - z$bar$fraction*(x1-x0); by <- y0 + .105*(y1-y0)
    rect(bx0-.006*(x1-x0), by-.070*(y1-y0), bx1+.006*(x1-x0), by+.095*(y1-y0), col=adjustcolor("black",.60), border=NA)
    segments(bx0,by,bx1,by,col="white",lwd=2.2)
    text((bx0+bx1)/2,by+.052*(y1-y0),if(z$bar$um>=1000) paste0(z$bar$um/1000," mm") else paste0(z$bar$um," µm"),col="white",font=2,cex=.62)
    text((x0+x1)/2,y0-.020,paste0(z$patient,"-",z$section),cex=.62)
    text((x0+x1)/2,y0-.044,paste0(format(z$n_qc,big.mark=",")," QC-retained spots"),cex=.62)
  }
}

save_one <- function(stem) {
  dir.create(dirname(stem), recursive=TRUE, showWarnings=FALSE)
  ragg::agg_png(paste0(stem,".png"),width=220/25.4,height=115/25.4,units="in",res=600,background="white"); draw_panel(); dev.off()
  grDevices::cairo_pdf(paste0(stem,".pdf"),width=220/25.4,height=115/25.4,family="Arial",onefile=TRUE); draw_panel(); dev.off()
}
save_one(out_stems[1])
for (stem in out_stems[-1]) for (ext in c("pdf","png")) stopifnot(file.copy(paste0(out_stems[1],".",ext),paste0(stem,".",ext),overwrite=TRUE))

write.table(data.frame(GSM=vapply(sections,`[[`,character(1),"gsm"), patient=vapply(sections,`[[`,character(1),"patient"), section=vapply(sections,`[[`,character(1),"section"), qc_retained_spots=vapply(sections,function(z)z$n_qc,integer(1)), official_in_tissue_spots=vapply(sections,function(z)z$n_official,integer(1)), scale_bar_um=vapply(sections,function(z)z$bar$um,numeric(1)), scale_bar_fraction=vapply(sections,function(z)z$bar$fraction,numeric(1)), text_pt=7.44), file.path(root,"panel_layout_revision_20261020/reports/FIG3A_VECTOR_SCALE_BAR_QA.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
cat("Fig3A min-7pt rerender complete\n")
