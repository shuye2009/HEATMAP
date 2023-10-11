#!/usr/bin/Rscript

## This script is used for generate genomic profiles for CLIP peaks in 5utr, cds, 3utr. Profiles from multiple
## CLIP peak files are to be combined into a matrix for heatmap plotting
## Author: Shuye Pu
## date created: May 02, 2022

#source("/home/greenblattlab/shuyepu/Rscripts_on_diamond/Rscripts/GenomicPlot/R/GenomicPlot.R")
library(GenomicPlot)

gtffile <- "/home/greenblattlab/shuyepu/genomic_feature/gencode.v19.annotation.gtf"

if(file.exists("/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")){
	txdb <- AnnotationDbi::loadDb("/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")
	#print(class(txdb))
}else{
	txdb <-  makeTxDbFromGFF(gtffile)
	AnnotationDbi::saveDb(txdb, file="/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")
}

if(file.exists("/home/greenblattlab/shuyepu/genomic_feature/3metaFeatures.rds")){
        metaFeatures <- readRDS("/home/greenblattlab/shuyepu/genomic_feature/3metaFeatures.rds")
	#print(class(metaFeatures))
}else{
	metaFeatures <- prepare_5parts_genomic_features(txdb=txdb, longest=TRUE, meta=TRUE, nbins=100, useIntron=FALSE, fiveP=0, threeP=0)
	saveRDS(metaFeatures, file="/home/greenblattlab/shuyepu/genomic_feature/3metaFeatures.rds")
}

args <- commandArgs(trailingOnly = T)
queryfile <- args[1]
querylabel <- args[2]
outdir <- args[3]

if(length(args) < 2) stop("query file and query label are mandatory!")
if(is.na(outdir)) outdir <- "."

print(queryfile)
print(querylabel)
print(outdir)

names(queryfile) <- querylabel

handleBamparams <- list(fix_width=21, fix_point="center", useScore=FALSE, outRle=TRUE, CLIP_reads=FALSE,
norm=FALSE, useSizeFactor=FALSE, genome="hg19")

df <- plot_5parts_metagene(queryfiles=queryfile,
                           gFeatures=metaFeatures,
                           handleInputParams=handleBamparams,
                           inputfiles=NULL,
                           smooth=TRUE,
                           scale=FALSE,
                           stranded=TRUE,
                           outPrefix=file.path(outdir, paste0(querylabel, "_metagene")),
                           heatmap=TRUE,
                           rmOutlier=FALSE)

#print(summary(df))
write.table(df, file.path(outdir, paste0(querylabel, "_plot_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)

## END
