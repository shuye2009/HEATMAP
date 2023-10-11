#!/usr/bin/Rscript

## This script is used for generate peak annotation stats for peaks in bed file
## stats from multiple peak files will be used to plot stacked barchart
## Author: Shuye Pu
## date created: Oct , 2022

library(GenomicPlot)

gtffile <- "/home/greenblattlab/shuyepu/genomic_feature/gencode.v19.annotation.gtf"

args <- commandArgs(trailingOnly = T)
queryfile <- args[1]
querylabel <- args[2]
outdir <- args[3]
op <- args[4]

if(is.na(op)) op <- NULL
if(length(args) < 2) stop("query file and query label are mandatory!")
if(is.na(outdir)) outdir <- "."

print(queryfile)
print(querylabel)
print(outdir)

names(queryfile) <- querylabel

handleBedparams <- setImportParams(fix_width=0, fix_point="center", useScore=FALSE, outRle=FALSE, offset=0,
norm=FALSE, useSizeFactor=FALSE, genome="hg19")

pa <- plot_peak_annotation(peakFile=queryfile, 
				gtfFile=gtffile, 
				importParams=handleBedparams, 
				fiveP=0,
				dsTSS=0,
                                threeP=0, 
				simple=FALSE,
				verbose=FALSE, 
				outPrefix=op
				nc=5)

write.table(pa$annotation, file.path(outdir, paste0(querylabel, "_biotype_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)
write.table(pa$stat, file.path(outdir, paste0(querylabel, "_barchart_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)
write.table(pa$simplified, file.path(outdir, paste0(querylabel, "_simplified_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)

## END
