#!/usr/bin/Rscript

## This script is used for generate genomic profiles for CLIP peaks in 5utr, cds, 3utr. Profiles from multiple
## CLIP peak files are to be combined into a matrix for heatmap plotting
## Author: Shuye Pu
## date created: May 02, 2022

source("/home/greenblattlab/shuyepu/Rscripts_on_diamond/Rscripts/GenomicPlot/R/GenomicPlot.R")

gtffile <- "/home/greenblattlab/shuyepu/genomic_feature/gencode.v19.annotation.gtf"

if(file.exists("/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")){
	txdb <- AnnotationDbi::loadDb("/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")
	#print(class(txdb))
}else{
	txdb <-  makeTxDbFromGFF(gtffile)
	AnnotationDbi::saveDb(txdb, file="/home/greenblattlab/shuyepu/genomic_feature/txdb.sql")
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
df <- plot_start_end_feature(queryfiles=queryfile, 
				txdb=txdb, 
				inputfiles=NULL,
				featureName="intron",
				binsize=10,
				longest=TRUE,
				insert=500,
				ext=c(-100,500,-500,100),
				hl=c(-0,0,-0,0),
				randomize=FALSE, 
				norm=FALSE,
                                fix_width=21, 
				CLIP_reads=FALSE,
                                smo=TRUE, 
				scale=FALSE, 
				stranded=TRUE, 
				outPrefix=NULL, 
				genome="hg19",
                                heatmap=FALSE, 
				rm.outlier=FALSE, 
				useScore=FALSE, 
				useSizeFactor=FALSE)
#print(summary(df))
write.table(df, file.path(outdir, paste0(querylabel, "_plot_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)


## END
