#!/usr/bin/Rscript

## This script is used for generate genomic profiles for CLIP peaks around a fixed center points such as m6A sites. 
# Profiles from multiple CLIP peak files are to be combined into a matrix for heatmap plotting
## Author: Shuye Pu
## date created: May 24, 2022

#source("/home/greenblattlab/shuyepu/Rscripts_on_diamond/Rscripts/GenomicPlot/R/GenomicPlot.R")
library(GenomicPlot)

args <- commandArgs(trailingOnly = T)
queryfile <- args[1]
querylabel <- args[2]
centerfile <- args[3]
centerlabel <- args[4]
outdir <- args[5]

if(length(args) < 4) stop("named queryfile and centerfile  are mandatory!")
if(is.na(outdir)) outdir <- "."

print(queryfile)
print(querylabel)
print(centerfile)
print(centerlabel)
print(outdir)

names(queryfile) <- querylabel
names(centerfile) <- centerlabel

importParams <- setImportParams(fix_point="center", fix_width=20)

df <- plot_locus(queryFiles=queryfile, 
		centerFiles=centerfile, 
		ext=c(-500,500), 
		hl=c(0,0), 
		shade=FALSE, 
		smooth=TRUE, 
		importParams=importParams,  
		verbose=FALSE, 
		binSize=10, 	
		refPoint="center", 
		Xlab=centerlabel, 
		inputFiles=NULL, 
		stranded=TRUE,
                heatmap=FALSE, 
		scale=FALSE, 
		outPrefix=NULL, 
		rmOutlier=0, 
		nc=5, 
		statsMethod="wilcox.test")$plot
#print(summary(df))
write.table(df, file.path(outdir, paste0(querylabel, "_", centerlabel, "_plot_df.tab")), sep="\t", quote=F, col.names=T, row.names=T)


## END
