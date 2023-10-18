## This script has to be called after "run_plot_5parts_metagene.R" and "run_plot_intron_junctions.R"

library(ComplexHeatmap)
library(circlize)
library(RColorBrewer)

args <- commandArgs(trailingOnly=TRUE)
wd <- args[1] ## "/home/greenblattlab/shuyepu/Nabeel/ncRBP"
setwd(wd)
outd <- args[2]
filePattern <- args[3] #"3parts", "intron", or "m6A", "R_loop"
featureTypes <- unlist(strsplit(args[4], split=",", fixed=TRUE)) # "utr5", "cds", "utr3" or "none"
peak <- args[5]  # merged_CITS
all <- args[6] # enter "all" or "selected"
limits <- as.numeric(unlist(strsplit(args[7], split=",", fixed=T))) # for m6A or R_loop, up stream and downstrem boundary, eg. -100,100
if(length(limits)==2){
	uplimit <- limits[1]
	downlimit <- limits[2]
	print(c(uplimit, downlimit))
}

print(filePattern)
print(featureTypes)
print(peak)
print(all)

dat_files <- list.files(pattern=paste0("_", filePattern, "_plot_df.tab"))
print(length(dat_files))

intensity_dat <- lapply(dat_files, function(f){
	df <- as.data.frame(read.table(f, header=TRUE, sep="\t"))
})
#print(intensity_dat[[1]])
intensity_mat <- lapply(intensity_dat, function(x){return(x$Intensity)})
intensity_mat <- as.data.frame(do.call(rbind, intensity_mat))
#print(intensity_mat)
genes <- sapply(intensity_dat, function(x){unique(gsub(paste0("_",filePattern), "", x$Query))})
print("Duplicated genes")
print(genes[duplicated(genes)])
rownames(intensity_mat) <- genes
colnames(intensity_mat) <- intensity_dat[[1]]$Position
print(dim(intensity_mat))

if(all == "all"){ 
	exclude <- c("")
}else{
	exclude <- c("GFP", "m6AGFP", "m6AZBTB48", "m6AFTO", "m6AHek", "Input", "CBLL1")
}
print(paste(c("excluding:", exclude)))

intensity_mat <- intensity_mat[!rownames(intensity_mat) %in% exclude,]
intensity_mat[intensity_mat < 0] <- 0 ## remove values turned into negative by spline smoothing

## normalize matrix
#intensity_matN <- t(scale(t(intensity_mat))) ## zscore normalization
#intensity_matN <- normalize.quantiles(intensity_mat) ## quantile normalization
intensity_matN <- t(apply(intensity_mat, 1, scales::rescale)) ## scale between 0 and 1
rownames(intensity_matN) <- rownames(intensity_mat)
colnames(intensity_matN) <- intensity_dat[[1]]$Position
positions <- intensity_dat[[1]]$Position
print(dim(intensity_mat))
#print(colnames(intensity_mat))
if (filePattern %in% c("m6A", "GLORIm6A", "GLORIm6Aclust", "m6Am", "TSScgt",  "R_loop", "ChIPoverlap")){
	selected_columns <- as.character(intensity_dat[[1]]$Position[intensity_dat[[1]]$Position >= uplimit & intensity_dat[[1]]$Position <= downlimit])
	intensity_mat <- intensity_mat[, selected_columns]
	intensity_matN <- intensity_matN[, selected_columns]
	cols <- rep("#c6dcff", length(selected_columns))
	positions <- as.numeric(selected_columns)
	#print(paste(c("selected", selected_columns), collapse=" "))
}

#print(colnames(intensity_mat))
print(dim(intensity_mat))
print("plotting heatmap")
pdffile <- file.path(outd, paste0(peak, "_peak_", filePattern, "_heatmap_", all, ".pdf"))
if(filePattern %in% c("3parts", "ChIPparts")){
	features <- factor(intensity_dat[[1]]$Feature, levels=featureTypes)
}else if(filePattern == "intron"){
	features <- factor(intensity_dat[[1]]$Feature, levels=c("UE", featureTypes, "DE"))
	features[features == "Start of intron" & positions < 0] <- "UE"
	features[features == "End of intron" & positions >= 0] <- "DE"
}else if (filePattern %in% c("m6A", "GLORIm6A", "GLORIm6Aclust", "m6Am", "TSScgt", "R_loop", "ChIPoverlap")){
	features <- factor(rep(filePattern, length(selected_columns)), levels=featureTypes)
        features[positions < 0] <- "upstream"
        features[positions > 0] <- "downstream"
	pdffile <- file.path(outd, paste0(peak, "_peak_", filePattern, "_heatmap_", all, "_", downlimit, ".pdf"))
}else{
	stop("file pattern is not supported")
}

print("preparing horizontal annotation")
cols <- rep("#c6dcff", length(features))
names(cols) <- features
#print(cols)
ha <- HeatmapAnnotation(df = data.frame(feature = features), col=list(feature=cols), which="column", show_legend=FALSE)

vannotFile <- "/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_list.txt"
vdf <- read.delim(vannotFile, header=T, stringsAsFactors=F)
print(paste(c("annotation", dim(vdf))))

vannot <- vdf$Annotation
names(vannot) <- vdf$Gene

print("Genes in annotation table but not in data:")
print(setdiff(vdf$Gene, rownames(intensity_mat)))
print("Genes in data but not in annotation table:")
print(setdiff(rownames(intensity_mat), vdf$Gene))

vannot <- vannot[rownames(intensity_mat)]  ## Genes in row names get annotations, genes not in vdf$Gene will get NA values as annotation
names(vannot) <- rownames(intensity_mat)
vannot[is.na(vannot)] <- "Other"

vcolors <- brewer.pal(length(unique(vdf$Annotation)), name = "Set1")
names(vcolors) <- sort(unique(vdf$Annotation))
vcols <- vcolors[vannot] ## annotations get colors, NA annotations get NA values as color
names(vcols) <- vannot

print("preparing vertical annotation")
#print(vcols)
va <- HeatmapAnnotation(df = data.frame(Gene=vannot), col=list(Gene=vcols), which="row", na_col="grey")

	limmin <- min(intensity_mat)
        limmax <- max(intensity_mat)
        limmed <- (limmax+limmin)/2

        print(c(limmin, limmed, limmax))


intensity_list <- list("Raw"=intensity_mat, "Scaled"=log(intensity_mat*1000+1), "Normalized"=intensity_matN)
pdf(pdffile, width=8, height=8)
lapply(names(intensity_list), function(x){
	intensity_x <- intensity_list[[x]]
	limmin <- min(intensity_x)
	limmax <- max(intensity_x)
	limmed <- (limmax+limmin)/2

	print(c(limmin, limmed, limmax))
	h <- Heatmap(intensity_x, 
		name=paste(x,"Density"), 
		col = colorRamp2(c(min(0, limmin), limmed, limmax), c("#3e26a7", "#13beb7", "#f9fb15")), 
		top_annotation = ha,
		left_annotation = va, 
		clustering_distance_rows = "euclidean",
		clustering_method_rows = "average",
		show_row_names = TRUE, 	
		show_column_names = FALSE,
		#column_names_side = "bottom", 
		cluster_columns=FALSE,
	        row_dend_reorder = TRUE,	
		row_names_side="left", 	
		row_names_gp = gpar(fontsize = 7), 	
		column_split = features,
		column_gap = unit(0.5, "mm"))
	#draw(h)
})
dev.off()

print(paste("finish", filePattern, "heatmap of CITS peaks"))





