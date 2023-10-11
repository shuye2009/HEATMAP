## This script depends on 'compute_heatmap_CITS.sh' which in turn depends on 
## "run_plot_5parts_metagene.R",
## "run_plot_intron_junctions.R",
## "run_plot_reference_locus.R",
## "run_plot_peak_annotation.R"
## This script is depended on by 'run_plot_CITS_heatmap.sh'

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggtree)
library(ggsci)

args <- commandArgs(trailingOnly=TRUE)

wd <- args[1]
outd <- args[2]
filePattern <- args[3] #"peak_annotation"
peak <- args[4]  # CITS_merged_filtered
all <- args[5] # enter "all" or "selected"
setwd(wd)

print(filePattern)
print(peak)
print(all)

print("Grouping genes")
vannotFile <- "/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/CLIP_list.txt"
vdf <- read.delim(vannotFile, header=T, stringsAsFactors=F)
print(paste(c("annotation", dim(vdf))))

vannot <- vdf$Annotation
names(vannot) <- vdf$Gene

if(all == "all"){
        exclude <- c("")
}else{
        exclude <- c("GFP", "m6AGFP", "m6AZBTB48", "m6AFTO", "m6AHek", "Input", "CBLL1")
}
print(paste(c("excluding:", exclude)))

pdffile <- file.path(outd, paste0(peak, "_peak_", filePattern, "_barchart_", all, ".pdf"))
pdf(pdffile, width=12, height=12)

lvs <- list("barchart"=c("Intron", "3'UTR", "CDS", "5'UTR"), 
	    "biotype"=c("protein_coding", "lincRNA", "rRNA", "snRNA", "snoRNA", "antisense", "pseudogene", "intergenic", "other"), 
	    "simplified"=c("3'UTR", "CDS", "5'UTR"))

for(g in c("barchart", "biotype", "simplified")){
	lv <- lvs[[g]]
	dat_files <- list.files(pattern=paste0(filePattern, "_", g, "_df.tab"))
	print(length(dat_files))
	
	genes <- sapply(dat_files, function(x){unlist(strsplit(x, split="_"))[1]})

	stat_dat <- lapply(dat_files, function(f){
		df <- as.data.frame(read.delim(f, header=TRUE, sep="\t"))
	})
	print(stat_dat[[1]])
	measures <-  c("count", "percent", "norm_count", "norm_percent")
	Ylabs <- c("Number of peaks", "Percent of peaks", "Number of peaks per Kb", "Percent of peak per Kb")
	names(Ylabs) <- measures
	if(g == "biotype")  measures <- c("count", "percent")
	for(measure in measures){
		print(paste("processing column", measure))
		stat_mat <- lapply(stat_dat, function(x){return(x[[measure]])})
		stat_mat <- as.data.frame(do.call(rbind, stat_mat))

		rownames(stat_mat) <- genes
        	colnames(stat_mat) <- stat_dat[[1]]$feature
        	stat_mat <- stat_mat[!rownames(stat_mat) %in% exclude,]
        	print(dim(stat_mat)); print(head(stat_mat))

		of <- max(apply(stat_mat, 1, sum))/100
		Ylab <- Ylabs[measure]
	
		hc <- hclust(dist(stat_mat, method="euclidean"), method="average")
	
		ncluster <- 5
		clus <- cutree(hc, ncluster)
        	g <- split(names(clus), clus)

        	p <- ggtree(hc, linetype='solid')
        	clades <- sapply(g, function(n) MRCA(p, n))

        	p <- groupClade(p, clades, group_name='subtree') # + aes(color=subtree)

        	d <- data.frame(label = names(clus), annot = vannot[names(clus)])
		d[is.na(d)] <- "Other"
        	atree <- p %<+% d +
                	#layout_dendrogram() +
                	#geom_tippoint(aes(shape=20, color=factor(annot)), size=2) +
                	geom_tiplab(aes(color=factor(annot)), size=3, offset=of, hjust=0) +
                	#scale_color_brewer(palette='Set1', breaks=1:ncluster) +
                	guides(color = guide_legend(title = "GeneType")) +
			scale_color_d3() +
                	hexpand(0.15)

		print("Genes in annotation table but not in data:"); print(setdiff(vdf$Gene, rownames(stat_mat)))
		## rotate matrix and add a column named 'Feature'
		stat_mat <- as.data.frame(t(stat_mat)) %>%
			mutate(Feature = rownames(.))
		stat_mat_long <- pivot_longer(data=stat_mat, cols = 1:(ncol(stat_mat)-1), names_to = 'Gene', values_to = 'Count')
		stat_mat_long <- stat_mat_long %>%
			mutate(Gene = as.factor(Gene), Annotation = as.factor(vannot[.data[["Gene"]]])) %>%
			mutate(Feature=factor(Feature, levels=lv)) %>%
			select(Gene, Count, Feature, Annotation)

		# Stacked bar plot with hierarchical clustering tree
		p <- facet_plot(atree, data=stat_mat_long, geom=geom_bar, panel=Ylab,
			   mapping=aes(x=Count, fill=Feature), stat="identity", color="black", orientation="y") +
			   guides(fill = guide_legend(title = "Feature")) +
			   scale_fill_d3() +
			   theme_tree2()

		p <- facet_widths(p, c(Tree=0.6))
		print(p);
	}
}

dev.off()
print(paste("finish", filePattern, "barchart of CITS peaks"))





