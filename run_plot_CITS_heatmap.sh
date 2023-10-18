#!/usr/bin/bash

peak=$1 # crosslink_recurring or CITS_merged or CITS_merged_filtered
all=$2 # enter "all" or "selected"
plotPattern=$3
if [[ $# -lt 3 ]]; then
	echo "argument1: crosslink_recurring, crosslink_merged, CITS_recurring, CITS_merged or CITS_merged_filtered"
	echo "argument2: all or selected"
	echo "argument3: 3parts, intron, m6A, m6Am, TSScgt, R_loop, ChIPoverlap, ChIPparts, peak_annotation or all"
	exit 1
fi

wd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/${peak}_heatmap
outd=/home/greenblattlab/shuyepu/Nabeel/clip_analysis/data_input/${peak}_heatmap

rcmd=$HOME/HEATMAP/plot_CITS_heatmap.R

if [[ "$plotPattern" == +(3parts|all) ]]; then
	filePattern="3parts"
	featureTypes="utr5,cds,utr3"
	Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all"
fi

if [[ "$plotPattern" == +(ChIPparts|all) ]]; then
        filePattern="ChIPparts"
	featureTypes="promoter,utr5,cds,utr3,TTS"
        Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all"
fi

if [[ "$plotPattern" == +(intron|all) ]]; then
	filePattern="intron"
	featureTypes="Start of intron,Center of intron,End of intron"
	Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all"
fi

if [[ "$plotPattern" == +(m6A|all) ]]; then
	filePattern="m6A"
	featureTypes="upstream,m6A,downstream"
	Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-1000,1000"
fi

if [[ "$plotPattern" == +(GLORIm6A|all) ]]; then
        filePattern="GLORIm6A"
        featureTypes="upstream,GLORIm6A,downstream"
        Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-1000,1000"
fi

if [[ "$plotPattern" == +(GLORIm6Aclust|all) ]]; then
        filePattern="GLORIm6Aclust"
        featureTypes="upstream,GLORIm6Aclust,downstream"
        Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-1000,1000"
fi

if [[ "$plotPattern" == +(m6Am|all) ]]; then
        filePattern="m6Am"
        featureTypes="upstream,m6Am,downstream"
        Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-500,500"
fi

if [[ "$plotPattern" == +(TSScgt|all) ]]; then
        filePattern="TSScgt"
        featureTypes="upstream,TSScgt,downstream"
        Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-500,500"
fi

if [[ "$plotPattern" == +(R_loop|all) ]]; then
	filePattern="R_loop"
	featureTypes="upstream,R_loop,downstream"
	Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-1000,1000"
fi

if [[ "$plotPattern" == +(ChIPoverlap|all) ]]; then
	filePattern="ChIPoverlap"
	featureTypes="upstream,ChIPoverlap,downstream"
	Rscript $rcmd $wd $outd $filePattern "$featureTypes" $peak "$all" "-1000,1000"
fi

if [[ "$plotPattern" == +(peak_annotation|all) ]]; then
	filePattern="peak_annotation"
	rcmd=$HOME/HEATMAP/plot_CITS_barchart.R
	Rscript $rcmd $wd $outd $filePattern $peak $all
fi

## END
