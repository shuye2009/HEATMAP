#!/usr/bin/bash

## script for computing metagene profile matrix for plotting heatmap of ${peakType} peaks

peak=$1 # recurring or merged
peakType=$2 # crosslink or CITS
replace=$3 # replace existing file or not
factor=$4 # specific gene name of all genes
if [[ $# -lt 2 ]]; then
	echo "Enter: merged or merged_filtered or recurring for the first argument"
	echo "Enter: CITS or crosslink for the second argument"
	echo "Enter: true or type of plot for the third argument, if you want to override the data"
	echo "Enter: gene name like SP1 for specific gene, or all for all genes"
	exit 1
fi
if [[ -z $replace ]]; then replace=false; fi # default: do not replace
if [[ $factor == "all" ]]; then factor="*" ; fi

rcmd_metagene3=$HOME/HEATMAP/run_plot_3parts_metagene.R
rcmd_metagene5=$HOME/HEATMAP/run_plot_5parts_metagene.R
rcmd_intron=$HOME/HEATMAP/run_plot_intron_junctions.R
rcmd_locus=$HOME/HEATMAP/run_plot_locus.R
rcmd_peakAnnotation=$HOME/HEATMAP/run_plot_peak_annotation.R

GLORI_m6A_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/GLORI_identified_m6A_sites_in_HEK293_cells_all_hg19.bed
TSScgt_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/TSS_start_with_CGT.bed
m6A_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/putative_m6A_sites_from_m6AHek_CITS.bed
m6Am_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/putative_m6Am_sites_from_m6AHek_CITS.bed
R_loop_sites=$HOME/Nabeel/R_loop/HEK293_R_loop.narrowPeak.bed
Chip_peak_dir=$HOME/Nabeel/clip_analysis/data_input/CLIP_ChIP_overlap/chip_seq_files_hg19

wd=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits
outd=$HOME/Nabeel/clip_analysis/data_input/${peakType}_${peak}_heatmap
mkdir -p $outd

cd $wd
pv=0.01 # default 0.01
shopt -s extglob ## to allow using of regex in ls
for d in $(ls -d ${factor}_CITS_crosslinkSite); do
	factor=$(basename $d "_CITS_crosslinkSite")
	echo "entering $d"
	if [[ -f $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ]]; then
		echo "   found combined_${peakType}_${pv}_${factor}.${peak}.bed"

		if [[ ! -f $outd/${factor}_5parts_plot_df.tab ]] || [[ $replace == +(5parts|true) ]]; then

                        echo "      processing $factor for metagene5"
                        submitjob -m 10 Rscript $rcmd_metagene5 $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor}_5parts $outd
                else
                        echo "      $factor metagene has been processed"
                fi

		if [[ ! -f $outd/${factor}_3parts_plot_df.tab ]] || [[ $replace == +(3parts|true) ]]; then

			echo "      processing $factor for metagene3"
			submitjob -m 10 Rscript $rcmd_metagene3 $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor}_3parts $outd
		else
			echo "      $factor metagene has been processed"
		fi

                if [[ ! -f $outd/${factor}_intron_plot_df.tab ]] || [[ $replace == +(intron|true) ]]; then
                        echo "      processing $factor for intron junctions"
                        submitjob -m 10 Rscript $rcmd_intron $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor}_intron $outd
                else
                        echo "      $factor intron junction has been processed"
                fi

		if [[ ! -f $outd/${factor}_m6A_plot_df.tab ]] || [[ $replace == +(m6A|true) ]]; then
                        echo "      processing $factor for m6A sites"
                        submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${m6A_sites} "m6A" $outd
                else
                        echo "      $factor m6A site has been processed"
                fi


		if [[ ! -f $outd/${factor}_GLORIm6A_plot_df.tab ]] || [[ $replace == +(GLORIm6A|true) ]]; then
                        echo "      processing $factor for GLORIm6A sites"
                        submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${GLORI_m6A_sites} "GLORIm6A" $outd
                else
                        echo "      $factor GLORIm6A site has been processed"
                fi


		if [[ ! -f $outd/${factor}_m6Am_plot_df.tab ]] || [[ $replace == +(m6Am|true) ]]; then
                        echo "      processing $factor for m6Am sites"
                        submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${m6Am_sites} "m6Am" $outd
                else
                        echo "      $factor m6Am site has been processed"
                fi

		 if [[ ! -f $outd/${factor}_TSScgt_plot_df.tab ]] || [[ $replace == +(TSScgt|true) ]]; then
                        echo "      processing $factor for TSScgt sites"
                        submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${TSScgt_sites} "TSScgt" $outd
                else
                        echo "      $factor TSScgt site has been processed"
                fi

		if [[ ! -f $outd/${factor}_R_loop_plot_df.tab ]] || [[ $replace == +(R_loop|true) ]]; then
                        echo "      processing $factor for R_loop sites"
                        submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${R_loop_sites} "R_loop" $outd
                else
                        echo "      $factor R_loop site has been processed"
                fi

		if [[ ! -f $outd/${factor}_peak_annotation_barchart_df.tab ]] || [[ $replace == +(annotation|true) ]]; then
                        echo "      processing $factor for peak annotation"
                        submitjob -m 10 Rscript $rcmd_peakAnnotation $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor}_peak_annotation $outd
                else
                        echo "      $factor peak annotation has been processed"
                fi

		if [[ ! -f $outd/${factor}_ChIPoverlap_plot_df.tab ]] || [[ $replace == +(ChIPoverlap|true) ]]; then
                        echo "      processing $factor for ChIPseq peak overlap"
			if [[ -f ${Chip_peak_dir}/${factor}_summits_hg19.bed ]]; then
                        	submitjob -m 10 Rscript $rcmd_locus $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ${factor} ${Chip_peak_dir}/${factor}_summits_hg19.bed "ChIPoverlap" $outd
			else
				echo "    	   ${factor}_summits_hg19.bed does not exit!"
			fi
                else
                        echo "      $factor ChIPoverlap has been processed"
                fi

		if [[ ! -f $outd/${factor}_ChIPparts_plot_df.tab ]] || [[ $replace == +(ChIPparts|true) ]]; then
                        echo "      processing $factor for ChIPseq peak metagene"
                        if [[ -f ${Chip_peak_dir}/${factor}_summits_hg19.bed ]]; then
                                submitjob -m 10 Rscript $rcmd_metagene5 ${Chip_peak_dir}/${factor}_summits_hg19.bed ${factor}_ChIPparts $outd
                        else
                                echo "             ${factor}_summits_hg19.bed does not exit!"
                        fi
                else
                        echo "      $factor ChIP metagene has been processed"
                fi



	else
		echo "   combined_${peakType}_${pv}_${factor}.${peak}.bed does not exist yet"
	fi
done		




## END
