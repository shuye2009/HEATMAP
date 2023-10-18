#!/usr/bin/bash

## script for computing metagene profile matrix for plotting heatmap of ${peakType} peaks

peak=$1 # recurring or merged
peakType=$2 # crosslink or CITS
subject=$3 # type of plot
replace=$4 # replace existing file or not
factor=$5 # specific gene name of all genes
if [[ $# -lt 2 ]]; then
	echo "Enter: merged or merged_filtered or recurring for the first argument"
	echo "Enter: CITS or crosslink for the second argument"
	echo "Enter: type of plot for the third argument, valid valuses are: \
		5parts, 3parts, intron, m6A, GLORIm6A, m6Am, R_loop, TSScgt, ChIPoverlap, ChIPparts, annotation"
	echo "Enter: true or false, depending on whether you want to override the data"
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
GLORI_m6A_clustered=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/GLORI_identified_m6A_sites_in_HEK293_cells_clustered_hg19.bed
TSScgt_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/TSS_start_with_CGT.bed
m6A_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/putative_m6A_sites_from_m6AHek_CITS.bed
m6Am_sites=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits/m6AHek_CITS_crosslinkSite/putative_m6A/putative_m6Am_sites_from_m6AHek_CITS.bed
R_loop_sites=$HOME/Nabeel/R_loop/HEK293_R_loop.narrowPeak.bed
Chip_peak_dir=$HOME/Nabeel/clip_analysis/data_input/CLIP_ChIP_overlap/chip_seq_files_hg19

wd=$HOME/Nabeel/clip_analysis/data_input/CLIP_cits
outd=$HOME/Nabeel/clip_analysis/data_input/${peakType}_${peak}_heatmap
mkdir -p $outd
shopt -s extglob ## to allow using of regex in ls

case "$subject" in
	"5parts")
		rcmd=${rcmd_metagene5}
		;;
	"ChIPparts") 
		rcmd=${rcmd_metagene5}
		;;
	"3parts") 
		rcmd=${rcmd_metagene3}
		;;
        "intron") 
		rcmd=${rcmd_intron}
		;;
        "annotation") 
		rcmd=${rcmd_peakAnnotation}
		;;
	"m6A")
		rcmd=${rcmd_locus}
		sites=${m6A_sites}
		;;
	"GLORIm6A")
		rcmd=${rcmd_locus}
		sites=${GLORI_m6A_sites}
		;;
	"GLORIm6Aclust")
                rcmd=${rcmd_locus}
                sites=${GLORI_m6A_clustered}
                ;;
	"m6Am")
		rcmd=${rcmd_locus}
		sites=${m6Am_sites}
		;;
	"R_loop")
		rcmd=${rcmd_locus}
		sites=${R_loop_sites}
		;;
	"TSScgt")
		rcmd=${rcmd_locus}
		sites=${TSScgt_sites}
		;;
	"ChIPoverlap") 
		rcmd=${rcmd_locus}
		;;
	*)
		echo "Not a valid subject term, valid valuses are: \
                	5parts, 3parts, intron, m6A, GLORIm6A, m6Am, R_loop, TSScgt, ChIPoverlap, ChIPparts, annotation"
		exit 1
		;;
esac

cd $wd
pv=0.01 # default 0.01

for d in $(ls -d ${factor}_CITS_crosslinkSite); do
	factor=$(basename $d "_CITS_crosslinkSite")
	echo "entering $d"
	if [[ -f $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed ]]; then
		echo "   found combined_${peakType}_${pv}_${factor}.${peak}.bed"

		if [[ ! -f $outd/${factor}_${subject}_plot_df.tab ]] || [[ $replace == true ]]; then
                        echo "      processing $factor for $subject"
			if [[ $subject == +(ChIPparts|ChIPoverlap) ]]; then
				if [[ ! -f ${Chip_peak_dir}/${factor}_summits_hg19.bed ]]; then
					echo "             ${factor}_summits_hg19.bed does not exit!"
					exit 0
				fi
                        fi
			case "$subject" in
				+("5parts"|"3parts"|"intron"|"annotation"))
                        		submitjob -c 5 -m 20 Rscript $rcmd $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed \
									   ${factor}_${subject} \
									   $outd
					;;
				"ChIPparts")
                                        submitjob -c 5 -m 20 Rscript $rcmd ${Chip_peak_dir}/${factor}_summits_hg19.bed \
                                                                           ${factor}_${subject} \
                                                                           $outd
                                        ;;
				"ChIPoverlap")
                                        submitjob -c 5 -m 20 Rscript $rcmd $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed \
                                                                           ${factor} \
                                                                           ${Chip_peak_dir}/${factor}_summits_hg19.bed  \
                                                                           ${subject} \
                                                                           $outd
					;;
				*)
					submitjob -c 5 -m 20 Rscript $rcmd $wd/$d/combined_${peakType}_${pv}_${factor}.${peak}.bed \
									   ${factor} \
									   ${sites}  \
									   ${subject} \
									   $outd
					;;
			esac

                else
                        echo "      $factor $subject has been processed"
                fi
	else
		echo "   combined_${peakType}_${pv}_${factor}.${peak}.bed does not exist yet"
	fi
done		




## END
