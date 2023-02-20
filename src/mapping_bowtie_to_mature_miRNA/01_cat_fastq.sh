#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/small_rna_mapping/fastq/"
outputdir="/pine/scr/m/i/mikelaff/small_rna_mapping/fastq/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping_bowtie_to_mature_miRNA/RNAIDnumbers.txt"

source ~/.bashrc

mkdir -p ${outputdir}

for sample in `cat ${sampleList}`
do
	echo ${sample}
	
	# submit one job per cat
	#ssub --notify=OFF --wrap=\"cat ${datadir}${sample}*.fastq.gz \> ${outputdir}${sample}.fastq.gz\"
	# use login node or one job submission
	cat ${datadir}${sample}*.fastq.gz > ${outputdir}${sample}.fastq.gz
done

