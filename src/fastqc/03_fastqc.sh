#!/bin/bash
#Mike Lafferty
#Jason Stein Lab
#University of North Carolina at Chapel Hill

dataroot_fq="/pine/scr/m/i/mikelaff/small_rna_mapping/fastq/"
dataroot_trimmed_fq="/pine/scr/m/i/mikelaff/small_rna_mapping/cutadapt/"

outputroot_fq="/pine/scr/m/i/mikelaff/small_rna_mapping/fastQC/fq/"
outputroot_trimmed_fq="/pine/scr/m/i/mikelaff/small_rna_mapping/fastQC/trimmed_fq/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc
module add fastqc

for sample in `cat ${sampleList}`
do
	echo ${sample}

	#directory for output for fq
	outputdir_fq="${outputroot_fq}/${sample}"

	mkdir -p ${outputdir_fq}
	cd ${outputdir_fq}
	
	if [ ! -f ${sample}_fastqc.html ]; then
		ssub --notify=OFF --wrap=\"fastqc --outdir ${outputdir_fq} ${dataroot_fq}/${sample}.fastq.gz\"
	fi

	# skip for not finished cutadapt jobs
	#if [ "${sample}" = "R00RL0023" ] || [ "${sample}" = "R00RL0024" ] || [ "${sample}" = "R00RL0070" ] || [ "${sample}" = "R00RL0071" ]; then
		#echo "skipping"
		#continue
	
	#directory for output for trimmed fq
	outputdir_trimmed_fq="${outputroot_trimmed_fq}/${sample}"

	mkdir -p ${outputdir_trimmed_fq}
	cd ${outputdir_trimmed_fq}

	if [ ! -f ${sample}.trimmed_fastqc.html ]; then
		ssub --notify=OFF --wrap=\"fastqc --outdir ${outputdir_trimmed_fq} ${dataroot_trimmed_fq}/${sample}/${sample}.trimmed.fastq.gz\"
	fi
done

