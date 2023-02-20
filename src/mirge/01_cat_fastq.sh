#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/mirge/fastq/"
outputdir="/pine/scr/m/i/mikelaff/mirge/fastq/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirge/RNAIDnumbers.txt"

source ~/.bashrc

mkdir -p ${outputdir}

for sample in `cat ${sampleList}`
do
	echo ${sample}
	
	ssub --notify=OFF --wrap=\"cat ${datadir}${sample}*.fastq.gz \> ${outputdir}${sample}.fastq.gz\"
done

