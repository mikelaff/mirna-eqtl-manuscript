#!/bin/bash

source ~/.bashrc
module add picard

datadir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"
outputdir="/pine/scr/m/i/mikelaff/smallrna/duplicate_marked_bams/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

mkdir -p ${outputdir}
cd ${outputdir}

##Loop over each sample
for sample in `cat ${sampleList}`
do
	echo ${sample}

	if [ -f ${sample}_duplication_metrics.txt ]; then
		echo "FINISHED"
	else
		ssub --notify=OFF\
			--mem=32g\
			--wrap=\"java -Xmx16g -jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar MarkDuplicates\
			I=${datadir}/${sample}.sortedByCoord.bam\
			O=${outputdir}/${sample}.duplicateMarked.sortedByCoord.bam\
			M=${outputdir}/${sample}_duplication_metrics.txt\"
	fi
done























