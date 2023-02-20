#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/smallrna/duplicate_marked_bams/"
outputdir="/pine/scr/m/i/mikelaff/smallrna/duplicate_marked_bams/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc
module add samtools

cd ${datadir}

for sample in `cat ${sampleList}`
do
	echo ${sample}
	
	if [ -f ${sample}.duplicateMarked.sortedByCoord.bai ]; then
		echo "Finished"
	else
		ssub --notify=OFF\
			--mem=16g\
			--wrap=\"samtools index\
				${datadir}/${sample}.duplicateMarked.sortedByCoord.bam\
				${datadir}/${sample}.duplicateMarked.sortedByCoord.bai\"
	fi
done

