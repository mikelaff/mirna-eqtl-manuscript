#!/bin/bash

safFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping_bowtie_to_mature_miRNA/miRs.saf"

dataroot="/pine/scr/m/i/mikelaff/small_rna_mapping/bowtie/"
outputroot="/pine/scr/m/i/mikelaff/small_rna_mapping/featureCounts/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping_bowtie_to_mature_miRNA/RNAIDnumbers.txt"

outputSuffix=".mir4707.counts"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}

	outputdir="${outputroot}${sample}"
	mkdir -p ${outputdir}
	cd ${outputdir}
	
	if grep --quiet --no-messages 'Status' ${sample}${outputSuffix}.summary; then
		echo "FINISHED"
	else
		rm -f *
	
		ssub --notify=OFF -p steinlab\
	     	--wrap=\"featureCounts\
				    -s 0\
				    -M\
				    -O\
				    -a ${safFile}\
				    -F SAF\
				    -o ${sample}${outputSuffix}\
				    ${dataroot}/${sample}/${sample}_all.sortedByCoord.bam\"
	fi
done

