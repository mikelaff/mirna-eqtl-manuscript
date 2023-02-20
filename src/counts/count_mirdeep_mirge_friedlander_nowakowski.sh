#!/bin/bash

gtfFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/gtf_and_granges/20190430_mirdeep_mirge_friedlander_nowakowski_mirna.gtf"

dataroot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"
outputroot="/pine/scr/m/i/mikelaff/counts/mirdeep_mirge_friedlander_nowakowski/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/counts/RNAIDnumbers.txt"

outputSuffix=".mirdeep.mirge.friedlander.nowakowski.miRNA.counts"

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
	
		ssub --notify=OFF\
	     	--wrap=\"featureCounts\
				    -s 0\
				    -M\
				    -f\
				    -O\
				    -t miRNA\
				    -g Name\
				    --extraAttributes ID,Alias,Derives_from,seq_type\
				    -a ${gtfFile}\
				    -o ${sample}${outputSuffix}\
				    ${dataroot}/${sample}.sortedByCoord.bam\"
	fi
done

