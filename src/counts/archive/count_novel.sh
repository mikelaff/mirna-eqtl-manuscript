#!/bin/bash

gtfFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/20190421_novel.gtf"

dataroot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"
outputroot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/novel_counts/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/count_novel/RNAIDnumbers.txt"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}

	outputdir="${outputroot}${sample}"
	mkdir -p ${outputdir}
	cd ${outputdir}
	
	if grep --quiet --no-messages 'Status' ${sample}.miRNA.counts.summary; then
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
				    --extraAttributes seq_type\
				    -a ${gtfFile}\
				    -o ${sample}.miRNA.counts\
				    ${dataroot}/${sample}.sortedByCoord.bam\"
	fi
done

