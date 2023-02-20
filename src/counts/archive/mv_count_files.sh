#!/bin/bash

source ~/.bashrc

dataroot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/novel_counts/"

outputdir_gene_counts="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/novel_counts2/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/count_novel/RNAIDnumbers.txt"

##Loop over each sample
for sample in `cat ${sampleList}`
do
	echo ${sample}
	
	# copy gene count files
	mv ${dataroot}/${sample}/${sample}.miRNA.counts ${outputdir_gene_counts}
done























