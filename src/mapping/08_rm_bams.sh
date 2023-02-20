#!/bin/bash

dir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}
	rm ${dir}/${sample}.sortedByCoord.bam
	rm ${dir}${sample}.sortedByCoord.bai
done

