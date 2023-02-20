#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/small_rna_mapping/bowtie/"
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}
	mv ${dataroot}/${sample}/${sample}.sortedByCoord.* ${outputdir}
done

