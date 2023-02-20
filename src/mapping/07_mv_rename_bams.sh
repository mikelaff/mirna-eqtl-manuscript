#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/smallrna/duplicate_marked_bams/"
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}
	mv ${dataroot}/${sample}.duplicateMarked.sortedByCoord.bam ${outputdir}/${sample}.smallRNAseq.duplicateMarked.sortedByCoord.bam
	mv ${dataroot}/${sample}.duplicateMarked.sortedByCoord.bai ${outputdir}/${sample}.smallRNAseq.duplicateMarked.sortedByCoord.bai
done

