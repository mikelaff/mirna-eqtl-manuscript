#!/bin/bash

source ~/.bashrc

dataroot="/pine/scr/m/i/mikelaff/counts/"

outputdir_gene_counts="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/counts/small_rna_seq/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/counts/RNAIDnumbers.txt"

##Loop over each sample
for sample in `cat ${sampleList}`
do
	echo ${sample}
	
	# copy gene count files
	mv ${dataroot}/all_known_and_novel/${sample}/${sample}.all.known.and.novel.miRNA.counts ${outputdir_gene_counts}/all_known_and_novel/
	mv ${dataroot}/mirdeep_mirge/${sample}/${sample}.mirdeep.mirge.novel.miRNA.counts ${outputdir_gene_counts}/mirdeep_mirge/
	mv ${dataroot}/mirdeep_mirge_friedlander_nowakowski/${sample}/${sample}.mirdeep.mirge.friedlander.nowakowski.miRNA.counts ${outputdir_gene_counts}/mirdeep_mirge_friedlander_nowakowski/
done























