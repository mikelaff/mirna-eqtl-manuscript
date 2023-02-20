#!/bin/bash

vcfDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/vcf/"
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/"

source ~/.bashrc
module add plink/1.90b3

for chr in {1..22} X Y
do

	ssub -p general --mem=16g --notify=OFF --wrap=\"plink\
		--vcf ${vcfDir}ALL.chr${chr}_GRCh38.genotypes.20170504.vcf.gz\
		--vcf-half-call missing\
		--make-bed\
		--out ${outputDir}ALL.chr${chr}_GRCh38\"

done


