#!/bin/bash

# imputation file directory
dataDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/raw/"

# output directory: filtered vcf files
vcfDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/vcf_r2g03/"

# QC filters

source ~/.bashrc
module add samtools

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(X)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	echo $chr

	# bcftools to filter vcf files for R2 threshold and MAF
	ssub --notify=ON\
		--mem=8g\
		--wrap=\"bcftools view -i \'R2\>0.3\' -Oz -o ${vcfDir}${chr}.topmed.dose.qc.vcf.gz ${dataDir}${chr}.dose.vcf.gz\"


done

