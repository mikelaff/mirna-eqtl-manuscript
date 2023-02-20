#!/bin/bash

# imputation file directory
dataDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/raw/"

# output directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/plink.hardcall.qc/"

# output directory: plink2.0 dosages
dosageDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/plink2.dosage.qc/"

# output directory: filtered vcf files
vcfDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/vcf.qc/"

# fam file for sex information
famFile="/proj/steinlab/projects/FetalTissueQTL/genotypes/genotyped/AllSamplesQC.fam"

# QC filters

source ~/.bashrc
module add samtools

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHROMS=(X)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	echo $chr

	# bcftools to filter vcf files for R2 threshold and MAF
	ssub --notify=ON\
		--mem=8g\
		--wrap=\"bcftools view -i \'R2\>0.3 \& MAF\>0.01\' -Oz -o ${vcfDir}${chr}.topmed.dose.qc.vcf.gz ${dataDir}${chr}.dose.vcf.gz\"


done

