#!/bin/bash

# update chr codes to include "chr"

dataDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/1000G_phase3v5/raw/"
outputDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/1000G_phase3v5/hg38/raw_hg38/"

source ~/.bashrc

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHROMS=(22)
# loop over each vcf file
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	echo ${chr}

	ssub -p steinlab --mem=8g --wrap=\"gzip -dc ${dataDir}${chr}.dose.vcf.gz \| awk \'{if\(\\\$0 \!\~ /\^\#/\) print \"chr\"\\\$0\; else print \\\$0}\' \| gzip \> ${outputDir}tmp.${chr}.dose.mod.vcf.gz\"

done



