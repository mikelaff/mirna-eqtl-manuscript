#!/bin/bash

inputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/ALL.plink/"
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/EUR.plink/"

# 1000Genomes unrelated EUR individuals
keepFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/1000G.EUR.QC.1.fam"

mkdir -p ${outputDir}

source ~/.bashrc
module add plink/1.90b3

for chr in {1..22} X Y
do

	ssub -p general --mem=8g --notify=OFF --wrap=\"plink\
		--bfile ${inputDir}ALL.chr${chr}_GRCh38\
		--keep ${keepFile}\
		--make-bed\
		--out ${outputDir}EUR.chr${chr}_GRCh38\"

done


