#!/bin/bash

# EMMAX eQTL Analysis
# Mike Lafferty
# Stein Lab
# Nov 2019

# use EMMAX to make a series of kinship matrices, one for each chrom (excluding that chrom in creation)

# /proj/steinlab/projects/R00/atac-qtl/EMMAXResult/KinshipMatrix/MakeKinshipMatrix.sh
# by Dan Liang used as a template.

source ~/.bashrc
module add emmax/07
module add plink/1.90b3

# plink bfile for non-imputed genotypes, all chromosomes
bfile="/proj/steinlab/projects/FetalTissueQTL/genotypes/genotyped/AllSamplesQC"

# plink tfile and kinship mat output dir
tfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/kinship_mat/mirQTLor/"

# files to keep in the miRNA-eQTL analysis
keepFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/samples/20200120_mirQTLor_DonorID_DNAID.txt"

mkdir -p ${tfileDir}

chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)

for chr in "${chroms[@]}"
do

	# run plink to exclude chr and create tfiles
	jid1=$(ssub -p general --notify=OFF --mem=8g --wrap=\"plink --bfile ${bfile}\
							--out ${tfileDir}chr${chr}.mirQTLor\
							--keep ${keepFile}\
							--not-chr ${chr}, 23-26\
							--output-missing-genotype 0\
							--recode 12 transpose\")

	# run emmax for kinship mat creation
	ssub -p general --notify=OFF --mem=8g -d afterok:${jid1} --wrap=\"emmax-kin -v -d 10 ${tfileDir}chr${chr}.mirQTLor\"

done

# run plink for chrX and create tfiles
jid1=$(ssub -p general --notify=OFF --mem=8g --wrap=\"plink --bfile ${bfile}\
						--out ${tfileDir}chrX.mirQTLor\
						--keep ${keepFile}\
						--not-chr 23-26\
						--output-missing-genotype 0\
						--recode 12 transpose\")

# run emmax for kinship mat creation
ssub -p general --notify=OFF --mem=8g -d afterok:${jid1} --wrap=\"emmax-kin -v -d 10 ${tfileDir}chrX.mirQTLor\"

