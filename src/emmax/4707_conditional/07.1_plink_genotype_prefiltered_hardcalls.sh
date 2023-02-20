#!/bin/bash

# EMMAX eQTL Analysis
# Mike Lafferty
# Stein Lab
# Nov 2019

# Prepare input hardcall genotype files using plink.
# tped/tfam files are needed for emmax

# Genotpyes files have already been filtered for mirQTLor samples, maf > 0.01, and hom/het sample thresholds.

source ~/.bashrc
module add plink/1.90b3

# Plink imputed genotype file directory
bfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# genes tsv file
genesFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/phenotype_files/20201012_mirQTLor_VST_miRNA_expression_residual_4707_edu_attain/20201012_mirQTLor_genes_hg38.txt"

# Output directory
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/"

# make output dir
mkdir -p ${outputDir}

# loop over each chromosome and make directory
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X XPAR)
for chr in "${chroms[@]}"
do
	mkdir -p ${outputDir}chr${chr}
done

while IFS=$'\t' read -r -a line
do

	uni="${line[0]}"
	chr="${line[1]}"
	str=$((${line[2]}-1000000))
	stp=$((${line[3]}+1000000))

	if [ $str -lt 0 ]
	then
		str=0
	fi

	bfile="${bfileDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"
	outFile="${outputDir}${chr}/${chr}.hardcall.prefiltered.mirQTLor.${uni}"

	# create bed/bim/fam files
	plink --bfile ${bfile}\
		--out ${outFile}\
		--chr ${chr:3}\
		--from-bp ${str}\
		--to-bp ${stp}\
		--make-bed

	# create frqx files
	plink --bfile ${outFile}\
		--out ${outFile}\
		--freqx

	# transpose for emmax
	plink --bfile ${outFile}\
		--out ${outFile}\
		--output-missing-genotype 0\
		--recode 12 transpose

	# transpose for .traw file
	plink --bfile ${outFile}\
		--out ${outFile}\
		--recode A-transpose

	#echo ${bfile}
	#echo ${outFile}

	#echo "${str}	${stp}"

done < ${genesFile}

