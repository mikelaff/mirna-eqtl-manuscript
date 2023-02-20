#!/bin/bash

# EMMAX eQTL Analysis
# Mike Lafferty
# Stein Lab

# Prepare input dosage genotype files using plink.
# tped/tfam files are needed for emmax

# Genotpyes files have already been filtered for mirQTLor samples, MAF > 0.01, and Hom. Minor Samples != 1 & Het Samples > 1.

source ~/.bashrc
module add plink/2.00a-20190527

# Plink imputed genotype file directory
pfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.dosage/"

# genes tsv file
genesFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/phenotype_files/20201012_mirQTLor_VST_miRNA_expression_residual_4707_edu_attain/20201012_mirQTLor_genes_hg38.txt"

# tfam file
tfamFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/mirQTLor.tfam"

# Output directory
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_dosages/"

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

	pfile="${pfileDir}${chr}.topmed.dosage.r2g03.maf001.mirQTLor.HomHetStrict"
	outFile="${outputDir}${chr}/${chr}.dosage.prefiltered.mirQTLor.${uni}"

	# create pgen/psam/pvar files
	plink2 --pfile ${pfile}\
		--out ${outFile}\
		--chr ${chr:3}\
		--from-bp ${str}\
		--to-bp ${stp}\
		--make-pgen

	# create .acount files
	plink2 --pfile ${outFile}\
		--out ${outFile}\
		--freq counts

	# transpose for emmax
	plink2 --pfile ${outFile}\
		--out ${outFile}\
		--export A-transpose

	echo "Creating .tped and .tfam files..."

	# check for successful .traw or skip
	if [ -f ${outFile}.traw ]; then

		# convert to tped
		cut -f-4,7- ${outFile}.traw | tail -n +2 > ${outFile}.tped
		# copy tfam from hardcalls
		cp ${tfamFile} ${outFile}.tfam
	else
		echo "No .traw file found"
	fi

	#echo ${bfile}
	#echo ${outFile}

	#echo "${str}	${stp}"

done < ${genesFile}

