#!/bin/bash

# EMMAX eQTL Analysis
# Mike Lafferty
# Stein Lab

# use EMMAX for eQTL analysis

source ~/.bashrc
module add emmax/07

#emmax -v -d 10 -Z -t [tped_prefix] -p [pheno_file] -k [kin_file] -c [cov_file] -o [out_prefix]

# Plink tped/tfam files directory
tfileDirRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_dosages/"

# Kinship Matrix directory
kfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/kinship_mat/mirQTLor/"

# Covariates file: autosomes
#cfileAuto="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/covariates/20191204_chrAutosomes.mirQTLor.genotypePCA_10.expressionPCA_10.pool.purMethod.rin.sex.gestWeek.cov"

# Covariates file: chromosome X (non-PAR)
#cfileChrX="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/covariates/20191204_chrX.mirQTLor.genotypePCA_10.expressionPCA_10.pool.purMethod.rin.gestWeek.cov"

# Phenotype files directory root
phenoDirRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/secondary/20200120_mirQTLor_fdr5percent_noEigenMT_VST_miRNA_expression_residual_after_secondary_eqtl/"

# Output directory
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/tertiary/association_results/"

# create output directory
mkdir -p ${outputDir}

# loop over each chromosome
chroms=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X XPAR)
for chr in "${chroms[@]}"
do
	echo "Submitting jobs for chr${chr}..."

	# kinship mat for this chrom (no XPAR kinship file, so use chrX)
	if [[ "${chr}" == "XPAR" ]]
	then
		kfile="${kfileDir}chrX.mirQTLor.BN.kinf"
		echo "Using kinship matrix: ${kfile}"
	else
		kfile="${kfileDir}chr${chr}.mirQTLor.BN.kinf"
		echo "Using kinship matrix: ${kfile}"
	fi

	# phenotype directory for this chrom
	phenoDir="${phenoDirRoot}chr${chr}/"

	# tfile directory for this chrom
	tfileDir="${tfileDirRoot}chr${chr}/"

	# loop over all phenotype files in phenotype directory
	for phenoFile in `ls ${phenoDir}`
	do
		# unique name for this mir
		uniName="$(echo ${phenoFile} | cut -d'.' -f3)"

		#echo ${uniName}

		#echo ${phenoFile}

		tfile="chr${chr}.dosage.prefiltered.mirQTLor.${uniName}"

		#echo ${tfile}

		outputPrefix="chr${chr}.dosage.prefiltered.mirQTLor.tertiary.${uniName}"

		#echo ${outputPrefix}

		# ssub -p general --notify=OFF --mem=4g --wrap=\"emmax -v -d 10 -Z\
		# 						-t ${tfileDir}${tfile}\
		# 						-k ${kfile}\
		# 						-o ${outputDir}${outputPrefix}\
		# 						-p ${phenoDir}${phenoFile}\"

		emmax -v -d 10 -Z\
				-t ${tfileDir}${tfile}\
				-k ${kfile}\
				-o ${outputDir}${outputPrefix}\
				-p ${phenoDir}${phenoFile}
	done
done

