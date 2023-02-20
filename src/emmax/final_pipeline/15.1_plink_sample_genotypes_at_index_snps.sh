#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# name of association results folder
resultsName="20200120_mirQTLor"

# directory root of plink hardcalls (by chrom, by phenotype file)
hardcallBfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/"

# directory root of plink dosages
dosageBfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_dosages/"

# parents association results directory
resultsRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/association_results/${resultsName}/"

# clumped file directory
clumpDir="${resultsRoot}/clumped/"

# output directory for sample genotypes
outputDir="${resultsRoot}/sample_genotypes/"

#files=`find ${psFileDir} -type f -name '*.ps'`

mkdir -p ${outputDir}

for file in `find ${clumpDir} -type f -name '*.indexSNPs'`
do
	#echo ${file}

	# name for this mir
	uniName=`basename ${file} | cut -d '.' -f 5`

	#echo ${uniName}

	# chromosome number
	chr=`basename ${file} | cut -d '.' -f 1`

	#echo ${chr}

	# plink bfile for this ps file
	bfile="${hardcallBfileRoot}/${chr}/${chr}.hardcall.prefiltered.mirQTLor.${uniName}"

	#echo ${bfile}

	# loop over index snps for this mir
	for snp in `cat ${file}`
	do
		echo ${snp}

		plink --bfile ${bfile}\
			--recode AD\
			--snp ${snp}\
			--out ${outputDir}${chr}.hardcall.prefiltered.mirQTLor.${snp}
	done

done

module add plink/2.00a-20190527

for file in `find ${clumpDir} -type f -name '*.indexSNPs'`
do
	#echo ${file}

	# name for this mir
	uniName=`basename ${file} | cut -d '.' -f 5`

	#echo ${uniName}

	# chromosome number
	chr=`basename ${file} | cut -d '.' -f 1`

	#echo ${chr}

	# plink pfile for this ps file
	pfile="${dosageBfileRoot}/${chr}/${chr}.dosage.prefiltered.mirQTLor.${uniName}"

	#echo ${bfile}

	# loop over index snps for this mir
	for snp in `cat ${file}`
	do
		echo ${snp}

		plink2 --pfile ${pfile}\
			--export AD\
			--snp ${snp}\
			--out ${outputDir}${chr}.dosage.prefiltered.mirQTLor.${snp}
	done

done






