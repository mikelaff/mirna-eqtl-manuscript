#!/bin/bash

# loop over EMMAX result files (.ps)
# 2. clump based on nominal p-value threshold
# 3. make list of index snps per miR

module add plink/1.90b3

# name of association results folder
resultsName="20200120_mirQTLor"

# directory root of plink hardcalls (by chrom, by phenotype file)
bfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/"

# parents association results directory
resultsRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/association_results/${resultsName}/"

# ps file directory
psDir="${resultsRoot}/ps/"

# output directory for clumped results
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/${resultsName}/clumped/"

# nominal pvalue file
nomPvalFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/${resultsName}/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt"

PVALUE_THRESHOLD=`cat ${nomPvalFile}`

echo ${PVALUE_THRESHOLD}

#files=`find ${psFileDir} -type f -name '*.ps'`

mkdir -p ${outputDir}

for file in `find ${psDir} -type f -name '*.ps'`
do
	#echo ${file}

	# name for this mir
	uniName=`basename ${file} | cut -d '.' -f 5`

	#echo ${uniName}

	# chromosome number
	chr=`basename ${file} | cut -d '.' -f 1`

	#echo ${chr}

	# output prefix
	outPrefix=`basename ${file} .ps`

	echo ${outPrefix}

	# plink bfile for this ps file
	bfile="${bfileRoot}/${chr}/${chr}.hardcall.prefiltered.mirQTLor.${uniName}"

	#echo ${bfile}

	# clump using plink
	plink --bfile ${bfile}\
		--clump ${file}\
		--clump-kb 250\
		--clump-p1 ${PVALUE_THRESHOLD}\
		--clump-p2 0.01\
		--clump-r2 0.2\
		--out ${outputDir}${outPrefix}

	# list of index SNPs
	if [ -f ${outputDir}${outPrefix}.clumped ]
	then
		awk '{if (NR!=1) {print $3}}' ${outputDir}${outPrefix}.clumped | awk NF > ${outputDir}${outPrefix}.indexSNPs
	fi
done






