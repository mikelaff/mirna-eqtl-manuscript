#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# name of association results folder
resultsName="20200120_mirQTLor"

# directory root of plink hardcalls (by chrom, by phenotype file)
bfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/"

# parents association results directory
resultsRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/association_results/${resultsName}/"

# clumped file directory
clumpDir="${resultsRoot}/clumped/"

# output directory for LD
outputDir="${resultsRoot}/ld/"

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
	bfile="${bfileRoot}/${chr}/${chr}.hardcall.prefiltered.mirQTLor.${uniName}"

	#echo ${bfile}

	# loop over index snps for this mir
	for snp in `cat ${file}`
	do
		echo ${snp}

		plink --bfile ${bfile}\
			--r2\
			--ld-snp ${snp}\
			--ld-window 2000\
			--ld-window-kb 1000\
			--ld-window-r2 0.2\
			--out ${outputDir}${chr}.hardcall.prefiltered.mirQTLor.${snp}
	done

done






