#!/bin/bash

module add plink/1.90b3

# list of index snps
snpList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_compiled/index.snps.txt"

# directory root of plink hardcalls
hardcallBfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# directory root of plink dosages
dosageBfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.dosage/"

# parents results directory
resultsRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_compiled/"

# output directory for sample genotypes
outputDir="${resultsRoot}/sample_genotypes/"

mkdir -p ${outputDir}

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# plink bfile
	bfile="${hardcallBfileRoot}/${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"

	plink --bfile ${bfile}\
			--recode AD\
			--snp ${snp}\
			--out ${outputDir}${chr}.hardcall.prefiltered.mirQTLor.${snp}

done

module add plink/2.00a-20190527

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# plink pfile
	pfile="${dosageBfileRoot}/${chr}.topmed.dosage.r2g03.maf001.mirQTLor.HomHetStrict"

    plink2 --pfile ${pfile}\
			--export AD\
			--snp ${snp}\
			--out ${outputDir}${chr}.dosage.prefiltered.mirQTLor.${snp}

done








