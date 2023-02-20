#!/bin/bash

# create freqx files for hardcall data in order to create list of snps which pass hom.minor and het sample filters

# output directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/plink.hardcall/"

module add plink/1.90b3

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(22)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	# plink1.9 for hard calls
	plink\
		--bfile ${hardDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTL\
		--freqx\
		--out ${hardDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTL
	
done

