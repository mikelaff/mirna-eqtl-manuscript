#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# directory root of plink hardcalls mirQTL
bfileDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# list of snps to find ld overlap index snps
snpList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_compiled/index.snps.txt"

# output dir mirQTL ld
outputDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_compiled/ld/"

mkdir -p ${outputDir_mirQTL}

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_mirQTL}/${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mirQTL}/conditional.mirQTLor.index.${snp}

done
