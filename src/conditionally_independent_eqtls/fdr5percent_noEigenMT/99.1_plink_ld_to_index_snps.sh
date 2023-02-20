#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# directory root of plink hardcalls (by chrom, by phenotype file)
bfileRoot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# list of index snps
snpList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/index.snps.txt"

# output directory for LD
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/ld/"

#files=`find ${psFileDir} -type f -name '*.ps'`

mkdir -p ${outputDir}

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileRoot}/${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir}/mirQTL.index.${snp}

done
