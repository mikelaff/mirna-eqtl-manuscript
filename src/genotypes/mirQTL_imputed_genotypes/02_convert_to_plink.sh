#!/bin/bash

# filtered vcf file directory
dataDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/vcf_r2g03/"

# output directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/plink.hardcall/"

# output directory: plink2.0 dosages
dosageDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/plink.dosage/"

# fam file for sex information
famFile="/proj/steinlab/projects/FetalTissueQTL/genotypes/genotyped/AllSamplesQC.fam"

# QC filters
#R2_THRESHOLD=0.3
#MISSING_GENOTYPE_RATE=0.05
HWE=0.000001
#MAF=0.01

#module add samtools
#module add plink/2.00a-20190527

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(22)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	# plink1.9 for hard calls
	module add plink/1.90b3
	plink\
		--vcf ${dataDir}${chr}.topmed.dose.qc.vcf.gz\
		--id-delim '_'\
		--update-sex ${famFile} 3\
		--make-bed\
		--out ${hardDir}${chr}.topmed.hardcall.r2g03

	# remove nosex file
	rm ${hardDir}${chr}.topmed.hardcall.r2g03.nosex

	# plink2.0 for dosages
	module add plink/2.00a-20190527
	plink2\
		--vcf ${dataDir}${chr}.topmed.dose.qc.vcf.gz 'dosage=HDS'\
		--id-delim '_'\
		--update-sex ${famFile} 'col-num=5'\
		--make-pgen\
		--out ${dosageDir}${chr}.topmed.dosage.r2g03

done

