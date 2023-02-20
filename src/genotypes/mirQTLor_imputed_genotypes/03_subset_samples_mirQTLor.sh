#!/bin/bash

# input directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/plink.hardcall/"

# input directory: plink2.0 dosages
dosageDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTL_imputed_genotypes/plink.dosage/"

# output directory: plink1.9 hard-calls
outHardDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# output directory: plink2.0 dosages
outDosageDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.dosage/"

# samples file
keepFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/mirQTLor_samples.txt"

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
		--bfile ${hardDir}${chr}.topmed.hardcall.r2g03\
		--keep ${keepFile}\
		--hwe ${HWE}\
		--make-bed\
		--out ${outHardDir}${chr}.topmed.hardcall.r2g03.mirQTLor

	# plink2.0 for dosages
	module add plink/2.00a-20190527
	plink2\
		--pfile ${dosageDir}${chr}.topmed.dosage.r2g03\
		--keep ${keepFile}\
		--hwe ${HWE}\
		--make-pgen\
		--out ${outDosageDir}${chr}.topmed.dosage.r2g03.mirQTLor

done

