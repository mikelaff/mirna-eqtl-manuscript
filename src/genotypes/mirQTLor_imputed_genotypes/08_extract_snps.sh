#!/bin/bash

# extract just snps that pass hom. minor != 1 and het > 1

# output directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# output directory: plink2.0 dosages
dosageDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.dosage/"

module add plink/1.90b3
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(22)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	# plink1.9 for hard calls
	plink\
		--bfile ${hardDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTLor\
		--extract ${hardDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict.snps\
		--make-bed\
		--out ${hardDir}${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict

done

module add plink/2.00a-20190527
CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(22)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	# plink2.0 for dosages
	plink2\
		--pfile ${dosageDir}${chr}.topmed.dosage.r2g03.maf001.mirQTLor\
		--extract ${dosageDir}${chr}.topmed.dosage.r2g03.maf001.mirQTLor.HomHetStrict.snps\
		--make-pgen\
		--out ${dosageDir}${chr}.topmed.dosage.r2g03.maf001.mirQTLor.HomHetStrict
done

