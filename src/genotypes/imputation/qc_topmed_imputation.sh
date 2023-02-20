#!/bin/bash

# imputation file directory
dataDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/raw/"

# output directory: plink1.9 hard-calls
hardDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/plink.hardcall.qc/"

# output directory: plink2.0 dosages
dosageDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/plink2.dosage.qc/"

# output directory: filtered vcf files
vcfDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/topmed_freeze5/vcf.qc/"

# fam file for sex information
famFile="/proj/steinlab/projects/FetalTissueQTL/genotypes/genotyped/AllSamplesQC.fam"

# QC filters
R2_THRESHOLD=0.3
#MISSING_GENOTYPE_RATE=0.05
HWE=0.000001
MAF=0.01

#module add samtools
module add plink/2.00a-20190527

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
#CHROMS=(X)
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	echo "Filtering ${chr} with bcftools"

	# bcftools to filter vcf files for R2 threshold and MAF
	#bcftools view -i "R2>${R2_THRESHOLD} & MAF>${MAF}" -Oz -o ${vcfDir}${chr}.topmed.dose.qc.vcf.gz ${dataDir}${chr}.dose.vcf.gz

	echo "Finished bcftools filtering"

	# plink1.9 for hard calls
#	module add plink/1.90b3
#	plink\
#		--vcf ${vcfDir}${chr}.topmed.dose.qc.vcf.gz\
#		--id-delim '_'\
#		--hwe ${HWE}\
#		--maf ${MAF}\
#		--update-sex ${famFile} 3\
#		--make-bed\
#		--out ${hardDir}${chr}.topmed.hardcall.qc

	# remove nosex file
#	rm ${hardDir}${chr}.topmed.hardcall.qc.nosex

	# plink2.0 for dosages
#	module add plink/2.00a-20190527
	plink2\
		--vcf ${vcfDir}${chr}.topmed.dose.qc.vcf.gz 'dosage=HDS'\
		--id-delim '_'\
		--hwe ${HWE}\
		--maf ${MAF}\
		--update-sex ${famFile} 'col-num=5'\
		--make-pgen\
		--out ${dosageDir}${chr}.topmed.dosage.qc

done

