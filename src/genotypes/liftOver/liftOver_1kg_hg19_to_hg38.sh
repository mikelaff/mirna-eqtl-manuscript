#!/bin/bash


#dataDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/1000G_phase3v5/raw/"
outputDir="/proj/steinlab/projects/FetalTissueQTL/genotypes/imputed/1000G_phase3v5/hg38/raw_hg38/"

chainFile="/nas/longleaf/home/mikelaff/tools/liftOver/hg19ToHg38.over.chain"
fastaReference="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/whole_genome_fasta/hg38.fa"

source ~/.bashrc
module add picard/2.18.22

CHROMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)
CHROMS=(22)
# loop over each vcf file
for chr in "${CHROMS[@]}"
do
	chr="chr${chr}"

	echo ${chr}

	ssub --notify=OFF\
		--mem=8g\
		-p steinlab\
		--wrap=\"java -Xmx16g -jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar\
		LiftoverVcf\
		I=${outputDir}tmp.${chr}.dose.mod.vcf.gz\
		O=${outputDir}${chr}.hg38.dose.vcf.gz\
		CHAIN=${chainFile}\
		REJECT=${outputDir}${chr}.rejected.vcf.gz\
		R=${fastaReference}\"
done



