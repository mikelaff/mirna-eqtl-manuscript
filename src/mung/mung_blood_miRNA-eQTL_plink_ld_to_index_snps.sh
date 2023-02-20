#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# output dir sQTL ld
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/huan_2015_blood_miRNA-eQTL/ld/"

# directory of plink filesets for 1000g eur data
bfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/EUR.plink/"

# list of snps to find ld
snpList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/huan_2015_blood_miRNA-eQTL/index.snps.txt"

mkdir -p ${outputDir}

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# rsid
	rsid=`basename ${snp} | cut -d ':' -f 2`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir}/EUR.${chr}_GRCh38.uniqueRSID"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${rsid}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir}/blood.index.${rsid}

done


