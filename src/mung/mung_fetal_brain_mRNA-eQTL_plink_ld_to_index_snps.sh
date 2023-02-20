#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# output dir mQTL ld
outputDir_mQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/ld/"

# directory of plink filesets mQTL
# chromosome files labeled as chr*_F.uniq.bim
bfileDir_mQTL="/proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/"
# merged x chrom PAR and nonPAR into scratch folder for X chrom
bfileXDir_mQTL="/pine/scr/m/i/mikelaff/"

# list of snps to find ld mQTL
snpList_mQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/index.snps.txt"

mkdir -p ${outputDir_mQTL}

for snp in `cat ${snpList_mQTL}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_mQTL}/${chr}_F.uniq"

	# only for X chrom
	if [[ "${chr}" == "chrX" ]]
	then
	    bfile="${bfileXDir_mQTL}/chrX_F.uniq"
	fi

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mQTL}/mQTL.index.${snp}

done


