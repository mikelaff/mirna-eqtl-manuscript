#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# output dir sQTL ld
outputDir_sQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_sQTL/ld/"

# directory of plink filesets sQTL
# chromosome files labeled as chr*_F.uniq.bim
bfileDir_sQTL="/proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/"
# merged x chrom PAR and nonPAR into scratch folder for X chrom
bfileXDir_sQTL="/pine/scr/m/i/mikelaff/"

# list of snps to find ld sQTL
snpList_sQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_sQTL/index.snps.txt"

mkdir -p ${outputDir_sQTL}

for snp in `cat ${snpList_sQTL}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_sQTL}/${chr}_F.uniq"

	# only for X chrom
	if [[ "${chr}" == "chrX" ]]
	then
	    bfile="${bfileXDir_sQTL}/chrX_F.uniq"
	fi

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_sQTL}/sQTL.index.${snp}

done


