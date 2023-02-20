#!/bin/bash

# find LD to colocalized snps

module add plink/1.90b3

# output dir mQTL ld
outputDir_mQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/ld/mQTL/"
# output dir mirQTL ld
outputDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/ld/mirQTL/"

# directory of plink filesets mQTL
# "/proj/steinlab/projects/R00/eQTLanalysis/genofiles/Genotypehg38/22.dose.R2g03.QC"
bfileDir_mQTL="/proj/steinlab/projects/R00/eQTLanalysis/genofiles/Genotypehg38/"
# directory of plink filesets mirQTL
# "/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/chr8.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"
bfileDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# sample keep file for mQTL
keepFile_mQTL="/proj/steinlab/projects/R00/FetalBraineQTLs/new.donorlist.txt"

# list of snps to find ld mQTL
snpList_mQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalized_snps_mQTL.txt"
# list of snps to find ld mirQTL
snpList_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalized_snps_mirQTL.txt"

mkdir -p ${outputDir_mQTL}
mkdir -p ${outputDir_mirQTL}

for snp in `cat ${snpList_mQTL}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_mQTL}/${chrNum}.dose.R2g03.QC"

	plink --bfile ${bfile}\
		--keep ${keepFile_mQTL}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 2000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mQTL}/mQTL.${snp}

done

for snp in `cat ${snpList_mirQTL}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# plink bfile
	bfile="${bfileDir_mirQTL}/${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 2000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mirQTL}/mirQTL.${snp}

done
