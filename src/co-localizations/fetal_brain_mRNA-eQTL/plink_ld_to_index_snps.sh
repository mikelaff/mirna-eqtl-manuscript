#!/bin/bash

# find LD to index snps

module add plink/1.90b3

# directory root of plink hardcalls mirQTL
bfileDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# directory of plink filesets mQTL
# chromosome files labeled as chr*_F.uniq.bim
bfileDir_mQTL="/proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/"

# list of snps to find ld overlap index snps
snpList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/co-localization/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_overlap.index.snps.txt"

# output dir mQTL ld
outputDir_mQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/co-localization/fetal_brain_mRNA-eQTL/ld/mQTL/"
# output dir mirQTL ld
outputDir_mirQTL="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/co-localization/fetal_brain_mRNA-eQTL/ld/mirQTL/"

mkdir -p ${outputDir_mQTL}
mkdir -p ${outputDir_mirQTL}

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_mQTL}/${chr}_F.uniq"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mQTL}/mQTL.overlap.index.${snp}

done

for snp in `cat ${snpList}`
do
	# chromosome
	chr=`basename ${snp} | cut -d ':' -f 1`

	# chromosome number
	#chrNum=${chr#"chr"}

	# plink bfile
	bfile="${bfileDir_mirQTL}/${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict"

	plink --bfile ${bfile}\
		--r2\
		--ld-snp ${snp}\
		--ld-window 4000\
		--ld-window-kb 1000\
		--ld-window-r2 0.2\
		--out ${outputDir_mirQTL}/mirQTL.overlap.index.${snp}

done
