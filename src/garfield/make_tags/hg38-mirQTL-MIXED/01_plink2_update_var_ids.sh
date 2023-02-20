#!/bin/bash

# input dir: plink1 filesets
input_dir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/"

# output dir: plink1 filesets 
output_dir="/pine/scr/m/i/mikelaff/garfield-mirQTL-MIXED/"

module add plink/2.00a-20190527

for chr in {1..22} X
do
	plink2 --bfile ${input_dir}chr${chr}.topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict\
		--set-all-var-ids '@:#'\
		--make-bed\
		--out ${output_dir}chr${chr}_updated_var_ids

done


