#!/bin/bash

# input dir: plink1 filesets
input_dir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/1000genomes_phase3_hg38/EUR.plink/"

# output dir: plink1 filesets 
output_dir="/pine/scr/m/i/mikelaff/garfield/"

module add plink/2.00a-20190527

for chr in {1..22} X Y
do
	plink2 --bfile ${input_dir}EUR.chr${chr}_GRCh38\
		--set-all-var-ids '@:#'\
		--make-bed\
		--out ${output_dir}chr${chr}_updated_var_ids

done


