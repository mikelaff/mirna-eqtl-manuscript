#!/bin/bash

# working dir: plink1 filesets 
working_dir="/pine/scr/m/i/mikelaff/garfield-mirQTL-MIXED/"

module add plink/1.90b3

for chr in {1..22} X
do
	plink --bfile ${working_dir}chr${chr}_updated_var_ids\
		--maf 0.0001\
		--make-bed\
		--out ${working_dir}chr${chr}_updated_var_ids_maf

done


