#!/bin/bash

# working dir: plink1 filesets 
working_dir="/pine/scr/m/i/mikelaff/garfield-mirQTL-MIXED/"

source ~/.bashrc
module add plink/1.90b3

for chr in {1..22} X
do

	ssub -p general --mem=16g --notify=OFF --wrap=\"plink\
		--bfile ${working_dir}chr${chr}_updated_var_ids_maf\
		--show-tags all\
		--tag-kb 1000\
		--tag-r2 0.8\
		--out ${working_dir}chr${chr}_08\"

done


