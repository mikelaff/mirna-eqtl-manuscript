#!/bin/bash

source ~/.bashrc
module add python/3.6.6

# path to eigenMT executable
eigen="~/tools/eigenMTwithTestData/eigenMT.py"

# qtl files dir
dir_qtl="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/20200120_mirQTLor/qtl_files/"

# gen files dir
dir_gen="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/20200120_mirQTLor/gen_files/"

# genpos files dir
dir_genpos="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/20200120_mirQTLor/genpos_files/"

# phepos files dir
dir_phepos="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/20200120_mirQTLor/phepos_files/"

# output directory
dir_output="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/eigenmt/20200120_mirQTLor/output/"

mkdir -p ${dir_output}

for chr in `seq 1 22` X
do

	ssub -p steinlab --notify=OFF --mem=16g --wrap=\"python ${eigen}\
								--CHROM chr${chr}\
								--QTL ${dir_qtl}chr${chr}.qtl.txt\
								--GEN ${dir_gen}chr${chr}.gen.txt\
								--GENPOS ${dir_genpos}chr${chr}.genpos.txt\
								--PHEPOS ${dir_phepos}chr${chr}.phepos.txt\
								--cis_dist 1000000\
								--OUT ${dir_output}chr${chr}\"

done
