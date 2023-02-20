#!/bin/bash

# EMMAX eQTL Analysis
# Mike Lafferty
# Stein Lab
# July 2019

# use Plink to calculate MDS and PCA components
# 212 unique and verified donors with small-rna-seq data, exclude sex and MT chromosomes
# without HapMap3 samples

source ~/.bashrc
module add plink/1.90b3

# Plink genotype file prefix
bfile="/proj/steinlab/projects/FetalTissueQTL/genotypes/genotyped/AllSamplesQC"

# Plink keep file
keepFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/samples/20200120_mirQTLor_DonorID_DNAID.txt"

# Output directory
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/pca_mds/mirQTLor/"

# Output prefix
prefix="chrAutosomes.mirQTLor"

mkdir -p ${outputDir}

ssub --notify=ON --mem=8g --wrap=\"plink --bfile ${bfile}\
					--out ${outputDir}${prefix}\
					--keep ${keepFile}\
					--not-chr 23-26\
					--cluster\
					--pca 20\
					--mds-plot 20\"

