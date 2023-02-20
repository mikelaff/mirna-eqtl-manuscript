#!/bin/bash
#compile expression data (.csv) files from mirdeep2 quantification of known miRNA in mirbase

#fastq prefixes which are the directory names
prefixlist="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/lists/fastq_prefixes.txt"
#mirdeep2 output root
mirdeepdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirdeep2/mirbase_and_novel_expression/mirdeep2_runs_by_rnaid/"
#dir for compiled expression files
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirdeep2/mirbase_and_novel_expression/mirbase_and_novel_expression_by_rnaid/"

for dir in `cat ${prefixlist}`
do
  #move to directory
  cd ${mirdeepdir}${dir}
  #echo "current working dir: `pwd`"
  #copy expression file to output dir
  cp miRNAs_expressed_all_samples_15* "${outputdir}${dir}_mirbase_expression.csv"
done
