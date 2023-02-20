#!/bin/bash

datadir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirge2.0/annotate_miRBase_v22/"
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirge2.0/counts/miRBase_v22/"

source ~/.bashrc

cd ${datadir}

for rnaid in `ls ${datadir}`
do
  echo ${rnaid}
  cd ${datadir}${rnaid}
  file=`find ${datadir}${rnaid} -name 'miR.Counts.csv'`
  echo ${file}
  cp ${file} ${outputdir}${rnaid}.miR.counts.csv
  #base=`basename ${file} _R1_001.fq`
  #sample=`echo ${base} | cut -d '_' -f 1`
done

