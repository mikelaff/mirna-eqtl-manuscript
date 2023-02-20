#!/bin/bash

datadir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirge2.0/predict_miRBase_v22_individually/"
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/mirge2.0/novel/miRBase_v22/"

source ~/.bashrc

cd ${datadir}

for rnaid in `ls ${datadir}`
do
  echo ${rnaid}
  cd ${datadir}${rnaid}
  file=`find ${datadir}${rnaid} -name "${rnaid}_novel_miRNAs_report.csv"`
  echo ${file}
  cp ${file} ${outputdir}
  #base=`basename ${file} _R1_001.fq`
  #sample=`echo ${base} | cut -d '_' -f 1`
done

