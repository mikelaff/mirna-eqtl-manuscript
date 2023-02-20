#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/mirdeep/quant_mirbaseV22/"
outputdir="/pine/scr/m/i/mikelaff/mirdeep/counts_mirbaseV22/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirdeep2/quant/RNAIDnumbers.txt"

mkdir -p ${outputdir}
cd ${outputdir}

for sample in `cat ${sampleList}`
do
	echo ${sample}
	cp ${dataroot}/${sample}/miRNAs_expressed_all_samples_*.csv ${sample}.mirdeep2.mirbaseV22.counts.csv
done

