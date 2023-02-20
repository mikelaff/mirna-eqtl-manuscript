#!/bin/bash

dataroot="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/verifyBamID/smallRNAseq_samples/"
outputdir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/verifyBamID/bestSM_files/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/verifyBamID/RNAIDnumbers.txt"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}

	cp ${dataroot}/${sample}/${sample}.VerifyBamID.bestSM ${outputdir}
done

