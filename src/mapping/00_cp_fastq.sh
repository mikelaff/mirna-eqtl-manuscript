#!/bin/bash

dataroot="/proj/steinlab/raw/FetalTissueQTL/smallRNA/CorticalWall/FASTQ/"
outputdir="/pine/scr/m/i/mikelaff/small_rna_mapping/fastq/"

source ~/.bashrc

mkdir -p ${outputdir}

for file in `find ${dataroot} -name *.fastq.gz`
do
	echo ${file}

	# submit one job per copy
	#ssub --notify=OFF --wrap=\"cp ${file} ${outputdir}\"
	# use login node or submitting one job for coping
	cp ${file} ${outputdir}
done

