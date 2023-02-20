#!/bin/bash

dataroot="/proj/steinlab/raw/FetalTissueQTL/smallRNA/CorticalWall/FASTQ/"
outputdir="/pine/scr/m/i/mikelaff/mirge/fastq/"

source ~/.bashrc

mkdir -p ${outputdir}

for file in `find ${dataroot} -name *.fastq.gz`
do
	ssub --notify=OFF --wrap=\"cp ${file} ${outputdir}\"
done

