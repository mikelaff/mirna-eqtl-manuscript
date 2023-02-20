#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/mirge/fastq/"
outputroot="/pine/scr/m/i/mikelaff/mirge/predict/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirge/RNAIDnumbers.txt"

mirgeLibs="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirge2.0/miRge.Libs"
bowtieBin="/nas/longleaf/apps/bowtie/1.1.2/bin"
samtoolsBin="/nas/longleaf/apps/samtools/1.8/bin"
rnafoldBin="/nas/longleaf/home/mikelaff/tools/ViennaRNA-2.3.5/install_dir/bin"

source ~/.bashrc
module add python/2.7.12
module add bowtie/1.1.2
module add samtools
export PYTHONPATH=/nas/longleaf/home/mikelaff/.local/lib/python2.7/site-packages:$PYTHONPATH

for sample in `cat ${sampleList}`
do
	echo ${sample}

	outputdir="${outputroot}${sample}/"
	mkdir -p ${outputdir}
	cd ${outputdir}

	if [ ! -f miRge.*/annotation.report.csv ]; then
		echo "NOT FINISHED"
		rm -rf miRge.*

		ssub --notify=OFF --mem=32g -n 4 --wrap=\"miRge2.0 predict\
							-s ${datadir}${sample}.fastq.gz\
							-d miRBase\
							-pb ${bowtieBin}\
							-lib ${mirgeLibs}\
							-ps ${samtoolsBin}\
							-pr ${rnafoldBin}\
							-sp human\
							-ad illumina\
							-ai\
							-gff\
							-cpu 4\"
	else
		echo "finished"
	fi
done
