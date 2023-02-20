#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/mirge/fastq/"
outputroot="/pine/scr/m/i/mikelaff/mirge/annotate_known_and_novel/"

sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirge/RNAIDnumbers_1.txt"

mirgeLibs="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirge2.0/miRge.Libs"
bowtieBin="/nas/longleaf/apps/bowtie/1.1.2/bin"

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

		ssub --notify=ON --mem=32g -n 4 --wrap=\"miRge2.0 annotate\
							-s ${datadir}${sample}.fastq.gz\
							-d miRBase\
							-pb ${bowtieBin}\
							-lib ${mirgeLibs}\
							-sp novel\
							-ad illumina\
							-cpu 4\"
	else
		echo "finished"
	fi
done
