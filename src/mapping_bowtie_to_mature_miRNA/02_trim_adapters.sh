#!/bin/bash

datadir="/pine/scr/m/i/mikelaff/small_rna_mapping/fastq/"
outputroot="/pine/scr/m/i/mikelaff/small_rna_mapping/cutadapt/"

# small RNA seq Illumina 3' adapter
readAdapter="TGGAATTCTCGGGTGCCAAGG"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping_bowtie_to_mature_miRNA/RNAIDnumbers.txt"

source ~/.bashrc
module add cutadapt/1.18

for sample in `cat ${sampleList}`
do
	echo ${sample}
	mkdir -p ${outputroot}/${sample}
	cd ${outputroot}/${sample}

	# check for "Finished" in cutadapt log file. only do job if not present
	if grep --quiet --no-messages 'Finished' ${sample}.log; then
		echo "FINISHED"
	else
		# clear directory in case this is a redo job
		rm *
		# submit job	
		ssub --notify=OFF\
			     --wrap=\"cutadapt -a ${readAdapter}\
			     			--minimum-length 15\
				                -o ${sample}.trimmed.fastq.gz\
				   	    	${datadir}/${sample}.fastq.gz \>\
					    	${sample}.log\"
	fi
done

