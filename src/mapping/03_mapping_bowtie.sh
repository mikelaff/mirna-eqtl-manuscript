#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/small_rna_mapping/cutadapt/"
outputroot="/pine/scr/m/i/mikelaff/small_rna_mapping/bowtie/"

# bowtie index
bowtieIndex="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/bowtie_index/hg38"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping/RNAIDnumbers.txt"

source ~/.bashrc
module add bowtie/1.2.2
module add picard/2.10.3

for sample in `cat ${sampleList}`
do
	echo ${sample}
	mkdir -p ${outputroot}/${sample}
	cd ${outputroot}/${sample}

	# check for "Finished" in cutadapt log file. only do job if not present
	#if grep --quiet --no-messages 'Finished' ${sample}.log; then
	#	echo "FINISHED"
	#else
		# submit job	
		jid1=$(ssub --notify=OFF\
			    --mem=32g\
			    --wrap=\"bowtie ${bowtieIndex}\
			     			-n 2\
						-l 15\
						-k 50\
						--best\
						--strata\
						-S\
						--chunkmbs 512\
						${dataroot}/${sample}/${sample}.trimmed.fastq.gz\
						1\>${sample}.sam\
						2\>${sample}.log\")
		ssub --notify=OFF\
			-d afterok:$jid1\
			--mem=32g\
			--wrap=\"java -Xmx16g\
			-jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar\
			SortSam\
				INPUT=${sample}.sam\
				OUTPUT=${sample}.sortedByCoord.bam\
				SORT_ORDER=coordinate\
				CREATE_INDEX=TRUE\"
	#fi
done

