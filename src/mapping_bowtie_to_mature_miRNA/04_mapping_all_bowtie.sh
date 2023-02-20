#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/small_rna_mapping/cutadapt/"
outputroot="/pine/scr/m/i/mikelaff/small_rna_mapping/bowtie/"

# bowtie index
bowtieIndex="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/bowtie_index/miRge_human_mirna_SNP_pseudo_miRBase/miRge"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mapping_bowtie_to_mature_miRNA/RNAIDnumbers.txt"

source ~/.bashrc
module add bowtie/1.2.3
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
			    --mem=8g\
			    --wrap=\"bowtie ${bowtieIndex}\
			     			-n 0\
						--norc\
						--all\
						-S\
						${dataroot}/${sample}/${sample}.trimmed.fastq.gz\
						1\>${sample}_all.sam\
						2\>${sample}_all.log\")
		ssub --notify=OFF\
			-d afterok:$jid1\
			--mem=8g\
			--wrap=\"java -Xmx16g\
			-jar /nas/longleaf/apps/picard/2.10.3/picard-2.10.3/picard.jar\
			SortSam\
				INPUT=${sample}_all.sam\
				OUTPUT=${sample}_all.sortedByCoord.bam\
				SORT_ORDER=coordinate\
				CREATE_INDEX=TRUE\"
	#fi
done

