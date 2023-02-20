#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/mirdeep/cutadapt/"
outputroot="/pine/scr/m/i/mikelaff/mirdeep/quant_known_and_novel/"

# bowtie index
#bowtieIndex="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/bowtie_index/hg38"

# precursor (hairpin) file
precursorFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/20190208_known_and_novel_hairpin.fa"

# mature mirna file
matureFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/known_and_novel_mirna/20190208_known_and_novel_mature.fa"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirdeep2/quant/RNAIDnumbers.txt"

source ~/.bashrc
module add bowtie/1.2.0
module add python/2.7.12

export PERL5LIB=$PERL5LIB:/nas/longleaf/home/mikelaff/tools/mirdeep2_0_0_8/lib
export PERL5LIB=$PERL5LIB:/nas/longleaf/home/mikelaff/perl5/lib/perl5

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
			    --wrap=\"gunzip -c ${dataroot}/${sample}/${sample}.trimmed.fastq.gz \>\
			    			${sample}.trimmed.fastq\")
		jid2=$(ssub --notify=OFF\
			    -d afterok:$jid1\
			    --mem=16g\
			    --wrap=\"mapper.pl ${sample}.trimmed.fastq\
			    			-e -h -j -l 18 -m\
						-s ${sample}.collapsed.fa\
						-v\")
		ssub --notify=OFF\
			-d afterok:$jid2\
			--mem=16g\
			--wrap=\"quantifier.pl -p ${precursorFile}\
					       -m ${matureFile}\
					       -r ${sample}.collapsed.fa\"
	#fi
done

