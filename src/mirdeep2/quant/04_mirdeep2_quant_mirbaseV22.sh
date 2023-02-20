#!/bin/bash

dataroot="/pine/scr/m/i/mikelaff/mirdeep/cutadapt/"
outputroot="/pine/scr/m/i/mikelaff/mirdeep/quant_mirbaseV22/"

# bowtie index
#bowtieIndex="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/genome/hg38_UCSC/bowtie_index/hg38"

# precursor (hairpin) file
precursorFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirbase22/hsa_hairpin.fa"

# mature mirna file
matureFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirbase22/hsa_mature.fa"

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
		ssub --notify=OFF\
			--mem=16g\
			--wrap=\"quantifier.pl -p ${precursorFile}\
					       -m ${matureFile}\
					       -r ../../quant_known_and_novel/${sample}/${sample}.collapsed.fa\"
	#fi
done

