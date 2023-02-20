#!/bin/bash

datadir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/small_rna_seq_bams/"
outputroot="/pine/scr/m/i/mikelaff/VarifyBamID/"

# list of RNAID numbers
sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/VarifyBamID/RNAIDnumbers.txt"

# VCF file with all QCed genotype data
vcfFile="/proj/steinlab/projects/R00/atac-qtl/GenotypeData/PlinkData_original/QCMergedData/allmgergedGenoMafHweQCnewVCFhg38.vcf"

# verifyBamID executable
verifyBamID="/proj/steinlab/projects/sharedApps/VerifyBamID/verifyBamID_1.1.3/verifyBamID/bin/verifyBamID"

source ~/.bashrc

for sample in `cat ${sampleList}`
do
	echo ${sample}
	mkdir -p ${outputroot}/${sample}
	cd ${outputroot}/${sample}

	# check for "Finished"
	if grep --quiet --no-messages 'finished' ${sample}.VerifyBamID.log; then
		echo "FINISHED"
	else
		# clear directory in case this is a redo job
		rm -f *
		# submit job	
		ssub --notify=OFF\
			--mem=32g\
			--wrap=\"${verifyBamID}\
				--vcf ${vcfFile}\
				--bam ${datadir}/${sample}.smallRNAseq.duplicateMarked.sortedByCoord.bam\
				--out ${sample}.VerifyBamID\
				--best\"
	fi
done

