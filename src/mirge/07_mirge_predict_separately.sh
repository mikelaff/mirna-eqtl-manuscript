#!/bin/bash

# script moved to working directory (scratch space) before running
# working directory: /pine/scr/m/i/mikelaff/mirge/separately_predict/

#datadir="/pine/scr/m/i/mikelaff/mirge/fastq/"
#outputroot="/pine/scr/m/i/mikelaff/mirge/predict/"

#sampleList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirge/RNAIDnumbers.txt"
fileList="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/src/mirge/fastqFiles.txt"

mirgeLibs="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/data/mirge2.0/miRge.Libs"
bowtieBin="/nas/longleaf/apps/bowtie/1.1.2/bin"
samtoolsBin="/nas/longleaf/apps/samtools/1.8/bin"
rnafoldBin="/nas/longleaf/home/mikelaff/tools/ViennaRNA-2.3.5/install_dir/bin"

source ~/.bashrc
module add python/2.7.12
module add bowtie/1.1.2
module add samtools
export PYTHONPATH=/nas/longleaf/home/mikelaff/.local/lib/python2.7/site-packages:$PYTHONPATH

ssub --notify=ON --mem=500g -n 8 -N 1 --wrap=\"miRge2.0 predict\
					-s ${fileList}\
					-d miRBase\
					-pb ${bowtieBin}\
					-lib ${mirgeLibs}\
					-ps ${samtoolsBin}\
					-pr ${rnafoldBin}\
					-sp human\
					-ad illumina\
					-ai\
					-gff\
					-cpu 8\"

