#!/bin/bash

# EMMAX result files (.ps) need to be labeled for downstream analysis.
# label columns
# change permissions to prevent modifying these files

resultsDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/association_results/20201012_mirQTLor_4707_cond/"

# make dir for log files
mkdir -p ${resultsDir}logs/
# move logs
echo "Moving log files..."
mv ${resultsDir}*.log ${resultsDir}logs/

# make dir for .reml files
mkdir -p ${resultsDir}reml/
# move reml files
echo "Moving reml files..."
mv ${resultsDir}*.reml ${resultsDir}reml/

# loop over .ps files and modify header
echo "Modifying column headers..."
for file in `find ${resultsDir} -type f -name '*.ps'`
do
	sed -i '1 i\SNP\tBETA\tP' ${file}
done

# make dir for .ps files
mkdir -p ${resultsDir}ps/
# move ps files
echo "Moving ps files..."
mv ${resultsDir}*.ps ${resultsDir}ps/

# change permissions
echo "Changing permissions..."
chmod 440 ${resultsDir}ps/*.ps
