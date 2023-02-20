#!/bin/bash

module add plink/1.90b3

# compile genotype data across the various genotyping runs
# convert Final Reports to plink files: .lgen, .map, .fam
# make .bed .bim files with plink

# 10 genotype runs, 10 sets of data to process

# Set 2: 2012-272
setName="set2"

# Raw Data Directory (Final Reports)
dataDir="/proj/steinlab/raw/R00CellCultureQTL/Genotypes/genotypes/RawData/2012-272/FinalReports/"

# Output Directory
outputDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/compiled_raw_genotypes/set2/"
mkdir -p ${outputDir}

# plink file outputs
lgenFile="${outputDir}/${setName}.lgen"
mapFile="${outputDir}/${setName}.map"
famFile="${outputDir}/${setName}.fam"

# temp file
tmpFile="${outputDir}/tmp.txt"

# loop over files and create .lgen file
for file in `ls ${dataDir}2012-272_*.bz2`
do
	echo ${file}

	# unzip to temp
	bzip2 -cd ${file} > ${tmpFile}

	# parse text file to .lgen
	awk 'BEGIN {FS=" ";} {if (NR > 11) {if ($5 == "-") {print $1,$1,$2,0,0} else {print $1,$1,$2,$5,$6}}}' ${tmpFile} >> ${lgenFile}

	# remove temp file
	rm ${tmpFile}
done

# use first file to create .map file
bzip2 -cd `ls ${dataDir}2012-272_*.bz2 | awk '{if (NR==1) print $1}'` > ${tmpFile}
awk 'BEGIN {FS=" ";} {if (NR > 11) print $3,$2,0,$4}' ${tmpFile} > ${mapFile}
rm ${tmpFile}

# get unique IDs to create .fam file
awk '{print $1}' ${lgenFile} | uniq | awk '{print $1,$1,0,0,0,-9}' > ${famFile}

echo ${lgenFile}
echo ${mapFile}
echo ${famFile}

# use plink to make .bed .bim .fam combination
plink --lfile ${outputDir}${setName} --out ${outputDir}${setName} --make-bed --allow-no-sex
