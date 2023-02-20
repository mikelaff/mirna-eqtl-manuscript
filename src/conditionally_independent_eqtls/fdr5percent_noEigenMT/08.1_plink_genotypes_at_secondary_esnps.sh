#!/bin/bash

# get genotype dosages at primary eSNPs

module add plink/2.00a-20190527

# name of association results folder
resultsName="20200120_mirQTLor_fdr5percent_noEigenMT"

# Plink imputed genotype file directory
pfileDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/genotypes/mirQTLor_imputed_genotypes/plink.dosage/"

# results directory
resultsDir="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/conditional_eqtls/${resultsName}/secondary/"

# samples file
keepFile="/proj/steinlab/projects/miRNA-eQTL-project/mirna-eqtl/results/emmax/samples/20200120_mirQTLor_DonorID_DNAID.txt"

# list of primary eSNPs
esnpsFile="${resultsDir}/${resultsName}_secondary_eSNPs.txt"

# output directory for sample genotypes
outputDir="${resultsDir}/secondary_eSNP_genotypes/"

mkdir -p ${outputDir}

while read esnp
do
    echo ${esnp}

    chr="$(echo ${esnp} | cut -d ":" -f 1)"

    pfile="${pfileDir}/${chr}.topmed.dosage.r2g03.maf001.mirQTLor.HomHetStrict"
    outFile="${outputDir}/${esnp}"

    plink2 --pfile ${pfile}\
	    --out ${outFile}\
	    --keep ${keepFile}\
	    --snp ${esnp}\
	    --export A

done < ${esnpsFile}







