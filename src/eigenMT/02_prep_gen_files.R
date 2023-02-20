# prep "GEN" files for eigenMT
# plink transposed, missing genotypes NA, 0 1 2
# one file per chrom

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES #########################################################################################################
# directory for gen files
output.dir <- paste0(here("results/eigenmt/"), results.name, "/gen_files/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTLor prefiltered genotypes, hardcalls
bfile.dir <- here("results/genotypes/mirQTLor_imputed_genotypes/plink.hardcall/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Export Files #########################################################################################################

for (chrom in CHROMS) {

    # use PLINK to transpose
    # prefiltered genotypes bfile
    bfile <- paste0(bfile.dir, chrom, ".topmed.hardcall.r2g03.maf001.mirQTLor.HomHetStrict")
    # output file prefix
    outputPrefix <- paste0(output.dir, chrom)
    # plink command
    command <- sprintf("plink --bfile %s --out %s --recode A-transpose",
                       bfile, outputPrefix)
    # call PLINK
    system(command = command)

    # import .traw file
    df.traw <- read_table2(paste0(output.dir, chrom, ".traw"))

    # select columns
    df.traw %<>%
        select(ID = SNP, starts_with("D"))

    # write tab sep file
    write_tsv(df.traw,
              path = paste0(output.dir, chrom, ".gen.txt"))

}

