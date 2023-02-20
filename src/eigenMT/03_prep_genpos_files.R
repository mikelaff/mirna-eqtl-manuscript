# prep "GENPOS" files for eigenMT
# 'snp', 'chr_snp', 'pos'
# one file per chrom

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES #########################################################################################################
# directory for genpos files
output.dir <- paste0(here("results/eigenmt/"), results.name, "/genpos_files/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# .traw files created in 02_prep_gen_files.R
traw.dir <- paste0(here("results/eigenmt/"), results.name, "/gen_files/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Export Files #########################################################################################################

for (chrom in CHROMS) {

    print(paste("Importing chr:", chrom))

    # import .traw file
    df.traw <- read_table2(paste0(traw.dir, chrom, ".traw"))

    # select columns
    df.traw %<>%
        select(snp = SNP,
               pos = POS) %>%
        mutate(chr_snp = chrom) %>%
        select(snp, chr_snp, pos)

    # write tab sep file
    write_tsv(df.traw,
              path = paste0(output.dir, chrom, ".genpos.txt"))

}

