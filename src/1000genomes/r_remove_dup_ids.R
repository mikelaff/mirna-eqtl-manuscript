# find and remove duplicate rsids in 1000 genomes EUR files

library(here)
library(dplyr)
library(readr)

# OUTPUT FILES #########################################################################################################
dir.1kg.eur.plink <- here("data/1000genomes_phase3_hg38/EUR.plink/")

# INPUT FILES ##########################################################################################################

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

# Find Duplicates ######################################################################################################

# loop over each chr, find duplicates, write to file
# for (chr in CHROMS) {
#    print(chr)
#
#    # 1000Genomes EUR .bim file
#    bimFile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.bim")
#
#    # duplicate rsid file
#    dupFile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.duplicatedRSIDs")
#
#    # load .bim plink file
#    df.bim <- read_tsv(bimFile, col_names = c("chr", "rsid", "cm", "pos", "a1", "a2"))
#
#    # duplicate rsids
#    dup.rsids <- unique(df.bim$rsid[duplicated(df.bim$rsid)])
#
#    # write duplicate rsids to file
#    write_lines(dup.rsids, dupFile)
#
# }

# Exclude Duplicates ###################################################################################################

# loop over each chr, remove duplicates using PLINK

for (chr in CHROMS) {

    # 1000Genomes EUR bfile
    bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")

    # duplicate rsid file
    dupFile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.duplicatedRSIDs")

    # output file prefix
    outputPrefix <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.uniqueRSID")

    # plink command
    command <- sprintf("plink --bfile %s --exclude %s --make-bed --out %s",
                       bfile, dupFile, outputPrefix)
    # call PLINK
    system(command = command)

}

