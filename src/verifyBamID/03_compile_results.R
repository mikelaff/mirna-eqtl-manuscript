# import small rna seq verify bam id best files

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaff)


# OUTPUT FILES ########################################################################################################
smallRNAseq.verifybamid.output.tsv <- here("results/verifyBamID/20190724_smallRNAseq_verifyBamID_bestSM.tsv")
totalRNAseq.verifybamid.output.tsv <- here("results/verifyBamID/20190724_totalRNAseq_verifyBamID_bestSM.tsv")

# INPUT FILES #########################################################################################################
smallRNAseq.verifybamid.data.dir <- here("results/verifyBamID/smallRNAseq_verifyBamID_bestSM_files/")
totalRNAseq.verifybamid.data.dir <- here("results/verifyBamID/totalRNAseq_verifyBamID_bestSM_files/")

# Import smallRNAseq Files ############################################################################################

# get file names
files <- list.files(smallRNAseq.verifybamid.data.dir)

# import first file
df.small <- read_tsv(file.path(smallRNAseq.verifybamid.data.dir, files[1]))
df.small %<>% rename(SEQ_ID = `#SEQ_ID`, NUM_SNPS = `#SNPS`, NUM_READS = `#READS`)
# rnaid for this file
rnaid <- unlist(strsplit(files[1], "\\."))[1]
# add rnaid to df
df.small$SEQ_ID[1] <- rnaid

for (i in 2:length(files)) {
    # import file
    tmp <- read_tsv(file.path(smallRNAseq.verifybamid.data.dir, files[i]))
    tmp %<>% rename(SEQ_ID = `#SEQ_ID`, NUM_SNPS = `#SNPS`, NUM_READS = `#READS`)
    # rnaid for this file
    rnaid <- unlist(strsplit(files[i], "\\."))[1]
    # add rnaid to df
    tmp$SEQ_ID[1] <- rnaid

    # append to bottom of df
    df.small <- bind_rows(df.small, tmp)
}

rm(tmp, i, files, rnaid)

# Import totalRNAseq Files ############################################################################################

# get file names
files <- list.files(totalRNAseq.verifybamid.data.dir)

# import first file
df.total <- read_tsv(file.path(totalRNAseq.verifybamid.data.dir, files[1]))
df.total %<>% rename(SEQ_ID = `#SEQ_ID`, NUM_SNPS = `#SNPS`, NUM_READS = `#READS`)
# rnaid for this file
rnaid <- unlist(strsplit(files[1], "\\."))[1]
# add rnaid to df
df.total$SEQ_ID[1] <- rnaid

for (i in 2:length(files)) {
    # import file
    tmp <- read_tsv(file.path(totalRNAseq.verifybamid.data.dir, files[i]))
    tmp %<>% rename(SEQ_ID = `#SEQ_ID`, NUM_SNPS = `#SNPS`, NUM_READS = `#READS`)
    # rnaid for this file
    rnaid <- unlist(strsplit(files[i], "\\."))[1]
    # add rnaid to df
    tmp$SEQ_ID[1] <- rnaid

    # append to bottom of df
    df.total <- bind_rows(df.total, tmp)
}

rm(tmp, i, files, rnaid)

# Write Output ########################################################################################################

write_tsv(df.small, smallRNAseq.verifybamid.output.tsv)
write_tsv(df.total, totalRNAseq.verifybamid.output.tsv)


# Scratch #############################################################################################################

df.small %<>%
    dplyr::select(rnaid = SEQ_ID, donor_dnaid = CHIP_ID, FREEMIX, CHIPMIX)

df.total %<>%
    dplyr::select(rnaid = SEQ_ID, donor_dnaid = CHIP_ID, FREEMIX, CHIPMIX)

df <- dplyr::full_join(df.small, df.total, by = "rnaid", suffix = c(".small", ".total"))

df$small <- !is.na(df$donor_dnaid.small)
df$total <- !is.na(df$donor_dnaid.total)

df$both <- df$small & df$total

tmp <- dplyr::filter(df, both)
tmp$match <- tmp$donor_dnaid.small == tmp$donor_dnaid.total
