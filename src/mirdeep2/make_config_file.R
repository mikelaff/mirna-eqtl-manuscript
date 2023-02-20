# script to create 3 character IDs for fastq file names used in miRDeep2
# updated for combined fastq files (one for each rnaid)
# 11 Dec 2018

library(tidyverse)
library(here)

# read in rnaids
df <- read_csv(here("src/mirdeep2/RNAIDnumbers.txt"), col_names = c("rnaid"))

# make file path and name
df$path <- paste("fastq/", df$rnaid, ".fastq", sep="")

# make seqid
df$seqid <- sprintf("%03d", seq(1,length(df$rnaid)))

# output fastq paths and seqids
write_tsv(select(df, path, seqid), here("src/mirdeep2/config.txt"), col_names = FALSE)

