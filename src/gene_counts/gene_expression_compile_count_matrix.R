# compile gene expression counts from fetal tissue samples
# write gene expression matrix

library(here)
library(readr)
library(magrittr)
library(dplyr)

# OUTPUT FILES ####################################################################################
outputFile <- paste(here("results/gene_counts/"), format(Sys.time(), "%Y%m%d"), "_fetalTissue_gene_counts.tsv", sep="")

# INPUT FILES #####################################################################################
# directory of gene counts, one csv per sample
directory <- here("results/gene_counts/counts/")

# Compile Count Data ##############################################################################
files <- list.files(directory)

#load first gene count file
countData <- read_tsv(paste(directory, files[1], sep="/"),
                      skip = 2,
                      col_names = c("ENSG", "Chr", "Start", "End", "Strand", "Length", strsplit(files[1], "[.]")[[1]][1]))
countData %<>% dplyr::select(1,7)

#loop through all files in count folder and add gene count data to countData object
for (i in 2:length(files)) {
  #load count file
  tempData <- read_tsv(paste(directory, files[i], sep="/"),
                       skip = 2,
                       col_names = c("ENSG", "Chr", "Start", "End", "Strand", "Length", strsplit(files[i], "[.]")[[1]][1]))
  tempData %<>% dplyr::select(1,7)
  #join tempData with countData
  countData <- left_join(countData, tempData, by="ENSG")
}

#remove tempData object
rm(tempData)

# write count matrix
write_tsv(countData, outputFile)
