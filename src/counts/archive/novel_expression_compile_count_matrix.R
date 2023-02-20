# compile novel expression counts from fetal tissue samples

library(here)
library(readr)
library(magrittr)
library(dplyr)

# OUTPUT FILES ####################################################################################
outputFile <- paste(here("results/known_and_novel_mirna/"), format(Sys.time(), "%Y%m%d"), "_novel_counts.tsv", sep="")

# INPUT FILES #####################################################################################
# directory of counts, one tsc per sample
directory <- here("results/known_and_novel_mirna/novel_counts2/")

# Compile Count Data ##############################################################################
files <- list.files(directory)

#load first gene count file
countData <- read_tsv(paste(directory, files[1], sep="/"),
                      skip = 2,
                      col_names = c("miRNA", "Chr", "Start", "End", "Strand", "Length", "Seq_Type", strsplit(files[1], "[.]")[[1]][1]))
countData %<>% dplyr::select(1,8)

#loop through all files in count folder and add gene count data to countData object
for (i in 2:length(files)) {
  #load count file
  tempData <- read_tsv(paste(directory, files[i], sep="/"),
                       skip = 2,
                       col_names = c("miRNA", "Chr", "Start", "End", "Strand", "Length", "Seq_Type", strsplit(files[i], "[.]")[[1]][1]))
  tempData %<>% dplyr::select(1,8)
  #join tempData with countData
  countData <- left_join(countData, tempData, by="miRNA")
}

#remove tempData object
rm(tempData)

# write count matrix
write_tsv(countData, outputFile)
