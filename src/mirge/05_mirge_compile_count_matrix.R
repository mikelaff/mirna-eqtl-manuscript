# compile mirge2.0 miRNA counts csv files into a count matrix
# from both miRBase v22 and MirGeneDB v2

library(here)
library(readr)
library(dplyr)
library(magrittr)

# INPUT FILES #####################################################################################
# counts directory roots
root_dir_mirbase <- here("results/mirge2.0/counts/miRBase_v22/")
root_dir_mirgene <- here("results/mirge2.0/counts/MirGeneDB_v2/")
# OUTPUT FILES ####################################################################################
# output files
outputFile_mirbase <- paste(here("results/mirge2.0/"), format(Sys.time(), "%Y%m%d"), "_mirbase_v22_mirge2.0_counts.tsv", sep="")
outputFile_mirgene <- paste(here("results/mirge2.0/"), format(Sys.time(), "%Y%m%d"), "_mirgenedb_v2_mirge2.0_counts.tsv", sep="")

# Compile miRBase v22 Results #####################################################################
# files
files <- list.files(root_dir_mirbase)

# load first count file
cts <- read_csv(paste(root_dir_mirbase, files[1], sep=""))

# loop through each file
for(i in 2:length(files)){
  # read in counts csv
  suppressMessages(
    df <- read_csv(paste(root_dir_mirbase, files[i], sep=""))
    )
  
  # join to cts df
  cts <- full_join(cts, df, by = "miRNA")
}
rm(df, i)

# modify column names
colnames(cts) <- sapply(strsplit(colnames(cts), "\\."), `[`, 1)

# remove total counts row
cts <- filter(cts, miRNA != "miRNAtotal")

# replace na values
cts[is.na(cts)] <- 0

cts %<>% mutate_if(is.double, as.integer)

# write counts data table
write_tsv(cts, outputFile_mirbase, col_names = TRUE)

# Compile MirGeneDB v2 Results ####################################################################
# files
files <- list.files(root_dir_mirgene)

# load first count file
cts <- read_csv(paste(root_dir_mirgene, files[1], sep=""))

# loop through each file
for(i in 2:length(files)){
  # read in counts csv
  suppressMessages(
    df <- read_csv(paste(root_dir_mirgene, files[i], sep=""))
  )
  
  # join to cts df
  cts <- full_join(cts, df, by = "miRNA")
}
rm(df, i)

# modify column names
colnames(cts) <- sapply(strsplit(colnames(cts), "\\."), `[`, 1)

# remove total counts row
cts <- filter(cts, miRNA != "miRNAtotal")

# replace na values
cts[is.na(cts)] <- 0

cts %<>% mutate_if(is.double, as.integer)

# write counts data table
write_tsv(cts, outputFile_mirgene, col_names = TRUE)

