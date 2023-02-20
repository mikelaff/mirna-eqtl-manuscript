# Analyze miRBase known miRNA expression in 480 smallRNA-seq samples
# from miRDeep2 expression analysis
# This scrip looks for duplicates in the mibase mirna expression counts created by mirdeep2
# There are 2 miRNA/precuror combinations that are duplicated, most have count 0 in each sample
# 2 samples have counts = 1, this script confirmed that
# in script 5.2_compile_expression_data.R these duplicate lines are removed so they dont 
# interfere with any downstream analyses

library(readr)
library(dplyr)
library(here)

##################################################
# set working directory to mirna-eqtl project root
print(paste("Working directory: ", getwd(), sep=""))

##################################################
# set global variables
##################################################
# relative path to expression analysis data directory
expression_data_dir <- here("results/mirdeep2/mirbase_and_novel_expression/mirbase_and_novel_expression_by_rnaid/")
# relative path to output directory
data_output_dir <- here("results/mirdeep2/mirbase_and_novel_expression/")
# config.txt file used by mirdeep to label fastq files
config_file <- here("data/lists/config.txt")
# file holding fastq prefixes
fastq_prefixes_file <- here("data/lists/fastq_prefixes.txt")

##################################################
# load and configure files data frame
##################################################
# load config.txt file
files <- read_tsv(config_file, col_names = FALSE)
# label columns
files <- rename(files, fastq_file = X1, seq_id = X2)
# remove 'fastq/' path prefix from fastq_file column
files$fastq_file <- sapply(strsplit(files$fastq_file, "/"), `[`, 2)
# make new column for RNAID
files$rnaid <- sapply(strsplit(files$fastq_file, "_"), `[`, 1)
# make new column for expression .csv file name
files$exp_data_file <- paste(sapply(strsplit(files$fastq_file, "\\."), `[`, 1), "_mirbase_expression.csv", sep="")

##################################################
# import data and compile into data frames
##################################################
dups <- list()
dups_df <- tibble()

for(i in 1:length(files$exp_data_file)) {
  # current filename to work with
  exp_file <- paste(expression_data_dir, files$exp_data_file[i], sep="")
  # current expression data imported
  exp_data_tmp <- read_tsv(exp_file, col_names = TRUE)
  exp_data_tmp <- rename(exp_data_tmp, miRNA = `#miRNA`, seq_norm = `seq(norm)`)
  # compute miRNA-precursor key
  exp_data_tmp$key <- paste(exp_data_tmp$miRNA, exp_data_tmp$precursor, sep="")
  
  keys <- exp_data_tmp$key
  dups <- keys[duplicated(keys)]
  
  dups_df <- bind_rows(dups_df, filter(exp_data_tmp, key %in% dups))
}












