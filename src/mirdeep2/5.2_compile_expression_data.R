# Analyze miRBase known miRNA expression in 480 smallRNA-seq samples
# from miRDeep2 expression analysis

# This script compiles the count data
# from the mirdeep2 quantification of expression .csv files.
# The data is compiled and saved as tsv files

library(readr)
library(dplyr)
library(magrittr)
library(here)
library(reshape2)

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
files %<>% rename(fastq_file = X1, seq_id = X2)
# remove 'fastq/' path prefix from fastq_file column
files$fastq_file <- sapply(strsplit(files$fastq_file, "/"), `[`, 2)
# make new column for RNAID
files$rnaid <- sapply(strsplit(files$fastq_file, "_"), `[`, 1)
# make new column for expression .csv file name
files$exp_data_file <- paste(sapply(strsplit(files$fastq_file, "\\."), `[`, 1), "_mirbase_expression.csv", sep="")

##################################################
# import data and compile into data frames
##################################################
# current filename to work with
exp_file <- paste(expression_data_dir, files$exp_data_file[1], sep="")

# current expression data imported
exp_data_reads <- read_tsv(exp_file, col_names = TRUE)
exp_data_reads %<>% rename(miRNA = `#miRNA`, seq_norm = `seq(norm)`)

# compute miRNA-precursor key
exp_data_reads %<>% mutate(key = paste(miRNA, precursor, sep = "_"))

# select only key and count
exp_data_reads %<>% select(key, count = seq)

# label each row with seq_id
exp_data_reads %<>% mutate(seq_id = files$seq_id[1])

# mark duplicate keys (from last to get initial key also marked as duplicated)
exp_data_reads %<>% mutate(dup_key = duplicated(key) | duplicated(key, fromLast = TRUE))

# # move duplicate keys into list
# #keys <- exp_data_tmp$key
# dups <- exp_data_tmp$key[duplicated(exp_data_tmp$key)]
# # remove data from duplicate keys
# #exp_data_tmp <- filter(exp_data_tmp, !(key %in% dups) )
# 
# # move expression data into compiled data frame
# exp_data_reads <- select(exp_data_tmp, miRNA, precursor, seq, key)
# #exp_data_rpm <- select(exp_data_tmp, miRNA, precursor, seq_norm, key)
# 
# # rename seq to seq_id
# names(exp_data_reads)[names(exp_data_reads) == "seq"] <- files$seq_id[1]
# #names(exp_data_rpm)[names(exp_data_rpm) == "seq_norm"] <- files$seq_id[1]


# loop through each expression files
#for(i in 2:length(files$exp_data_file)) {
for(i in 2:length(files$exp_data_file)) {
  # current filename to work with
  exp_file <- paste(expression_data_dir, files$exp_data_file[i], sep="")
  
  # current expression data imported
  exp_data_tmp <- read_tsv(exp_file, col_names = TRUE)
  exp_data_tmp %<>% rename(miRNA = `#miRNA`, seq_norm = `seq(norm)`)
  
  # compute miRNA-precursor key
  exp_data_tmp %<>% mutate(key = paste(miRNA, precursor, sep = "_"))
  
  # select only key and count
  exp_data_tmp %<>% select(key, count = seq)
  
  # label each row with seq_id
  exp_data_tmp %<>% mutate(seq_id = files$seq_id[i])
  
  # mark duplicate keys
  exp_data_tmp %<>% mutate(dup_key = duplicated(key) | duplicated(key, fromLast = TRUE))
  
  # bind rows into reads df
  exp_data_reads <- bind_rows(exp_data_reads, exp_data_tmp)
}
rm(exp_data_tmp)

# factor seq_id
exp_data_reads$seq_id <- factor(exp_data_reads$seq_id)
# integer count
exp_data_reads$count <- as.integer(exp_data_reads$count)

# reshape data to split counts across seq_id, aggrigate duplicated keys at floor(mean)
exp_data_reads <- dcast(exp_data_reads, key + dup_key ~ seq_id,
                        value.var = "count",
                        fun.aggregate = function(x) {as.integer(floor(mean(x)))})

# write files
write_tsv(exp_data_reads, paste(data_output_dir, "mirbase_and_novel_expression_reads.tsv", sep=""))
#write_tsv(exp_data_rpm, paste(data_output_dir, "mirbase_expression_rpm.tsv", sep=""))
