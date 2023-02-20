# build count matrix and summarized experiment for mirdeep2 quant data from mirbase v22
# needs summarized experiment from mirge data for phenotype data

library(readr)
library(dplyr)
library(magrittr)
library(here)
library(reshape2)


# OUTPUT ####################
# count matrix tsv
outputFile_tsv <- paste(here("results/mirdeep2/"), format(Sys.time(), "%Y%m%d"), "_mirbase_v22_mirdeep2_counts.tsv", sep="")
# summarized experiment with count data
output_se <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_mirbase_v22_mirdeep2_counts_se.rds", sep="")

# INPUT #####################
# mirge se for phenotype data
mirge_se <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# counts csv directory
csv_dir <- here("results/mirdeep2/mirbase_v22/counts_mirbaseV22/")

# Compile mirdeep2 miRBase v22 Results #############################################
# files
files <- list.files(csv_dir)

# load first count file
df <- read_tsv(paste(csv_dir, files[1], sep=""))

# mirnas annotated by mirna and precursor pair, to create a unique name for each mirna,
# cat together mirna and precursor. reads appear to be mapping to more than one mirna,
# precursor pair. also there are two mirna/precursor pairs that have duplicated,
# (not sure why): hsa-miR-3142_hsa-mir-3142, hsa-miR-4487_hsa-mir-4487
# just combining rows for these two
df %>%
  select(mirna = `#miRNA`, precursor, read_count) %>%
  mutate(mirna_precursor = paste(mirna, precursor, sep="_")) %>%
  select(mirna_precursor, read_count) %>%
  group_by(mirna_precursor) %>%
  summarise_all(sum) %>%
  rename(!!strsplit(files[1], "[.]")[[1]][1] := read_count) -> cts

# loop through each file
for(i in 2:length(files)){
  read_tsv(paste(csv_dir, files[i], sep="")) %>%
    select(mirna = `#miRNA`, precursor, read_count) %>%
    mutate(mirna_precursor = paste(mirna, precursor, sep="_")) %>%
    select(mirna_precursor, read_count) %>%
    group_by(mirna_precursor) %>%
    summarise_all(sum) %>%
    rename(!!strsplit(files[i], "[.]")[[1]][1] := read_count) -> df
  
  # join to cts df
  cts <- full_join(cts, df, by = "mirna_precursor")
}
rm(df, i)

# output count matrix to file
cts[2:241] <- lapply(cts[2:241], as.integer)
write_tsv(cts, outputFile_tsv)


# Create summarized experiment ##################
library(DESeq2)
se.mirge <- readRDS(mirge_se)

coldata <- as.data.frame(colData(se.mirge))
rm(se.mirge)

m <- as.matrix(cts[,-1])
rownames(m) <- cts$mirna_precursor
cts <- m
rm(m)

coldata <- coldata[colnames(cts),]
all(rownames(coldata) == colnames(cts))

se <- SummarizedExperiment(assays=SimpleList(counts=cts), colData=coldata)

saveRDS(se, output_se)
