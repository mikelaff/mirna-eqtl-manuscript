# build summarized experiment with mirna count data for both miRBase and novel mirnas
# cannot be ranged because mirge merged mirnas are produced in more than one location in the genome

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(SummarizedExperiment)
library(Biostrings)
library(RMySQL)

# OUTPUT FILES #########################################################################################################
output.se.rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"),
                        "_se_mirge_combined_and_non_overlapping_novel_mirna_counts.rds", sep="")

# INPUT FILES ##########################################################################################################
# mirge quantified expression counts for miRBase miRNAs
mirbase.counts.tsv <- here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv")
# featureCounts quantified counts from bowtie mapped bams for novel miRNAs
novel.counts.tsv <- here("results/counts/small_rna_seq/20190501_mirdeep_mirge_friedlander_nowakowski_counts.tsv")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# mirge mirna merges
mirge.merges.csv <- here("data/mirge2.0/miRge.Libs/human/annotation.Libs/human_merges_miRBase.csv")

# Non-overlapping miRNA annotations for known and novel miRNAs
known.novel.granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Import Non-overlapping miRNA #########################################################################################
gr.mirna <- readRDS(known.novel.granges.rds)

# Import Counts and Format #############################################################################################
# read in count data frame
df.counts <- read_tsv(mirbase.counts.tsv, col_types = cols(.default = "i", miRNA = "c"))

# create count matrix
mat.counts <- as.matrix(df.counts[,2:241])
rownames(mat.counts) <- df.counts$miRNA

# Import novel counts
df.counts.novel <- read_tsv(novel.counts.tsv)
# remove precursor counts
df.counts.novel %<>%
    filter(seq_type != "miRNA_putative_precursor")

# keep only novel counts in the non-overlapping list
df.counts.novel %<>%
    filter(name %in% gr.mirna$Name)

# create count matrix
mat.counts.novel <- as.matrix(df.counts.novel[,11:250])
rownames(mat.counts.novel) <- df.counts.novel$name

stopifnot(all(colnames(mat.counts) == colnames(mat.counts.novel)))

# bind count matrices
mat.counts <- rbind(mat.counts, mat.counts.novel)
rm(df.counts, df.counts.novel, mat.counts.novel)

# SteinLab Database ####################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")
rm(con, user_name, user_password, db_name, host_name)

# Load Metadata ########################################################################################################
df.metadata <- read_tsv(samples.metadata.tsv, col_types = cols(Lanes = "c"))

# combine with donor table to get gestation week and tissue acquisition date
df.metadata %<>%
    left_join(dplyr::select(donors, DonorID, TissueAcquisitionDate = AcquisitionDate, GestationWeek), by = "DonorID")

# rearrange columns
df.metadata %<>%
    dplyr::select(RNAID, DonorID, GestationWeek, TissueSection, TissueAcquisitionDate, everything())

rm(donors)

# format colData
# assign row names
col.data <- as.data.frame(df.metadata)
row.names(col.data) <- col.data$RNAID

# Import miRge Merged miRNAs ###########################################################################################
# Import mirge merged miRNAs to be added to row data
mirge.merges <- read_lines(mirge.merges.csv)

df.merges <- data.frame(miRNA = sapply(strsplit(mirge.merges, ","), `[`, 1),
                        merge.list = sapply(strsplit(mirge.merges, ","), function(x) paste(x[-1], collapse = ",")),
                        stringsAsFactors = FALSE)

stopifnot(all(df.merges$mirna %in% rownames(mat.counts)))

# format rowData
row.data <- data.frame(miRNA = rownames(mat.counts),
                       stringsAsFactors = FALSE)
row.data %<>%
    left_join(df.merges, by = "miRNA")
rownames(row.data) <- row.data$miRNA

# Build SummarizedExperiment ######################################################################

mode(mat.counts) <- "integer"

# format and check col.data
col.data <- col.data[colnames(mat.counts),]
stopifnot(all(rownames(col.data) == colnames(mat.counts)))


stopifnot(all(rownames(row.data) == rownames(mat.counts)))

# build se
se <- SummarizedExperiment(assays = SimpleList(counts = mat.counts),
                           rowData = row.data,
                           colData = col.data)

# save se
saveRDS(se, output.se.rds)


