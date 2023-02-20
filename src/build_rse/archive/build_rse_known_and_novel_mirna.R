# build ranged summarized experiment with mirna count data for both miRBase and novel mirnas

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(SummarizedExperiment)
library(Biostrings)
library(RMySQL)

# OUTPUT FILES ########################################################################################################
output.rse.rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"),
                        "_rse_all_known_and_novel_counts.rds", sep="")

# INPUT FILES #########################################################################################################
# mirge quantified expression counts for miRBase miRNAs
mirbase.counts.tsv <- here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv")
# featureCounts quantified counts from bowtie mapped bams for novel miRNAs
novel.counts.tsv <- here("results/counts/small_rna_seq/20190501_mirdeep_mirge_friedlander_nowakowski_counts.tsv")
# rna sample summary file
sample.summary.tsv <- here("data/metadata/rna_sample_summary.tsv")
# sequencing run summary file
run.summary.tsv <- here("data/metadata/sequencing_run_summary_small_rna.tsv")
# flowcell data
flowcells.tsv <- here("data/metadata/flowcells.tsv")
# config.txt file for file names and lane numbers
config.txt <- here("data/lists/config.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# genotype and sex data
genotype.mds.tsv <- here("results/HM3mds/sample_mds_sex_by_dnaid_and_donorid.tsv")
# male rnaids that are not labeled by genotype data
male.unlabeled.rnaids.txt <- here("data/metadata/male_unlabeled_rnaids.txt")
# female rnaids that are not labeled by gneotype data
female.unlabeled.rnaids.txt <- here("data/metadata/female_unlabeled_rnaids.txt")

# mirge mirna merges
mirge.merges.csv <- here("data/mirge2.0/miRge.Libs/human/annotation.Libs/human_merges_miRBase.csv")

# GRanges for mirbase and novel mirna
# mirbase.granges.rds <- here("data/gtf_and_granges/20190430_mirbase_v22_mirna_granges.rds")
# fried.granges.rds <- here("data/gtf_and_granges/20190430_friedlander_mirna_granges.rds")
# nowa.granges.rds <- here("data/gtf_and_granges/20190430_nowakowski_mirna_granges.rds")
# mirdeep.granges.rds <- here("data/gtf_and_granges/20190430_mirdeep_novel_mirna_granges.rds")
# mirge.granges.rds <- here("data/gtf_and_granges/20190430_mirge_novel_mirna_granges.rds")

# Non-overlapping miRNA annotations for known and novel miRNAs
known.novel.granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# SteinLab Database ###################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")
rm(con, user_name, user_password, db_name, host_name)

# Import Non-overlapping miRNAs ########################################################################################
gr.mirna <- readRDS(known.novel.granges.rds)

# Import Counts and Format ############################################################################################
# read in count data frame
df.counts <- read_tsv(mirbase.counts.tsv, col_types = cols(.default = "i", miRNA = "c"))

# some rows have SNP countifications, not sure why, they don't have any counts, so i'm removing them
df.counts %<>% filter(!grepl("SNP", miRNA))
# miRge merged some miRNAs because their sequences are the same, this includes different SNP versions of a miRNA
# I want to expand the merged miRNAs so that both genomic loci attributed to a similar sequence will both have counts
# I don't want to expand SNP version, those should stay collapsed.
merges <- read_lines(mirge.merges.csv)
# loop over merges and expand to both miRNA names
for (mir in merges) {
    # look for slash, otherwise merge is of SNPs and dont want to expand
    if (grepl("/", mir)) {
        # the mir which will be expanded
        expand.mir <- unlist(strsplit(mir, ","))[1]
        # a list of mirs to expand into
        merge.list <- unlist(strsplit(mir, ","))[-1]
        # check expand mir is in list only once
        if (sum(grepl(expand.mir, df.counts$miRNA)) == 1) {
            # row to be expanded
            copyRow <- filter(df.counts, miRNA == expand.mir)
            # new df with expanded rows
            df.new <- data.frame(miRNA = merge.list, copyRow[2:241], stringsAsFactors = FALSE)
            # bind rows to df.counts
            df.counts <- bind_rows(list(df.counts, df.new))
            # remove old mir that was expanded
            df.counts %<>% filter(miRNA != expand.mir)

        } else {
            stop("error")
        }
    } else {
        next
    }
}
rm(df.new, copyRow, mir, merge.list, merges, expand.mir)

# romove SNP suffix, get rid of repeats
df.counts$miRNA <- sapply(strsplit(df.counts$miRNA, "\\."), `[`, 1)
df.counts %<>% distinct()

# create count matrix
mat.counts <- as.matrix(df.counts[,2:241])
rownames(mat.counts) <- df.counts$miRNA

# Import novel counts
df.counts.novel <- read_tsv(novel.counts.tsv)
# remove precursor counts
df.counts.novel %<>% filter(seq_type != "miRNA_putative_precursor")
# create count matrix
mat.counts.novel <- as.matrix(df.counts.novel[,11:250])
rownames(mat.counts.novel) <- df.counts.novel$name

all(colnames(mat.counts) == colnames(mat.counts.novel))

# bind count matrices
mat.counts <- rbind(mat.counts, mat.counts.novel)
rm(df.counts, df.counts.novel, mat.counts.novel)

# Load Runs Data ##################################################################################
# load config.txt file
runs <- read_tsv(config.txt, col_names = FALSE)
# label columns
runs <- dplyr::rename(runs, fastq = X1, seq_id = X2)
# remove 'fastq/' path prefix from fastq_file column
runs$fastq <- sapply(strsplit(runs$fastq, "/"), `[`, 2)
# remove .fastq
runs$fastq <- sapply(strsplit(runs$fastq, "\\."), `[`, 1)
# make new column for RNAID
runs$rnaid <- sapply(strsplit(runs$fastq, "_"), `[`, 1)
# make new column for sequencing lane
runs$sequencing_lane <- as.integer(gsub("[^[:digit:]]", "", sapply(strsplit(runs$fastq, "_"), `[`, 3)))

# Load Run Summary Metadata #######################################################################
# load sequencing run summary file
run.summary <- read_tsv(run.summary.tsv)

# merge run_summary into runs table
runs <- left_join(runs, run.summary, by = c("rnaid", "sequencing_lane"))
rm(run.summary)

# Load Sample Summary Metadata ####################################################################
# load rna sample file
sample.summary <- read_tsv(sample.summary.tsv)

# merge rna sample data with runs data by RNAID to create column data
# remove barcode_id from sample_summary to prevent duplcated columns while merging
col.data <- left_join(runs, dplyr::select(sample.summary, -barcode_id), by = "rnaid")
rm(sample.summary, runs)

# load flowcell data
flowcells <- read_tsv(flowcells.tsv)
flowcells <- dplyr::rename(flowcells, fastq_file = fastq)
flowcells$fastq <- sapply(strsplit(flowcells$fastq_file, "\\."), `[`, 1)
# merge into colData
col.data <- left_join(col.data, flowcells, by = "fastq")
rm(flowcells)

# load genotype data
genotypes <- read_tsv(genotype.mds.tsv)
col.data <- dplyr::select(col.data, -dnaid, -sex)
col.data <- left_join(col.data, genotypes, by = "donor_id")
rm(genotypes)

# load sex data
# read male and female rnaids
males <- read_lines(male.unlabeled.rnaids.txt)
females <- read_lines(female.unlabeled.rnaids.txt)
# label unlabeled males and females
col.data$sex[which(col.data$rnaid %in% males)] <- "Male"
col.data$sex[which(col.data$rnaid %in% females)] <- "Female"
rm(males, females)

# combine sequencing lanes across rnaids
for (i in 1:length(col.data$rnaid)) {
    col.data$sequencing_lane[i] <- paste(col.data$sequencing_lane[col.data$rnaid == col.data$rnaid[i]], collapse = "_")
}

# removed duplicated rows
col.data <- col.data[!duplicated(col.data$rnaid),]
col.data %<>% dplyr::select(-fastq)

# change some character variables to factor
col.data <- mutate(col.data, sequencing_lane = factor(sequencing_lane))
col.data <- mutate(col.data, sequencing_date = factor(sequencing_date))
col.data <- mutate(col.data, sequencing_project = factor(sequencing_project))
col.data <- mutate(col.data, sequencing_pool = factor(sequencing_pool))
col.data <- mutate(col.data, sequencing_run = factor(sequencing_run))

col.data <- mutate(col.data, barcode = factor(barcode))
col.data <- mutate(col.data, barcode_id = factor(barcode_id))

col.data <- mutate(col.data, rna_sample_index = factor(rna_sample_index))
col.data <- mutate(col.data, rna_plate_number = factor(rna_plate_number))
col.data <- mutate(col.data, rna_plate_well = factor(rna_plate_well))
col.data <- mutate(col.data, rna_lane = factor(rna_lane))

col.data <- mutate(col.data, rna_extraction_date = factor(rna_extraction_date))
col.data <- mutate(col.data, acquisition_date = factor(acquisition_date))

col.data <- mutate(col.data, tissue_section = factor(tissue_section, levels = c("CW", "CP", "GZ")))
col.data <- mutate(col.data, flowcell = factor(flowcell))
col.data <- mutate(col.data, unique_lane = factor(paste(col.data$flowcell, col.data$sequencing_lane, sep="_")))

col.data <- dplyr::select(col.data, -seq_id, -name, -ancestry, -not_usable, -source, -freezer_rack_row,
                         -freezer_rack_col, -freezer_box_row, -freezer_box_col, -freezer_cell_row,
                         -freezer_cell_col, -freezer, -freezer_comments, -exp_condition, -incubator,
                         -dna_extraction_date, -dna_purification_method, -dna_concentration, -dna_260_280,
                         -dna_260_230, -dna_source, -dna_extracted_by, -dna_freezer, -dna_freezer_box,
                         -genotype_id)

# format colData
# assign row names
col.data <- as.data.frame(col.data)
row.names(col.data) <- col.data$rnaid

rm(i)

# Load Row Data as GRanges ##############################################################

gr.mirbase <- readRDS(mirbase.granges.rds)
gr.fried <- readRDS(fried.granges.rds)
gr.nowa <- readRDS(nowa.granges.rds)
gr.mirdeep <- readRDS(mirdeep.granges.rds)
gr.mirge <- readRDS(mirge.granges.rds)

gr.all <- c(gr.mirbase, gr.fried, gr.nowa, gr.mirdeep, gr.mirge)
rm(gr.mirbase, gr.fried, gr.nowa, gr.mirdeep, gr.mirge)

#names(gr.all) <- gr.all$Name

sum(!rownames(mat.counts) %in% gr.all$Name)

sum(duplicated(gr.all$Name))
sum(duplicated(rownames(mat.counts)))

# rows of count matrix not in granges object
# most likely 5p and 3p version of mirbase mirnas that are only annotated as
# miR w/o 5p and 3p destinctions in mirbase
df.mirs <- data.frame(row_name = rownames(mat.counts),
                              stringsAsFactors = FALSE)

# name w/o 5p or 3p for secondary match
df.mirs$base <- paste("hsa-miR", sapply(strsplit(df.mirs$row_name, "-"), `[`, 3), sep="-")
# name of hairpin for tertiary match
df.mirs$hairpin <- paste("hsa-mir", sapply(strsplit(df.mirs$row_name, "-"), `[`, 3), sep="-")
# primary, secondary, or tertiary index into gr.all
df.mirs$gr_index <- match(df.mirs$row_name, gr.all$Name)
df.mirs$gr_index2 <- match(df.mirs$base, gr.all$Name)
df.mirs$gr_index3 <- match(df.mirs$hairpin, gr.all$Name)
# select index to use
df.mirs$gr_ind <- ifelse(!is.na(df.mirs$gr_index),
                         df.mirs$gr_index,
                         ifelse(!is.na(df.mirs$gr_index2),
                                df.mirs$gr_index2,
                                df.mirs$gr_index3))
# mirs to remove from count matrix (only 4)
mat.counts <- mat.counts[df.mirs$row_name[!is.na(df.mirs$gr_ind)],]
# remove last 4 NAs
df.mirs <- df.mirs[!is.na(df.mirs$gr_ind),]
# final rowranges for construction of ranged summarized experiment
row.ranges <- gr.all[df.mirs$gr_ind]

rm(gr.all, df.mirs)

# Build SummarizedExperiment ######################################################################

mode(mat.counts) <- "integer"
# format and check col.data
col.data <- col.data[colnames(mat.counts),]
all(rownames(col.data) == colnames(mat.counts))

# most rownames of mat.counts should equal the name in the row.ranges, but some are no longer unique
# set names and check row.ranges
sum(rownames(mat.counts) == row.ranges$Name)
names(row.ranges) <- rownames(mat.counts)
all(names(row.ranges) == rownames(mat.counts))

# build rse
rse <- SummarizedExperiment(assays = SimpleList(counts = mat.counts),
                            rowRanges = row.ranges,
                            colData = col.data)

# save rse
saveRDS(rse, output.rse.rds)


