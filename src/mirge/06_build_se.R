# build summarized experiment with mirna count data for both miRBase v22 and MirGeneDB v2

library(here)
library(dplyr)
library(readr)
library(SummarizedExperiment)
library(Biostrings)

# OUTPUT FILES ####################################################################################
mirbase_output_file <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_mirbase_v22_mirge2.0_counts.rds", sep="")
mirgene_output_file <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_mirgenedb_v2_mirge2.0_counts.rds", sep="")

# INPUT FILES #####################################################################################
# mirge quantified expression counts
mirbase_counts_file <- here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv")
mirgene_counts_file <- here("results/mirge2.0/20180905_mirgenedb_v2_mirge2.0_counts.tsv")
# rna sample summary file
sample_summary_file <- here("data/metadata/rna_sample_summary.tsv")
# sequencing run summary file
run_summary_file <- here("data/metadata/sequencing_run_summary.tsv")
# flowcell data
flowcell_file <- here("data/metadata/flowcells.tsv")
# miRBase homo sapiens chrom coordiantes for miRNAs (GRCh38)
hsa_gff_file <- here("data/mirbase/hsa.gff3")
# config.txt file used by mirdeep to label fastq files
config_file <- here("data/lists/config.txt")
# genotype and sex data
genotype_file <- here("results/HM3mds/sample_mds_sex_by_dnaid_and_donorid.tsv")
# mirge mirna sequence data
mirge_mirbase_seqs_file <- here("data/mirge2.0/miRge.Libs/human/fasta.Libs/human_mirna_SNP_pseudo_miRBase.fa")
mirge_mirgene_seqs_file <- here("data/mirge2.0/miRge.Libs/human/fasta.Libs/human_mirna_SNP_pseudo_MirGeneDB.fa")

# Load runs data ##################################################################################
# load config.txt file
runs <- read_tsv(config_file, col_names = FALSE)
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

# Load run_summary metadata #######################################################################
# load sequencing run summary file
run_summary <- read_tsv(run_summary_file)

# merge run_summary into runs table
runs <- left_join(runs, run_summary, by = c("rnaid", "sequencing_lane"))
rm(run_summary)

# Load sample_summary metadata ####################################################################
# load rna sample file
sample_summary <- read_tsv(sample_summary_file)

# merge rna sample data with runs data by RNAID to create column data
# remove barcode_id from sample_summary to prevent duplcated columns while merging
colData <- left_join(runs, dplyr::select(sample_summary, -barcode_id), by = "rnaid")
rm(sample_summary, runs)

# load flowcell data
flowcells <- read_tsv(flowcell_file)
flowcells <- dplyr::rename(flowcells, fastq_file = fastq)
flowcells$fastq <- sapply(strsplit(flowcells$fastq_file, "\\."), `[`, 1)
# merge into colData
colData <- left_join(colData, flowcells, by = "fastq")
rm(flowcells)

# load genotype data
genotypes <- read_tsv(genotype_file)
colData <- dplyr::select(colData, -dnaid, -sex)
colData <- left_join(colData, genotypes, by = "donor_id")
rm(genotypes)

# format colData
# assign row names
colData <- as.data.frame(colData)
row.names(colData) <- colData$fastq

# change some character variables to factor
colData <- mutate(colData, sequencing_lane = factor(sequencing_lane))
colData <- mutate(colData, sequencing_date = factor(sequencing_date))
colData <- mutate(colData, sequencing_project = factor(sequencing_project))
colData <- mutate(colData, sequencing_pool = factor(sequencing_pool))
colData <- mutate(colData, sequencing_run = factor(sequencing_run))

colData <- mutate(colData, barcode = factor(barcode))
colData <- mutate(colData, barcode_id = factor(barcode_id))

colData <- mutate(colData, rna_sample_index = factor(rna_sample_index))
colData <- mutate(colData, rna_plate_number = factor(rna_plate_number))
colData <- mutate(colData, rna_plate_well = factor(rna_plate_well))
colData <- mutate(colData, rna_lane = factor(rna_lane))

colData <- mutate(colData, rna_extraction_date = factor(rna_extraction_date))
colData <- mutate(colData, acquisition_date = factor(acquisition_date))

colData <- mutate(colData, tissue_section = factor(tissue_section, levels = c("CW", "CP", "GZ")))
colData <- mutate(colData, flowcell = factor(flowcell))
colData <- mutate(colData, unique_lane = factor(paste(colData$flowcell, colData$sequencing_lane, sep="_")))

colData <- dplyr::select(colData, -seq_id, -name, -ancestry, -not_usable, -source, -freezer_rack_row,
                         -freezer_rack_col, -freezer_box_row, -freezer_box_col, -freezer_cell_row,
                         -freezer_cell_col, -freezer, -freezer_comments, -exp_condition, -incubator,
                         -dna_extraction_date, -dna_purification_method, -dna_concentration, -dna_260_280,
                         -dna_260_230, -dna_source, -dna_extracted_by, -dna_freezer, -dna_freezer_box,
                         -genotype_id)

# for mirge2.0 data, removed duplicated rows
colData <- colData[!duplicated(colData$rnaid),]
colData <- select(colData, -fastq)
row.names(colData) <- colData$rnaid

# load row data sequence information ##############################################################

sequence_mirbase <- readDNAStringSet(mirge_mirbase_seqs_file)
sequence_mirgene <- readDNAStringSet(mirge_mirgene_seqs_file)


mirbase_mirnas_df <- data.frame(sequence_mirbase)
mirgene_mirnas_df <- data.frame(sequence_mirgene)

#mirnas_df$name <- rownames(mirnas_df)
#mirna_df %<>% rename(seq = mirna_seqs)

mirbase_mirnas_df$gc <- as.numeric(letterFrequency(sequence_mirbase, "GC", as.prob = TRUE))
mirbase_mirnas_df$length <- width(sequence_mirbase)
mirgene_mirnas_df$gc <- as.numeric(letterFrequency(sequence_mirgene, "GC", as.prob = TRUE))
mirgene_mirnas_df$length <- width(sequence_mirgene)

rm(sequence_mirbase, sequence_mirgene)

# load counts data ################################################################################

cts_mirbase <- read_tsv(mirbase_counts_file)
cts_mirgene <- read_tsv(mirgene_counts_file)

# Build SummarizedExperiment ######################################################################

m <- as.matrix(cts_mirbase[,-1])
rownames(m) <- cts_mirbase$miRNA
cts_mirbase <- m
rm(m)

m <- as.matrix(cts_mirgene[,-1])
rownames(m) <- cts_mirgene$miRNA
cts_mirgene <- m
rm(m)

# build se for mirbase
rowData <- mirbase_mirnas_df[rownames(cts_mirbase),]
all(rownames(rowData) == rownames(cts_mirbase))

colData <- colData[colnames(cts_mirbase),]
all(rownames(colData) == colnames(cts_mirbase))

se_mirbase <- SummarizedExperiment(assays=SimpleList(counts=cts_mirbase), colData=colData, rowData=rowData)

# build se for mirgenedb
rowData <- mirgene_mirnas_df[rownames(cts_mirgene),]
all(rownames(rowData) == rownames(cts_mirgene))

colData <- colData[colnames(cts_mirgene),]
all(rownames(colData) == colnames(cts_mirgene))

se_mirgene <- SummarizedExperiment(assays=SimpleList(counts=cts_mirgene), colData=colData, rowData=rowData)

# save se objects
saveRDS(se_mirbase, mirbase_output_file)
saveRDS(se_mirgene, mirgene_output_file)
