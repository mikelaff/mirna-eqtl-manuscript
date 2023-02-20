# build gene expression ranged summarized experiment for fetal tissue samples

# this RSE is missing phenotype information (colData)

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)

# OUTPUT FILES ####################################################################################
output_file <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_fetalTissue_ranged_summarized_experiment_gene_counts.rds", sep="")

# INPUT FILES #####################################################################################
counts_file <- here("results/gene_counts/20180905_fetalTissue_gene_counts.tsv")
samples_file <- here("data/metadata/fetal_tissue_rnaseq.csv")
samples_file2 <- here("data/metadata/gene_expression_metadata.tsv")
rna_file <- here("data/metadata/rna_sample_summary.tsv")
# genotype and sex data
genotype_file <- here("results/HM3mds/sample_mds_sex_by_dnaid_and_donorid.tsv")
# male rnaids that are not labeled by genotype data
male_rnaids_file <- here("data/metadata/male_unlabeled_rnaids.txt")
# female rnaids that are not labeled by gneotype data
female_rnaids_file <- here("data/metadata/female_unlabeled_rnaids.txt")

# Row Data ########################################################################################
# the ensembl human database
edb <- EnsDb.Hsapiens.v86
rowranges <- genes(edb)
rm(edb)

# Column Data #####################################################################################
# create phenotype information (placeholder until more phenotype information is added)
samples <- read_csv(samples_file)
rounds <- read_tsv(samples_file2)
samples %<>% dplyr::select(rnaid, section)
rounds %<>% dplyr::select(rnaid, round)
samples$section[is.na(samples$section)] <- "CW"
samples$section[samples$section == "VZ"] <- "GZ"
samples %<>% dplyr::rename(tissue_section = section)
samples <- samples[!duplicated(samples$rnaid),]

samples <- left_join(samples, rounds, by = "rnaid")
rm(rounds)

# rna data
rna_data <- read_tsv(rna_file)
samples <- left_join(samples, dplyr::select(rna_data, -tissue_section), by = "rnaid")
rm(rna_data)

# load genotype data
genotypes <- read_tsv(genotype_file)
samples <- dplyr::select(samples, -dnaid, -sex)
samples <- left_join(samples, genotypes, by = "donor_id")
rm(genotypes)

# read male and female rnaids
males <- read_lines(male_rnaids_file)
females <- read_lines(female_rnaids_file)
# label unlabeled males and females
samples$sex[which(samples$rnaid %in% males)] <- "Male"
samples$sex[which(samples$rnaid %in% females)] <- "Female"
rm(males, females)

samples <- as.data.frame(samples)
rownames(samples) <- samples$rnaid

# Expression Data #################################################################################
counts <- read_tsv(counts_file)
# count matrix
cts <- as.matrix(dplyr::select(counts, -ENSG))
rownames(cts) <- counts$ENSG
rm(counts)

# Create RSE ######################################################################################
# order rowranges on rownames of cts
rowranges <- rowranges[rownames(cts)]
stopifnot(all(names(rowranges) == rownames(cts)))

# order count matrix
cts <- cts[, rownames(samples)]
stopifnot(all(rownames(samples) == colnames(cts)))

# build RangedSummarizedExperiment
rse <- SummarizedExperiment(assays = list(counts = cts),
                            rowRanges = rowranges,
                            colData = samples)

# save rse to rds file
saveRDS(rse, output_file)
