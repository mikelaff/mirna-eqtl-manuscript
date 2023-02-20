# look at diff gene expression between mir-4707 expressing and non-expressing donors

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)
library(DESeq2)

# OUTPUT FILES #########################################################################################################


# INPUT FILES ##########################################################################################################
# samples file
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# ranged summarized experiment for mirna expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# ranged summarized experiment for gene expression values
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# non-overlapping granges of known and novel mirnas
granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Autosomes covariate matrix
cov.mat.auto.file <- here("results/emmax/covariates/20200120_chrAutosomes.mirQTLor.columnLabeled.cov")

# directory with primary eqtl genotypes
dir.primary.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor/primary/primary_eSNP_genotypes/")

# GLOBALS ##############################################################################################################

# Import Samples #######################################################################################################
samples <- read_delim(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# Import RSE ###########################################################################################################
rse <- readRDS(rse.gene.rds)





