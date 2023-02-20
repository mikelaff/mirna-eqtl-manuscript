# Select samples for miRNA-eQTL analysis
# Write sample/dnaid list as plink readible file
# mirQTLor are 212 donors after outliers from expression PCA have been removed, see:
# src/genotypes/mirQTLor_imputed_genotypes/

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

date.prefix <- format(Sys.time(), "%Y%m%d")

# OUTPUT FILES ########################################################################################################
# output directory for sample information
output.dir <- here("results/emmax/samples/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# PLINK file for --keep filtering
plink.keep.txt <- paste0(output.dir, date.prefix, "_mirQTLor_DonorID_DNAID.txt")

# RNAIDs of usable cortical wall small-rna-seq samples along with corresponding DonorID and DNAID
rnaids.txt <- paste0(output.dir, date.prefix, "_mirQTLor_RNAID_DonorID_DNAID.txt")

# INPUT FILES #########################################################################################################
# this should be a database query, however, here is the table of all small rna-seq samples
samples.small.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")
samples.total.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_totalRNAseq_metadata.tsv")

# plink .fam file to get sexes by genotype
samples.fam <- here("results/genotypes/AllSamplesQC.fam")

# GLOBALS ##############################################################################################################
# label outliers by PCA (after batch effect removal of pool and rna purif. method, see below)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")

# Select Samples ######################################################################################################
# Select samples for fetal tissue cortical wall small-rna-seq eQTL analysis

# load small-rna-seq metadata
df.small.samples <- read_tsv(samples.small.tsv)

# subset for cortical wall samples
df.small.samples %<>%
    filter(TissueSection == "CW")

# load total-rna-seq metadata
df.total.samples <- read_tsv(samples.total.tsv)

# subset for cortical wall samples
df.total.samples %<>%
    filter(TissueSection == "CW")

# import .fam file
df.fam <- read.table(samples.fam, col.names = c("DonorID", "DNAID", "ID3", "ID4", "Sex", "Pheno"), stringsAsFactors = FALSE)
df.fam$index <- 1:length(df.fam$DonorID)

# convert sex and select columns
df.sex <- select(df.fam, DNAID, Sex, index)
df.sex$Sex <- ifelse(df.fam$Sex == 1, "M", ifelse(df.fam$Sex == 2, "F", NA))
df.sex %<>%
    select(DNAID, Sex.by.Genotype = Sex, index)

# tmp <- left_join(df.small.samples, df.total.samples, by = "RNAID", suffix = c(".small", ".total"))
# tmp %<>%
#     select(RNAID, starts_with("DonorID"), starts_with("VerifiedDNAID"))

# filter for samples with a VarifiedDNAID and no DonorID mismatches
df.small.samples %<>%
    filter(VerifiedDNAID != "MISSING" & VerifiedDNAID != "MIXTURE",
           DonorID.MISMATCH == FALSE)

df.total.samples %<>%
    filter(VerifiedDNAID != "MISSING" & VerifiedDNAID != "MIXTURE",
           DonorID.MISMATCH == FALSE)

# filter for small rna-seq samples in total rna-seq samples
df.small.samples %<>%
    filter(RNAID %in% df.total.samples$RNAID)

# join with genotype sexes
df.small.samples %<>%
    left_join(df.sex, by = c("VerifiedDNAID" = "DNAID"))

# label sex mismaches
df.small.samples$Sex.MISMATCH <- ! df.small.samples$Sex.by.XIST == df.small.samples$Sex.by.Genotype

# remove sex mismatches
df.small.samples %<>%
    filter(!Sex.MISMATCH)

# remove outliers
df.small.samples %<>%
    filter(!RNAID %in% outlier_rnaid)

# arrange by DonorID (this is how .fam files are organized)
df.small.samples %<>%
    arrange(index)

stopifnot(sum(duplicated(df.small.samples$RNAID)) == 0,
          sum(duplicated(df.small.samples$VerifiedDNAID)) == 0,
          sum(duplicated(df.small.samples$DonorID)) == 0)

# PLINK: Filtering with --keep ########################################################################################
#--keep accepts a space/tab-delimited text file with family IDs (DonorID) in the first column and within-family
# IDs (DNAID) in the second column, and removes all unlisted samples from the current analysis. --remove does the
# same for all listed samples.

write_lines(paste(df.small.samples$DonorID, df.small.samples$VerifiedDNAID), plink.keep.txt)


# Export RNAIDs #######################################################################################################

write_lines(paste(df.small.samples$RNAID, df.small.samples$DonorID, df.small.samples$VerifiedDNAID), rnaids.txt)

