
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(GenomicRanges)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES ########################################################################################################
dir.pdfs <- here("doc/eigenmt/pdfs/")

# INPUT FILES ##########################################################################################################
# EMMAX with BH 5% correction
bh.eqtl.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eQTLs_dataFrame.rds")
bh.emir.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_emiRs_dataFrame.rds")
bh.esnp.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eSNPs_dataFrame.rds")

# EMMAX with EigenMT-BH 5% correction
eigenmtbh.eqtl.df.rds <- paste0(here("results/eigenmt/"), results.name, "/compiled/", results.name, "_eQTLs_dataFrame.rds")
eigenmtbh.emir.df.rds <- paste0(here("results/eigenmt/"), results.name, "/compiled/", results.name, "_emiRs_dataFrame.rds")
eigenmtbh.esnp.df.rds <- paste0(here("results/eigenmt/"), results.name, "/compiled/", results.name, "_eSNPs_dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# non-overlapping granges of known and novel mirnas
#granges.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# vsd
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Seqinfo object
SEQINFO <- Seqinfo(genome = "hg38")[CHROMS]


# Load Results #########################################################################################################

df.bh.eqtl <- readRDS(bh.eqtl.df.rds)
df.bh.emir <- readRDS(bh.emir.df.rds)
df.bh.esnp <- readRDS(bh.esnp.df.rds)

df.eigenmtbh.eqtl <- readRDS(eigenmtbh.eqtl.df.rds)
df.eigenmtbh.emir <- readRDS(eigenmtbh.emir.df.rds)
df.eigenmtbh.esnp <- readRDS(eigenmtbh.esnp.df.rds)

df.bh.eqtl$MAF <- ( (df.bh.eqtl$A1.HOM.COUNT * 2) + df.bh.eqtl$HET.COUNT ) / df.bh.eqtl$OBS_CT
df.eigenmtbh.eqtl$MAF <- ( (df.eigenmtbh.eqtl$A1.HOM.COUNT * 2) + df.eigenmtbh.eqtl$HET.COUNT ) / df.eigenmtbh.eqtl$OBS_CT

sum(df.eigenmtbh.eqtl$eQTL %in% df.bh.eqtl$eQTL)

df.bh.eqtl %>%
    dplyr::filter(! eQTL %in% df.eigenmtbh.eqtl$eQTL) -> df.notpassing.eqtl



# MAF
pdf(paste0(dir.pdfs, "maf_eqtls.pdf"), width = 7, height = 4)

df.bh.eqtl %>%
    ggplot(aes(x = MAF)) +
    geom_histogram() +
    labs(title = "MAF of eQTLs passing BH 5%")

df.eigenmtbh.eqtl %>%
    ggplot(aes(x = MAF)) +
    geom_histogram() +
    labs(title = "MAF of eQTLs passing EigenMT-BH 5%")

df.notpassing.eqtl %>%
    ggplot(aes(x = MAF)) +
    geom_histogram() +
    labs(title = "MAF of eQTLs passing BH 5% but not passing EigenMT-BH 5%")

dev.off()



