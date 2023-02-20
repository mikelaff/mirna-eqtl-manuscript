# look at mirQTL distances and genomic context

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_eda/pdfs/")

# INPUT ################################################################################################################
# eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold
nom.p.val.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################
# compiled eqtls
df.eqtls <- readRDS(eqtls.dataframe.rds)












