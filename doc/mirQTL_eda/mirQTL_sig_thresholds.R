
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
# library(DESeq2)
# library(GenomicRanges)
# library(RColorBrewer)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_eda/pdfs/")
dir.create(dir.pdfs, showWarnings = FALSE, recursive = TRUE)

# INPUT ################################################################################################################
# eQTLs eigenMT-BH 5%
eqtls.eigenMT.fdr5percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_eQTLs_dataFrame.rds")
nomPvalue.eigenMT.fdr5percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_nomPvalue.txt")

# eQTLs eigenMT-BH 10%
eqtls.eigenMT.fdr10percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/20200120_mirQTLor_fdr10percent_conditional_eQTLs_dataFrame.rds")
nomPvalue.eigenMT.fdr10percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/20200120_mirQTLor_fdr10percent_conditional_nomPvalue.txt")

# eQTLs FDR 5%, no eigenMT
eqtls.fdr5percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_conditional_eQTLs_dataFrame.rds")
nomPvalue.fdr5percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_conditional_nomPvalue.txt")





# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Results ############################################################################################

df.eqtls.eigenMT.fdr5 <- read_rds(eqtls.eigenMT.fdr5percent.rds)
df.eqtls.eigenMT.fdr10 <- read_rds(eqtls.eigenMT.fdr10percent.rds)
df.eqtls.fdr5 <- read_rds(eqtls.fdr5percent.rds)

sum(duplicated(df.eqtls.fdr5$eQTL))

all(df.eqtls.eigenMT.fdr5$eQTL %in% df.eqtls.eigenMT.fdr10$eQTL)

all(df.eqtls.eigenMT.fdr5$eQTL %in% df.eqtls.fdr5$eQTL)

all(df.eqtls.eigenMT.fdr10$eQTL %in% df.eqtls.fdr5$eQTL)


df.eqtls.eigenMT.fdr5$SIGNIFICANCE <- "eigenMT_fdr5percent"
df.eqtls.eigenMT.fdr10$SIGNIFICANCE <- "eiganMT_fdr10percent"
df.eqtls.fdr5$SIGNIFICANCE <- "fdr5percent"

df.eqtls.eigenMT.fdr5 %>%
    bind_rows(filter(df.eqtls.eigenMT.fdr10, ! eQTL %in% df.eqtls.eigenMT.fdr5$eQTL)) -> df.eqtls

df.eqtls %<>%
    bind_rows(filter(df.eqtls.fdr5, ! eQTL %in% df.eqtls$eQTL))

sum(duplicated(df.eqtls$eQTL))





