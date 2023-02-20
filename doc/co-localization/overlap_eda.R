# look at mRNA-eQTL and sQTL overlaps


library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(mikelaffr)


# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/co-localization/pdfs/")
dir.create(dir.pdf, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# table of mRNA-eQTL overlaps
mrna.overlap.rds <- here("results/co-localization/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_overlaps_r2at0.8.rds")

# table of sQTL overlaps
sqtl.overlaps.rds <- here("results/co-localization/fetal_brain_sQTL/fetal_brain_sQTL_mirQTL_overlaps_r2at0.8.rds")

# mRNA expression
rse.gene.expression.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# miRNA expression from eQTL analysis
rse.mirna.expression.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import ###############################################################################################################

df.eqtls <- read_rds(eqtls.dataframe.rds)

df.mirna.overlaps <- read_rds(mrna.overlap.rds)

df.sqtl.overlaps <- read_rds(sqtl.overlaps.rds)

rse.gene.expression <- read_rds(rse.gene.expression.rds)

rse.mirna.expression <- read_rds(rse.mirna.expression.rds)


# Host Genes #################















