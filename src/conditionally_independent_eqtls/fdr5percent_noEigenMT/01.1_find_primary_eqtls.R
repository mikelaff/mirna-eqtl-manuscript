# find primary eQTLs:
# for each eGene/emiR (gene that has at least 1 eSNP below the eigenmt-bh corrected p-value)
#   find SNP with lowest p-value


library(here)
library(readr)
library(dplyr)
library(magrittr)
#library(DESeq2)
#library(GenomicRanges)
library(mikelaffr)

date.prefix <- format(Sys.time(), "%Y%m%d")
date.prefix <- "20200120"

# association results directory name
results.name <- "20200120_mirQTLor_fdr5percent_noEigenMT"

# OUTPUT FILES ########################################################################################################
# output directory
output.dir <- paste0(here("results/conditional_eqtls/"), results.name, "/primary/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# eQTL output file
eqtl.df.rds <- paste0(output.dir, results.name, "_primary_eQTLs_dataFrame.rds")

# primary eSNPs to get genotypes using plink
esnps.txt <- paste0(output.dir, results.name, "_primary_eSNPs.txt")

# emiR output file
#emir.df.rds <- paste0(output.dir, results.name, "_emiRs_dataFrame.rds")

# INPUT FILES ##########################################################################################################
# nominal p-value
nominal.p.value.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# vsd
#vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################
# CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Seqinfo object
#SEQINFO <- Seqinfo(genome = "hg38")[CHROMS]

# Import Summarized Data ###############################################################################################
# summarizedd results at each variant
df.results <- as_tibble(read_rds(summarized.results.dataframe.rds))

# nominal p-value
nominal.p.value <- as.numeric(read_lines(nominal.p.value.txt))

# Find Primary eQTLs ###################################################################################################
# find emiRs (eGenes)
df.results %>%
    group_by(UniName) %>%
    summarise(min.pval = min(P)) %>%
    filter(min.pval <= nominal.p.value) %>%
    pull(UniName) -> emirs

# find primary eSNPs: lowest P at each emiR, also lowest base position in event of 2+ variants with minimum p-value
df.results %>%
    filter(UniName %in% emirs) %>%
    group_by(UniName) %>%
    filter(P == min(P)) %>%
    filter(BP.hg38 == min(BP.hg38)) -> df.primary.eqtls

print(paste("Number of emiRs:", length(emirs)))
print(paste("Number of primary eSNPs:", sum(!duplicated(df.primary.eqtls$SNP))))
print(paste("Number of primary eQTLs:", nrow(df.primary.eqtls)))

# Export Data ##########################################################################################################

write_rds(df.primary.eqtls, eqtl.df.rds)

write_lines(unique(df.primary.eqtls$SNP), esnps.txt)


