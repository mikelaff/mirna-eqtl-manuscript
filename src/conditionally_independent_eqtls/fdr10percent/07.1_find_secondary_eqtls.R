# find secondary eQTLs:
# for each eGene/emiR (gene that has at least 1 eSNP below the eigenmt-bh corrected p-value)
#   find SNP with lowest p-value, after secondary analysis


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
results.name <- "20200120_mirQTLor_fdr10percent"

# OUTPUT FILES ########################################################################################################
# output directory
output.dir <- paste0(here("results/conditional_eqtls/"), results.name, "/secondary/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# eQTL output file
eqtl.df.rds <- paste0(output.dir, results.name, "_secondary_eQTLs_dataFrame.rds")

# secondary eSNPs to get genotypes using plink
esnps.txt <- paste0(output.dir, results.name, "_secondary_eSNPs.txt")

# emiR output file
#emir.df.rds <- paste0(output.dir, results.name, "_emiRs_dataFrame.rds")

# INPUT FILES ##########################################################################################################
# eigenMT-BH nominal p-value
nominal.p.value.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue_10percent.txt")

# Summarized association results for every variant within 1MB of each expressed miRNA: SECONDARY RESULTS
summarized.results.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/secondary/association_results/compiled/20200120_mirQTLor_fdr10percent_secondary_variants_dataFrame.rds")

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

# nominal p-value from eigenMT-BH proceedure
nominal.p.value <- as.numeric(read_lines(nominal.p.value.txt))

# Find Primary eQTLs ###################################################################################################
# find emiRs (eGenes)
df.results %>%
    group_by(UniName) %>%
    summarise(min.pval = min(P)) %>%
    filter(min.pval <= nominal.p.value) %>%
    pull(UniName) -> emirs

# find secondary eSNPs: lowest P at each emiR, also lowest base position in event of 2+ variants with minimum p-value
df.results %>%
    filter(UniName %in% emirs) %>%
    group_by(UniName) %>%
    filter(P == min(P)) %>%
    filter(BP.hg38 == min(BP.hg38)) -> df.secondary.eqtls

print(paste("Number of emiRs:", length(emirs)))
print(paste("Number of secondary eSNPs:", sum(!duplicated(df.secondary.eqtls$SNP))))
print(paste("Number of secondary eQTLs:", nrow(df.secondary.eqtls)))

# Export Data ##########################################################################################################

write_rds(df.secondary.eqtls, eqtl.df.rds)

write_lines(unique(df.secondary.eqtls$SNP), esnps.txt)


