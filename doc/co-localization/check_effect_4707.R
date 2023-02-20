# check co-localized effect alleles


library(here)
library(dplyr)
library(readr)
library(magrittr)
library(GenomicRanges)
library(ggplot2)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# dir.pdf <- here("doc/co-localization/pdfs/")
#
# # output folder for combined plots
# dir.output <- paste0(dir.pdf, "mir4707/")
# dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GWAS munged datasets
gr.edu.rds <- here("data/gwas_datasets/educational_attainment/educational_attainment.hg38.GRanges.rds")
gr.sa.global.rds <- here("data/gwas_datasets/enigma3/ENIGMA3_Global/enigma3.surface_area.global.hg38.GRanges.rds")

# GWAS colocalizations
df.coloc.edu.rds <- here("results/external_data/educational_attainment/educational_attainment_mirQTL_colocalizations.rds")

# fetal brain mRNA-eQTL clumped with ld buddies
clumped.ldbuddies.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs_and_LDbuddies.rds")
index.snps.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs.rds")

# 1000 Genomes EUR plink files
dir.1kg.eur.plink <- here("data/1000genomes_phase3_hg38/EUR.plink/")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

GWAS.PVAL.THRESHOLD <- 5e-8

# Import ##############

# 4707-results
df.eqtl.results <- as_tibble(read_rds(summarized.results.dataframe.rds))


df.eqtl.results %<>%
    dplyr::filter(UniName == "hsa-mir-4707_hsa-miR-4707-3p")

df.eqtl.results %<>%
    dplyr::select(SNP,
           BETA.mirQTL = BETA,
           P.mirQTL = P,
           UniName,
           CHR,
           BP.hg38,
           REF,
           ALT,
           EFFECT.ALLELE,
           ALT_CTS,
           OBS_CT)

# mRNA-eQTL results
# df.mrna.eqlt <- as_tibble(read_rds(clumped.ldbuddies.rds))
#
# df.mrna.eqlt %<>%
#     dplyr::filter(ENSG == "ENSG00000092036")
#
# df.mrna.eqlt <-
#
# df.mrna.eqlt %<>%
#     select(SNP.mQTL = SNP_A,
#            BETA.mQTL = BETA,
#            P.mQTL = P,
#            ENSG,
#            BP.hg38 = BP.hg38_A)

# educational attainment

gr.edu <- read_rds(gr.edu.rds)

gr.edu <- gr.edu[seqnames(gr.edu) == "chr14"]

df.edu <- as_tibble(gr.edu)

df.edu %<>%
    dplyr::select(BP.hg38 = start,
           RSID.edu = RSID,
           A1.effect.edu = A1.effect,
           A2.edu = A2,
           BETA.edu = Beta,
           SE.edu = SE,
           P.edu = Pval)

# combine
df.results <- left_join(df.eqtl.results, df.edu, by = "BP.hg38")

# global surface area
gr.sa <- read_rds(gr.sa.global.rds)

gr.sa <- gr.sa[seqnames(gr.sa) == "chr14"]

df.sa <- as_tibble(gr.sa)

df.sa %<>%
    dplyr::select(BP.hg38 = start,
           RSID.sa = RSID,
           A1.effect.sa = A1.effect,
           A2.sa = A2,
           BETA.sa = BETA,
           SE.sa = SE,
           P.sa = Pval)

# combine
df.results <- left_join(df.results, df.sa, by = "BP.hg38")





