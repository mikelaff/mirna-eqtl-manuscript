# co-localize mirQTL eSNPs with fetal brain mRNA-eQTLs

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(mikelaffr)


# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/co-localization/pdfs/")

# table of colocalizations
coloc.output.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalizations.rds")

# list of colocalized snps for ld calculations
coloc.snps.mQTL.output.txt <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalized_snps_mQTL.txt")
coloc.snps.mirQTL.output.txt <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalized_snps_mirQTL.txt")

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# mirQTL eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
esnps.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# mirQTL emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Directory for LD at each mirQTL index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# fetal brain mRNA-eQTL clumped with ld buddies
clumped.ldbuddies.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs_and_LDbuddies.rds")
index.snps.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs.rds")

# summarized experiment with mirna vst expression values used in this analysis
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import ###############################################################################################################
# mirQTL eqtls
df.mirQTL.eqtls <- readRDS(eqtls.dataframe.rds)

# mrna eqtl ld buddies
df.mQTL.ldbuddies <- readRDS(clumped.ldbuddies.rds)

# r2 > 0.8 ld buddies
df.mQTL.high.ldbuddies <- dplyr::filter(df.mQTL.ldbuddies, R2 > 0.8)

# label columns with mQTL
colnames(df.mQTL.high.ldbuddies) <- paste0(colnames(df.mQTL.high.ldbuddies), ".mQTL")

# position ID for matching to mirQTL ld buddies
df.mQTL.high.ldbuddies %<>%
    mutate(POSID_B = paste(CHR_B.mQTL, BP.hg38_B.mQTL, sep = "_"))

rm(df.mQTL.ldbuddies)

# miRNA expression data
vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd.mirna)]

# set seqnames for gene expression to UCSC
seqlevelsStyle(rse.gene) <- "UCSC"
genome(rse.gene) <- "hg38"
# remove non-standard chroms
seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

# Colocalizations ######################################################################################################


# loop over each eqtl (index SNP) and get ld buddies.
# check overlap of mirQTL ld buddies with mQTL ld buddies at r2 > 0.8
df.coloc <- tibble()

for (i in 1:nrow(df.mirQTL.eqtls)) {

    df.overlap <- NULL

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, df.mirQTL.eqtls$SNP.CHR[i], ".hardcall.prefiltered.mirQTLor.", df.mirQTL.eqtls$eSNP[i], ".ld")
    # ld for this index SNP
    df.mirQTL.ldbuddies <- read_table(ld.file, col_types = cols())
    df.mirQTL.ldbuddies %<>%
        mutate(CHR_A = paste0("chr", CHR_A),
               CHR_B = paste0("chr", CHR_B))
    df.mirQTL.ldbuddies %<>%
        dplyr::rename(BP.hg38_A = BP_A,
                      BP.hg38_B = BP_B)
    df.mirQTL.ldbuddies %<>%
        dplyr::mutate(eQTL = df.mirQTL.eqtls$eQTL[i],
                      emiR = df.mirQTL.eqtls$emiR[i],
                      eSNP = df.mirQTL.eqtls$eSNP[i],
                      BETA = df.mirQTL.eqtls$BETA[i],
                      P = df.mirQTL.eqtls$P[i])
    colnames(df.mirQTL.ldbuddies) <- paste0(colnames(df.mirQTL.ldbuddies), ".mirQTL")
    # position ID for matching to mQTL ld buddies
    df.mirQTL.ldbuddies %<>%
        mutate(POSID_B = paste(CHR_B.mirQTL, BP.hg38_B.mirQTL, sep = "_"))

    # r2 > 0.8
    df.mirQTL.high.ldbuddies <- dplyr::filter(df.mirQTL.ldbuddies, R2.mirQTL > 0.8)
    rm(df.mirQTL.ldbuddies)

    # ld buddy overlap by position id
    if (any(df.mirQTL.high.ldbuddies$POSID_B %in% df.mQTL.high.ldbuddies$POSID_B)) {
        df.overlap <- inner_join(df.mirQTL.high.ldbuddies, df.mQTL.high.ldbuddies, by = "POSID_B", suffix = c(".mirQTL", ".mQTL"))
        df.coloc <- bind_rows(df.coloc, df.overlap)
    }
}



length(unique(df.coloc$eQTL.mirQTL))

length(unique(df.coloc$SNP_ENSG.mQTL))

length(unique(df.coloc$eSNP.mirQTL))

# Expression Correlation ###############################################################################################

df.coloc$emiR_eGene_Corr <- NA

for (k in 1:nrow(df.coloc)) {

    emiR.expression <- assay(vsd.mirna)[df.coloc$emiR.mirQTL[k],]

    if (!df.coloc$ENSG.mQTL[k] %in% rownames(vsd.gene)) {
        df.coloc$emiR_eGene_Corr[k] <- 0
        next
    } else {
        ensg.expression <- assay(vsd.gene)[df.coloc$ENSG.mQTL[k],]
    }
    stopifnot(all(names(ensg.expression) == names(emiR.expression)))

    df.coloc$emiR_eGene_Corr[k] <- cor(emiR.expression, ensg.expression)

}

df.coloc %>%
    dplyr::filter(!duplicated(emiR_eGene_Corr)) %>%
    ggplot(aes(x = emiR_eGene_Corr)) +
    geom_histogram(bins = 10) +
    labs(title = "emiR-mRNA Expression Correlation for Co-Localized eQTLs")

ggsave(paste0(dir.pdf, "emir-mrna_expression_correlation_for_coloc_eqtls.pdf"))

# Export Table #########################################################################################################
# save colocalization data frame
saveRDS(df.coloc, coloc.output.rds)

# list of colocalized snps for ld calculations
write_lines(unique(df.coloc$SNP_B.mirQTL), coloc.snps.mirQTL.output.txt)
write_lines(unique(df.coloc$SNP_B.mQTL), coloc.snps.mQTL.output.txt)

