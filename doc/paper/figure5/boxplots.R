# plot boxplots of miRNA and mRNA expression
# mir 4707, HAUS4

library(here)
library(dplyr)
library(readxl)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(mikelaffr)

# OUTPUT ###############################################################################################################


# INPUT ################################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# genotypes directory
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# gene expression rse
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# vsd for the mirQTL
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")


# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

eQTL <- "hsa-mir-4707_hsa-miR-4707-3p+chr14:22953244:A:G"
emiR <- "hsa-mir-4707_hsa-miR-4707-3p"
eSNP <- "chr14:22953244:A:G"
chr <- "chr14"

rsid <- "rs4981455"
ensg <- "ENSG00000092036"
# Load RSE #############################################################################################################

# mirna expression
vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd.mirna)]

# # set seqnames for gene expression to UCSC
# seqlevelsStyle(rse.gene) <- "UCSC"
# genome(rse.gene) <- "hg38"
# # remove non-standard chroms
# seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS

# threshold for expression
rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

rm(rse.gene, dds.gene)


# miRNA eQTLs ##########################################################################################################
df.mirqtls <- read_rds(eqtls.dataframe.rds)



# get genotypes at this snp
df.genotypes <- read_table(paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", eSNP, ".raw"))
df.genotypes %<>%
    dplyr::select(DonorID = 1,
                  DNAID = 2,
                  genotypes_number_of_G = 7)


samples <- as_tibble(colData(vsd.gene))
# samples %<>%
#   mutate(donor_id = paste0("D", donor_id))

samples %<>%
    left_join(df.genotypes, by = "DonorID")

samples %<>%
    mutate(genotypes_number_of_G = factor(genotypes_number_of_G, levels = c(0,1,2), labels = c("AA", "AG", "GG"), ordered = TRUE))

samples %<>%
    dplyr::rename(genotype = genotypes_number_of_G)




# RNAID of genotypes
rnaid.aa <- samples$RNAID[samples$genotype == "AA"]
rnaid.ag <- samples$RNAID[samples$genotype == "AG"]
rnaid.gg <- samples$RNAID[samples$genotype == "GG"]







expression.vst <- as.data.frame(t(assay(vsd.mirna[emiR,], 1)))
colnames(expression.vst) <- "expression.vst"
expression.vst$RNAID <- rownames(expression.vst)

expression.vst.residual <- as.data.frame(t(assay(vsd.mirna[emiR,], 2)))
colnames(expression.vst.residual) <- "expression.vst.residual"
expression.vst.residual$RNAID <- rownames(expression.vst.residual)

expression.mrna <- as.data.frame(t(assay(vsd.gene[ensg,], 1)))
colnames(expression.mrna) <- "expression.haus4"
expression.mrna$RNAID <- rownames(expression.mrna)


df.expression <- as_tibble(left_join(expression.vst, expression.vst.residual, by = "RNAID"))
df.expression %<>%
    left_join(expression.mrna, by = "RNAID")

df.expression %<>%
    left_join(samples, by = "RNAID")


df.expression %>%
    ggplot(aes(x = genotype, y = expression.vst.residual)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/mir-4707-3p_boxplot.pdf", height = 2, width = 2)

df.expression %>%
    ggplot(aes(x = genotype, y = expression.vst)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/mir-4707-3p_vst_boxplot.pdf", height = 2, width = 2)

df.expression %>%
    ggplot(aes(x = genotype, y = expression.haus4)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/haus4_boxplot.pdf", height = 2, width = 2)

