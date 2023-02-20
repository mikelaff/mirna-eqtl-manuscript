# does miR-101-3p correlate with any covariates?

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_eda/pdfs/")
#dir.pngs <- here("doc/mirQTL_eda/pngs/")

# INPUT ################################################################################################################
# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20191204_mirQTLor_VST_miRNA_expression/20191204_mirQTLor_VST_miRNA_expression_rse.rds")

# Covariate matrix autosomes
cov.mat.auto.file <- here("results/emmax/covariates/20191204_chrAutosomes.mirQTLor.columnLabeled.cov")

# GLOBALS ##############################################################################################################

# Import Expression Data ###############################################################################################

# expression data
vsd <- readRDS(vsd.rds)

# Import Covariate Matrix ##############################################################################################

# autosome covariates file
cov.mat.auto <- read_delim(cov.mat.auto.file, delim = " ")
cov.mat.auto %<>%
    dplyr::rename(DonorID = X1, DNAID = X2)

# covariates
batchVars.auto <- c(paste0("PC", 1:10, ".genotype"),
                    paste0("PC", 1:10, ".expression"),
                    "PoolPool2", "PoolPool3", "PoolPool4",
                    "PoolPool5", "PoolPool6", "PoolPool7",
                    "PoolPool8","PurificationMethodmiRNeasy",
                    "PurificationMethodmiRNeasy_mini", "SexM", "RIN", "GestationWeek")

# formula strings
#formula.string.auto <- paste("expression ~", paste(batchVars.auto, collapse = " + "))

# Plot Correlation Heatmap #############################################################################################

# mir-101 expression
df.mir101.expr <- data.frame(mir_101_3p = assay(vsd)["hsa-mir-101-1_hsa-miR-101-3p",], rnaid = names(assay(vsd)["hsa-mir-101-1_hsa-miR-101-3p",]), stringsAsFactors = FALSE)

# format with donor id
df.mir101.expr <- dplyr::select(left_join(df.mir101.expr, as.data.frame(colData(vsd)), by = "rnaid"), mir_101_3p, DonorID = donor_id)
df.mir101.expr$DonorID <- paste0("D", df.mir101.expr$DonorID)

# combine with covariate matrix
df.cov <- left_join(cov.mat.auto, df.mir101.expr, by = "DonorID")
mat.cov <- as.matrix(df.cov[,4:36])

# correlation matrix
cormat <- cor(mat.cov, use = "pair")

# plot correlation matrix
pdf(paste0(dir.pdfs, "cov_mat_to_mir101_expr_heatmap.pdf"))
pheatmap(cormat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Expression vs PCs ####################################################################################################

pdf(paste0(dir.pdfs, "mir101_expr_pc3.pdf"))
df.cov %>%
    ggplot(aes(x = mir_101_3p, y = PC3.expression)) +
    geom_point() +
    plotTheme() +
    labs(x = "miR-101-3p Expression (VST Norm.)")
dev.off()

cor(df.cov$mir_101_3p, df.cov$PC3.expression)

df.cov %>%
    ggplot(aes(x = mir_101_3p, y = PC2.expression)) +
    geom_point()

df.cov %>%
    ggplot(aes(x = mir_101_3p, y = PC1.expression)) +
    geom_point()

df.cov %>%
    ggplot(aes(x = mir_101_3p, y = PC3.expression, color = GestationWeek)) +
    geom_point()

