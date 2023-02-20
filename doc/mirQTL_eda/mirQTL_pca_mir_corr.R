# does miR-101-3p correlate with any covariates?

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(pheatmap)
library(reshape2)
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

# use only miRBase expression for pca
vsd.mirbase <- vsd[mcols(vsd)$source == "miRBase_v22",]

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

# Expression PCA #######################################################################################################
expr.pca <- prcomp(t(assay(vsd.mirbase)))
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

loadings.pc1 <- sort(abs(expr.pca$rotation[,1]), decreasing = TRUE)[1:20]
loadings.pc2 <- sort(abs(expr.pca$rotation[,2]), decreasing = TRUE)[1:20]
loadings.pc3 <- sort(abs(expr.pca$rotation[,3]), decreasing = TRUE)[1:20]

df.expr <- as.data.frame(t(assay(vsd)))
df.expr$RNAID <- rownames(df.expr)

df.combo <- left_join(df.expr.pca, df.expr, by = "RNAID")

mat.combo <- as.matrix(df.combo[,c(1:10, 12:918)])

# correlation matrix
cormat <- cor(mat.combo, use = "pair")

mat.plotting <- cormat[1:10, 11:ncol(cormat)]

df.corrs <- melt(mat.plotting)

df.corrs %>%
    group_by(Var2) %>%
    summarise(max = max(abs(value))) %>%
    ggplot(aes(x = max)) +
    geom_histogram() +
    labs(y = "Number of miRs",
         x = "Abs(Corr. to any PC)",
         title = "Maximum Corr. between miR Expression and Any PC (PC1-10)")

ggsave(paste0(dir.pdfs, "number_mirs_max_corr_to_pcs.pdf"))

#sort(colSums(abs(mat.plotting)), decreasing = TRUE)[1:20]

column.maxs <- colMaxs(abs(mat.plotting))
names(column.maxs) <- colnames(mat.plotting)

sort(column.maxs, decreasing = TRUE)[1:100]

mat.plotting <- mat.plotting[,names(sort(column.maxs, decreasing = TRUE)[1:100])]


pdf(paste0(dir.pdfs, "top_100_mirna_expression_corr_to_pcs.pdf"), height = 5, width = 16)
pheatmap(mat.plotting,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         main = "Top 100 miRNA Expression Correlations to PCs")
dev.off()







