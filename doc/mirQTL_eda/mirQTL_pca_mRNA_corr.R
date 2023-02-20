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
# ranged summarized experiment of mRNA expression
rse.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# summarized experiment with mirna vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20191204_mirQTLor_VST_miRNA_expression/20191204_mirQTLor_VST_miRNA_expression_rse.rds")

# GLOBALS ##############################################################################################################

# Import Expression Data ###############################################################################################

# mirna data
vsd.mirna <- readRDS(vsd.rds)

# expression data
rse <- readRDS(rse.rds)

# build DESeq data set
dds <- DESeqDataSet(rse, design = ~1)
# threshold counts
dds <- dds[rowSums(counts(dds) > 10) > 10, ]
# filter for samples within mirna dataset
dds <- dds[,which(colnames(dds) %in% vsd.mirna$rnaid)]

rm(rse, vsd.mirna)

# variance stabilizing transformation
vsd <- vst(dds)

# Expression PCA #######################################################################################################
expr.pca <- prcomp(t(assay(vsd)))
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

loadings.pc1 <- sort(abs(expr.pca$rotation[,1]), decreasing = TRUE)[1:20]
loadings.pc2 <- sort(abs(expr.pca$rotation[,2]), decreasing = TRUE)[1:20]
loadings.pc3 <- sort(abs(expr.pca$rotation[,3]), decreasing = TRUE)[1:20]

df.expr <- as.data.frame(t(assay(vsd)))
df.expr$RNAID <- rownames(df.expr)

df.combo <- left_join(df.expr.pca, df.expr, by = "RNAID")

mat.combo <- as.matrix(df.combo[,c(1:10, 12:ncol(df.combo))])

# correlation matrix
cormat <- cor(mat.combo, use = "pair")

saveRDS(cormat, here("results/rdata_files/20200103_tmp_cormat.rds"))

mat.plotting <- cormat[1:10, 11:ncol(cormat)]

df.corrs <- melt(mat.plotting)

df.corrs %>%
    group_by(Var2) %>%
    summarise(max = max(abs(value))) %>%
    ggplot(aes(x = max)) +
    geom_histogram() +
    labs(y = "Number of mRNAs",
         x = "Abs(Corr. to any PC)",
         title = "Maximum Corr. between mRNA Expression and Any PC (PC1-10)")

ggsave(paste0(dir.pdfs, "number_mrnas_max_corr_to_pcs.pdf"))

#sort(colSums(abs(mat.plotting)), decreasing = TRUE)[1:20]

column.maxs <- colMaxs(abs(mat.plotting))
names(column.maxs) <- colnames(mat.plotting)

sort(column.maxs, decreasing = TRUE)[1:100]

mat.plotting2 <- mat.plotting[,names(sort(column.maxs, decreasing = TRUE)[1:500])]


pdf(paste0(dir.pdfs, "top_500_mrna_expression_corr_to_pcs.pdf"), height = 5, width = 16)
pheatmap(mat.plotting2,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = "Top 500 mRNA Expression Correlations to PCs")
dev.off()







