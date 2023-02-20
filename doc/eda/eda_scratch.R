# eda on gene expression of fetal tissue samples
# looking only at gz/cp

library(here)
library(dplyr)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(GGally)
library(limma)
library(reshape2)
library(mikelaff)

# OUTPUT FILES ########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/eda/pdfs/")

# INPUT FILES #########################################################################################################
# RangedSummarizedExperiment (rse) with gene expression
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")
# RangedSummarizedExperiment (rse) with known and novel counts. mirbaseV22 counts quantified by mirge2.0
# (mirbaseV22, friedlander, nowakowski, mirdeep2, mirge2.0)
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GLOBALS #############################################################################################################
WRITE.PLOTS <- TRUE

# miRNA Data ##########################################################################################################
# read in ranged summarized experiment
rse.mirna <- readRDS(rse.mirna.rds)

# Gene Expression Data ################################################################################################
rse.gene <- readRDS(rse.gene.rds)


genes <- rownames(rse.gene)

genes.proteinCoding <- rownames(rse.gene[rowData(rse.gene)$gene_biotype == "protein_coding",])

exprssed.rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10,]

genes.expressed <- rownames(exprssed.rse.gene)

genes.expressed.proteinCoding <- rownames(exprssed.rse.gene[rowData(exprssed.rse.gene)$gene_biotype == "protein_coding",])

df.gene.diff.expression.rds <- here("doc/diff_expression/rdata/20190519_diff_gene_expression_df.rds")
df.de.gene <- readRDS(df.gene.diff.expression.rds)
df.de.gene %<>%
  dplyr::filter(sig.gzcp | sig.gw)

genes.diffExpressed <- df.de.gene$gene_id

genes.diffExpressed.proteinCoding <- df.de.gene[df.de.gene$gene_biotype == "protein_coding",]$gene_id

genes.diffExpressed.NeuroGZ <- df.de.gene[df.de.gene$category == "NeuroGZ",]$gene_id

genes.diffExpressed.NeuroCP <- df.de.gene[df.de.gene$category == "NeuroCP",]$gene_id

genes.diffExpressed.MatEarly <- df.de.gene[df.de.gene$category == "MatEarly",]$gene_id

genes.diffExpressed.MatLate <- df.de.gene[df.de.gene$category == "MatLate",]$gene_id

genes.diffExpressed.NeuroGZ.proteinCoding <- df.de.gene[df.de.gene$category == "NeuroGZ" & df.de.gene$gene_biotype == "protein_coding",]$gene_id

genes.diffExpressed.NeuroCP.proteinCoding <- df.de.gene[df.de.gene$category == "NeuroCP" & df.de.gene$gene_biotype == "protein_coding",]$gene_id

genes.diffExpressed.MatEarly.proteinCoding <- df.de.gene[df.de.gene$category == "MatEarly" & df.de.gene$gene_biotype == "protein_coding",]$gene_id

genes.diffExpressed.MatLate.proteinCoding <- df.de.gene[df.de.gene$category == "MatLate" & df.de.gene$gene_biotype == "protein_coding",]$gene_id













# Build DDS gene expression data ######################################################################################
# read in ranged summarized experiment
rse <- readRDS(expression_rse_rds)

# build DESeq data set
dds <- DESeqDataSet(rse, design = ~1)
# remove rows with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]
# filter for samples within mirna dataset
dds <- dds[,which(colnames(dds) %in% se$rnaid)]

# subset for gz/cp samples
dds <- dds[,which(dds$tissue_section == "GZ" | dds$tissue_section == "CP") ]
# estimate size factors for normalization
dds <- estimateSizeFactors(dds)
# norm. transformed counts (DESeq size factors, log2 transformed)
#ntc <- assay(normTransform(dds))
# variance stabilizing tranformed data
vsd <- vst(dds)
rm(rse, se)

# PCA #################

pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
  geom_point() +
  geom_line() +
  labs(x="PC", y="Percent Total Variance",
       title="PCA: Gene Expression (VST Norm.)") +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, color=tissue_section)) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  scale_color_manual(values=cbPalette) +
  labs(title="Gene Expression (VST Norm.)")

df %>%
  ggplot(aes(PC1, PC2, color=round, shape=tissue_section)) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  scale_color_manual(values=cbPalette) +
  labs(title="Gene Expression (VST Norm.)")

df %>%
  ggplot(aes(PC1, PC2, color=round, shape=factor(donor_id))) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  scale_color_manual(values=cbPalette) +
  labs(title="Gene Expression (VST Norm.)")

df %>%
  ggplot(aes(PC1, PC2, color=sex, shape=factor(donor_id))) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  scale_color_manual(values=cbPalette) +
  labs(title="Gene Expression (VST Norm.)")
