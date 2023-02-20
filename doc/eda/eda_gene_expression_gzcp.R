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

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/eda/pdfs/")
# write plots
write_plots <- FALSE

# INPUT FILES #####################################################################################
# ranged summarized experiment with gene expression counts
expression_rse_rds <- here("results/rdata_files/20181018_fetalTissue_ranged_summarized_experiment_gene_counts.rds")
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")

# miRNA data #######################################################################################
# read in summarized experiment
se <- readRDS(mirbase_se_rds)

# Build DDS gene expression data #######################################################################################
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
