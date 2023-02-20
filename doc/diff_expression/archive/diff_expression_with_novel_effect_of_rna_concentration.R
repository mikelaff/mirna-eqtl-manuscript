# diff expression gz/cp

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(GGally)
library(limma)
library(reshape2)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/diff_expression/pdfs/")
png_dir <- here("doc/diff_expression/pngs")
# write plots
write_plots <- TRUE
# intermediate rdata files
dds_rds <- here("doc/diff_expression/rdata/dds_gz_cp_novel.rds")

# INPUT FILES #####################################################################################
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# dds with novel
dds_novel_rds <- here("results/rdata_files/tmp_novel.rds")

# Build DDS #######################################################################################
# read in summarized experiment
se <- readRDS(mirbase_se_rds)
# read in novels
dds.novel <- readRDS(dds_novel_rds)

cts <- assay(se)
coldata <- as.data.frame(colData(se))
rowdata <- as.data.frame(rowData(se))

cts.novel <- counts(dds.novel)
coldata.novel <- as.data.frame(colData(dds.novel))
rowdata.novel <- data.frame(names = rownames(dds.novel), sequence_mirbase = NA, gc = NA, length = NA, row.names = rownames(dds.novel))
rowdata.novel <- rowdata.novel[,2:4]

cts <- rbind2(cts, cts.novel)
rowdata <- rbind.data.frame(rowdata, rowdata.novel)

se.all <- SummarizedExperiment(assays=SimpleList(counts=cts), colData=coldata, rowData=rowdata)


# build DESeq data set
dds <- DESeqDataSet(se.all, design = ~1)

# subset for GZ/CP
dds.gzcp <- dds[,dds$tissue_section == "GZ" | dds$tissue_section == "CP"]
# remove GZ/CP samples
dds.gw <- dds[, dds$tissue_section == "CW"]
dds.gzcp$tissue_section <- droplevels(dds.gzcp$tissue_section)
dds.gw$tissue_section <- droplevels(dds.gw$tissue_section)

# remove rows with only zero counts
dds.gzcp <- dds.gzcp[rowSums(counts(dds.gzcp)) > 1, ]
dds.gw <- dds.gw[rowSums(counts(dds.gw)) > 1, ]
dds <- dds[rowSums(counts(dds)) > 1, ]

# estimate size factors
dds.gzcp <- estimateSizeFactors(dds.gzcp)
dds.gw <- estimateSizeFactors(dds.gw)
dds <- estimateSizeFactors(dds)

rm(se, se.all, dds.novel)

# variance stabilizing tranformed data
vsd.gzcp <- varianceStabilizingTransformation(dds.gzcp)
vsd.gw <- varianceStabilizingTransformation(dds.gw)
vsd <- varianceStabilizingTransformation(dds)


# PCA #############################################################################################

# # batch correct for sequencing date
# vsd.gzcp.trans <- vsd.gzcp
# assay(vsd.gzcp.trans) <- limma::removeBatchEffect(assay(vsd.gzcp.trans),
#                                                   batch = vsd.gzcp.trans$sequencing_date,
#                                                   design = model.matrix(~vsd.gzcp.trans$tissue_section))

# remove seq. pool and rna purification method
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             batch2 = vsd$rna_purification_method,
                                             design = model.matrix(~vsd$gestation_week))

pca <- prcomp(t(assay(vsd.trans)), scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, col=tissue_section, shape=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="batch corrected for seq. pool and pur. method",
       color="Tissue\nSection",
       shape="Sequencing\nDate") +
  plotTheme +
  scale_color_manual(values = c("grey50", "blue", "darkorange"))

if(write_plots){ggsave(file.path(png_dir,"pca_batch_corrected_novel.png"), width = 7, height = 6)}

df %>%
  ggplot(aes(PC1, PC2, col=gestation_week, shape=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="batch corrected for seq. pool and pur. method",
       color="Gestation\nWeek",
       shape="Tissue\nSection") +
  plotTheme +
  scale_color_gradientn(colors = c("blue", "grey60", "red"))

if(write_plots){ggsave(file.path(png_dir,"pca_batch_corrected_novel_by_gest_week.png"), width = 7, height = 6)}

df %>%
  ggplot(aes(PC1, PC2, col=log10(rna_concentration), shape=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="batch corrected for seq. pool and pur. method",
       color="Log10(RNA\nConc.)",
       shape="Tissue\nSection") +
  plotTheme +
  scale_color_gradientn(colors = c("blue", "grey60", "red"))

if(write_plots){ggsave(file.path(png_dir,"pca_batch_corrected_novel_by_rna_conc.png"), width = 7, height = 6)}

ggpairs(df,
        columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "rna_concentration"),
        upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
        lower = list(continuous = wrap(my_custom_smooth)),
        diag = list(continuous = wrap(my_custom_density))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, face = "bold"))

if(write_plots){ggsave(file.path(png_dir,"pca_batch_corrected_novel_by_rna_conc_pairs.png"), width = 10, height = 10)}

df %>%
  ggplot(aes(rna_concentration, PC1, col=tissue_section)) +
  geom_point(size=2) +
  labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       x="RNA Concentration",
       caption="batch corrected for seq. pool and pur. method",
       color="Tissue\nSection") +
  plotTheme +
  scale_color_manual(values = c("grey50", "blue", "darkorange"))

df %>%
  ggplot(aes(rna_concentration, PC4, col=tissue_section)) +
  geom_point(size=2) +
  labs(y=paste("PC4 (", round(percentVar[4], 1), "%)", sep=""),
       x="RNA Concentration",
       caption="batch corrected for seq. pool and pur. method",
       color="Tissue\nSection") +
  plotTheme +
  scale_color_manual(values = c("grey50", "blue", "darkorange"))

df %>%
  ggplot(aes(rna_concentration, PC5, col=tissue_section)) +
  geom_point(size=2) +
  labs(y=paste("PC5 (", round(percentVar[5], 1), "%)", sep=""),
       x="RNA Concentration",
       caption="batch corrected for seq. pool and pur. method",
       color="Tissue\nSection") +
  plotTheme +
  scale_color_manual(values = c("grey50", "blue", "darkorange"))


# GZ/CP #################################################################################
# factors for design equation:
# rin
# donor_id
# sequencing_date
# rna_concentration
# gestation_week is perfectly correlated with donor

levels(dds.gzcp$sequencing_date) <- c("2015_07_01", "2015_07_07")

dds.gzcp$donor_id <- factor(dds.gzcp$donor_id)

design(dds.gzcp) <- formula(~ sequencing_date + donor_id + rin + tissue_section)

dds.gzcp <- DESeq(dds.gzcp)

#saveRDS(dds.gzcp, dds_rds)

# load dds
#dds.gzcp <- readRDS(dds_rds)

# Results
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

summary(shrunkres.gzcp)

# MA Plot
shrunkres.gzcp.df <- as.data.frame(shrunkres.gzcp)
shrunkres.gzcp.df$mirna <- rownames(shrunkres.gzcp.df)
shrunkres.gzcp.df$sig <- "not_sig"
shrunkres.gzcp.df$sig[which(shrunkres.gzcp.df$padj < 0.1)] <- "sig"

shrunkres.gzcp.df <- shrunkres.gzcp.df[order(shrunkres.gzcp.df$padj, decreasing = TRUE),]

# add in mirna gc and length data
rows <- data.frame(length=mcols(dds.gzcp)[,"length"], mirna=rownames(dds.gzcp))
rows$novel <- ifelse(is.na(rows$length), "novel", "mirbase")
rows$novel <- factor(rows$novel)
shrunkres.gzcp.df <- left_join(shrunkres.gzcp.df, rows, by = "mirna")

#pdf("~/Desktop/gz_cp_maplot_shrunk_labs.pdf", height = 7, width = 10)

circColor <- "black"
labelColor <- "black"

ggplot(shrunkres.gzcp.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="miRNA Differential Expression GZ/CP (w/o RNA Conc.)") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.gzcp.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(png_dir,"de_gzcp_without_rna_conc.png"), width = 10, height = 6)}

# With RNA Conc.
design(dds.gzcp) <- formula(~ sequencing_date + donor_id + rin + rna_concentration + tissue_section)

dds.gzcp <- DESeq(dds.gzcp)

#saveRDS(dds.gzcp, dds_rds)

# load dds
#dds.gzcp <- readRDS(dds_rds)

# Results
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

summary(shrunkres.gzcp)

# MA Plot
shrunkres.gzcp.df <- as.data.frame(shrunkres.gzcp)
shrunkres.gzcp.df$mirna <- rownames(shrunkres.gzcp.df)
shrunkres.gzcp.df$sig <- "not_sig"
shrunkres.gzcp.df$sig[which(shrunkres.gzcp.df$padj < 0.1)] <- "sig"

shrunkres.gzcp.df <- shrunkres.gzcp.df[order(shrunkres.gzcp.df$padj, decreasing = TRUE),]

# add in mirna gc and length data
rows <- data.frame(length=mcols(dds.gzcp)[,"length"], mirna=rownames(dds.gzcp))
rows$novel <- ifelse(is.na(rows$length), "novel", "mirbase")
rows$novel <- factor(rows$novel)
shrunkres.gzcp.df <- left_join(shrunkres.gzcp.df, rows, by = "mirna")

#pdf("~/Desktop/gz_cp_maplot_shrunk_labs.pdf", height = 7, width = 10)

circColor <- "black"
labelColor <- "black"

ggplot(shrunkres.gzcp.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="miRNA Differential Expression GZ/CP (w RNA Conc.)") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.gzcp.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(png_dir,"de_gzcp_with_rna_conc.png"), width = 10, height = 6)}

# Gest. Week #################################################################################

# variable of interest: gestation_week
# sequencing_pool
# rna_purification_method
# rin

dds <- dds.gw

# relabel sequencing_pool factor labels
dds$sequencing_pool <- factor(dds$sequencing_pool)
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool1"] <- "Pool1"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool2"] <- "Pool2"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool3"] <- "Pool3"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool4"] <- "Pool4"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool5"] <- "Pool5"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool6"] <- "Pool6"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool7"] <- "Pool7"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool8"] <- "Pool8"

# relabel rna_purification_method factor labels
dds$rna_purification_method <- factor(dds$rna_purification_method)
dds$rna_purification_method_name <- dds$rna_purification_method
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "Trizol w Glycogen"] <- "trizol"

# set design equation
design(dds) <- formula(~ sequencing_pool + rna_purification_method_name + rin + gestation_week)

# Run DESeq2
dds <- DESeq(dds)

#saveRDS(dds, dds_rds)

# load dds
#dds <- readRDS(dds_rds)

# Results

# shrunken results
shrunkres <- lfcShrink(dds, coef="gestation_week", type="apeglm")

#plotMA(shrunkres)
summary(shrunkres)

# data frame from shrunken results
shrunkres.df <- as.data.frame(shrunkres)
shrunkres.df$mirna <- rownames(shrunkres.df)
shrunkres.df$sig <- "not_sig"
shrunkres.df$sig[which(shrunkres.df$padj < 0.1)] <- "sig"
# order data frame for plotting
shrunkres.df <- shrunkres.df[order(shrunkres.df$padj, decreasing = TRUE),]

rows <- data.frame(length=mcols(dds)[,"length"], mirna=rownames(dds))
rows$novel <- ifelse(is.na(rows$length), "novel", "mirbase")
rows$novel <- factor(rows$novel)
shrunkres.df <- left_join(shrunkres.df, rows, by = "mirna")

# create rank of log2foldchange
shrunkres.df$log2FoldChange.rank <- rank(shrunkres.df$log2FoldChange)


# MA Plot
circColor <- "black"
labelColor <- "black"

ggplot(shrunkres.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="miRNA Differential Expression GW (w/o RNA Conc.)") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(png_dir,"de_gw_without_rna_conc.png"), width = 10, height = 6)}

# set design equation
design(dds) <- formula(~ sequencing_pool + rna_purification_method_name +
                         rin + rna_concentration + gestation_week)

# Run DESeq2
dds <- DESeq(dds)

#saveRDS(dds, dds_rds)

# load dds
#dds <- readRDS(dds_rds)

# Results

# shrunken results
shrunkres <- lfcShrink(dds, coef="gestation_week", type="apeglm")

#plotMA(shrunkres)
summary(shrunkres)

# data frame from shrunken results
shrunkres.df <- as.data.frame(shrunkres)
shrunkres.df$mirna <- rownames(shrunkres.df)
shrunkres.df$sig <- "not_sig"
shrunkres.df$sig[which(shrunkres.df$padj < 0.1)] <- "sig"
# order data frame for plotting
shrunkres.df <- shrunkres.df[order(shrunkres.df$padj, decreasing = TRUE),]

rows <- data.frame(length=mcols(dds)[,"length"], mirna=rownames(dds))
rows$novel <- ifelse(is.na(rows$length), "novel", "mirbase")
rows$novel <- factor(rows$novel)
shrunkres.df <- left_join(shrunkres.df, rows, by = "mirna")

# create rank of log2foldchange
shrunkres.df$log2FoldChange.rank <- rank(shrunkres.df$log2FoldChange)


# MA Plot
circColor <- "black"
labelColor <- "black"

ggplot(shrunkres.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="miRNA Differential Expression GW (w RNA Conc.)") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(png_dir,"de_gw_with_rna_conc.png"), width = 10, height = 6)}






















