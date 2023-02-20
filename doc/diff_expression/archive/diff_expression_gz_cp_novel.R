# diff expression gz/cp

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(cqn)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/diff_expression/pdfs/")
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


#dds <- collapseReplicates(dds, groupby = dds$rnaid)
# subset for GZ/CP
dds.gzcp <- dds[,dds$tissue_section == "GZ" | dds$tissue_section == "CP"]
# remove GZ/CP samples
dds <- dds[, dds$tissue_section == "CW"]
dds.gzcp$tissue_section <- droplevels(dds.gzcp$tissue_section)
# remove rows with only zero counts
dds.gzcp <- dds.gzcp[rowSums(counts(dds.gzcp)) > 1, ]
# estimate size factors
dds.gzcp <- estimateSizeFactors(dds.gzcp)
rm(se)

# variance stabilizing tranformed data
vsd.gzcp <- varianceStabilizingTransformation(dds.gzcp)

# PCA #############################################################################################
# pca on uncorrected gz/cp data
pca <- prcomp(t(assay(vsd.gzcp)))
df <- data.frame(pca$x[,1:5], colData(vsd.gzcp))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, col=tissue_section, shape=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, col=sequencing_pool, shape=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme


# batch correct for sequencing date
vsd.gzcp.trans <- vsd.gzcp
assay(vsd.gzcp.trans) <- limma::removeBatchEffect(assay(vsd.gzcp.trans),
                                                  batch = vsd.gzcp.trans$sequencing_date,
                                                  design = model.matrix(~vsd.gzcp.trans$tissue_section))

pca <- prcomp(t(assay(vsd.gzcp.trans)))
df <- data.frame(pca$x[,1:5], colData(vsd.gzcp.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, col=tissue_section, shape=sequencing_date)) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="VST normalized, corrected for sequencing_date",
       color="Tissue\nSection",
       shape="Sequencing\nDate") +
  plotTheme +
  scale_color_manual(values = c("blue", "darkorange"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_gz_cp_batch_corrected_novel.pdf"), width = 7, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, col=factor(donor_id), shape=sex)) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="VST normalized, corrected for sequencing_date",
       color="Donor",
       shape="Sex") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(4,3,7)])

df %>%
  ggplot(aes(PC1, PC2, col=gestation_week, shape=tissue_section)) +
  geom_point(size=3) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="VST normalized, corrected for sequencing_date",
       color="Gestation\nWeek",
       shape="Tissue\nSection") +
  plotTheme +
  scale_color_gradientn(colors = c("blue", "grey50", "red"),breaks=c(17,18,19))

ggplot(df, aes(PC1, PC2, col=sequencing_pool, shape=sequencing_date)) +
  geom_point(size=2) +
  plotTheme

df$donor_id <- factor(df$donor_id)
ggplot(df, aes(PC1, PC2, col=donor_id, shape=tissue_section)) +
  geom_point(size=2) +
  plotTheme

# Design Equation #################################################################################
# factors for design equation:
# rin
# donor_id
# sequencing_date
# rna_concentration
# gestation_week is perfectly correlated with donor

levels(dds.gzcp$sequencing_date) <- c("2015_07_01", "2015_07_07")

dds.gzcp$donor_id <- factor(dds.gzcp$donor_id)

design(dds.gzcp) <- formula(~ sequencing_date + donor_id + rin + tissue_section)

# Run DESeq2 ######################################################################################
dds.gzcp <- DESeq(dds.gzcp)

#saveRDS(dds.gzcp, dds_rds)

# load dds
dds.gzcp <- readRDS(dds_rds)

# Results #########################################################################################
res.gzcp <- results(dds.gzcp)
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

plotMA(res.gzcp)
summary(res.gzcp)

plotMA(shrunkres.gzcp)
summary(shrunkres.gzcp)

# MA Plot #########################################################################################
res.gzcp.df <- as.data.frame(res.gzcp)
res.gzcp.df$mirna <- rownames(res.gzcp.df)
res.gzcp.df$sig <- "not_sig"
res.gzcp.df$sig[which(res.gzcp.df$padj < 0.1)] <- "sig"

res.gzcp.df <- res.gzcp.df[order(res.gzcp.df$padj, decreasing = TRUE),]

# add in mirna gc and length data
rows <- data.frame(mcols(dds.gzcp)[,c("sequence_mirbase", "gc", "length")], mirna=rownames(dds.gzcp))
res.gzcp.df <- left_join(res.gzcp.df, rows, by = "mirna")

#pdf("~/Desktop/gz_cp_maplot_hsas.pdf", height = 7, width = 10)

ggplot(res.gzcp.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Log2 Fold Change") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=1, color="blue") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue"))

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
       title="miRNA Differential Expression GZ/CP") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.gzcp.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

ggsave("~/Desktop/gzcp_de_novel.pdf", height = 6, width = 11, useDingbats=FALSE)

ggplot(shrunkres.gzcp.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change",
         title="miRNA Differential Expression GZ/CP") +
    geom_hline(yintercept = 0, size=.5, color="black") +
    theme(legend.position = "none") +
    plotTheme +
    scale_x_log10(limits = c(1,10000)) +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_point(data=subset(shrunkres.gzcp.df, novel == "novel" & sig == "sig"),
               aes(baseMean, log2FoldChange),
               shape=21, size=4, stroke=2, color=circColor) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "miRge_78_arm3"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "miRge_28_arm5"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "miRge_6_arm5"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=-.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "chr10_22893_mat"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=-.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "chr2_6912_mat"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=-.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "chr10_23310_mat"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=-.5) +
    geom_label(data=subset(shrunkres.gzcp.df, mirna == "miRge_111_arm3"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=-.5)

ggsave("~/Desktop/gzcp_de_novel_labels.pdf", height = 4, width = 5, useDingbats=FALSE)

dplyr::filter(shrunkres.gzcp.df, sig == "sig", novel == "novel")

#ggsave("~/Desktop/shrunken_l2fc_novel.pdf", width = 11, height = 6, useDingbats = FALSE)

#dev.off()

#pdf("~/Desktop/gz_cp_maplot_shrunk_labs2.pdf", height = 7, width = 10)

ggplot(shrunkres.gzcp.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=1, color="blue") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue")) +
  geom_label(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29a-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5) +
  geom_label(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29a-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5) +
  geom_label(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29c-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5) +
  geom_label(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29c-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5) +
  geom_label(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5) +
  geom_point(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29a-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color="darkorange") +
  geom_point(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29a-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color="darkorange") +
  geom_point(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29c-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color="darkorange") +
  geom_point(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29c-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color="darkorange") +
  geom_point(data=subset(shrunkres.gzcp.df, mirna == "hsa-miR-29b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color="darkorange")

#dev.off()

ggsave(plot=p1,height=6,width=10,dpi=300, filename="~/Desktop/ma_labels.pdf", useDingbats=FALSE)
 #ggsave(plot=p1,height=6,width=5.5,dpi=300, filename="~/Desktop/pca_tissue.pdf", useDingbats=FALSE)

ggplot(shrunkres.gzcp.df, aes(x=gc, y=log2FoldChange)) +
  geom_point(aes(col=sig)) +
  geom_smooth(method="lm") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue"))

ggplot(res.gzcp.df, aes(x=length, y=log2FoldChange)) +
  geom_point(aes(col=sig)) +
  geom_smooth(method="lm") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue"))

# tmp <- as.data.frame(results(dds.gzcp))
# tmp$mirna <- rownames(tmp)
# tmp %<>% select(mirna, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj)
# write_csv(tmp, here("results/20180327_gz_cp_mirna_diff_expression.csv"))

# cqn
cts <- counts(dds)
fit <- cqn(cts, x=mcols(dds)$gc, lengths=mcols(dds)$length, sizeFactors = dds$sizeFactor, verbose = TRUE, lengthMethod = "fixed")

cqnplot(fit, n=1)
#cqnplot(fit, n=2)

# cqn
cts <- counts(dds.gzcp)
fit <- cqn(cts, x=mcols(dds.gzcp)$gc, lengths=mcols(dds.gzcp)$length, sizeFactors = dds.gzcp$sizeFactor, verbose = TRUE, lengthMethod = "fixed")

cqnplot(fit, n=1, col=dds.gzcp$tissue_section, xlab = "GC content", lty = 1, main="tissue section")
legend("bottom", levels(dds.gzcp$tissue_section), fill=1:nlevels(dds.gzcp$tissue_section))

cqnplot(fit, n=1, col=dds.gzcp$sequencing_date, xlab = "GC content", lty = 1, main="sequencing date")
legend("bottom", levels(dds.gzcp$sequencing_date), fill=1:nlevels(dds.gzcp$sequencing_date))

cqnplot(fit, n=1, col=dds.gzcp$sequencing_date, xlab = "GC content", lty = 1, main="sequencing date")
legend("bottom", levels(dds.gzcp$sequencing_date), fill=1:nlevels(dds.gzcp$sequencing_date))
