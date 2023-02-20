# diff expression gestation week

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/diff_expression/pdfs/")
# write plots
write_plots <- TRUE
# intermediate rdata files
dds_rds <- here("doc/diff_expression/rdata/dds_gw_novel.rds")
#dds_outliers_removed_rds <- here("doc/diff_expression/rdata/dds_outliers_removed_gw_novel.rds")
dds_gz_cp_rds <- here("doc/diff_expression/rdata/dds_gz_cp_novel.rds")

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


# remove GZ/CP samples
dds <- dds[, dds$tissue_section == "CW"]
dds$tissue_section <- droplevels(dds$tissue_section)
# remove rows with only zero counts
dds <- dds[rowSums(counts(dds)) > 1, ]
# estimate size factors
dds <- estimateSizeFactors(dds)
rm(se)

# Design Equation #################################################################################

# variable of interest: gestation_week
# sequencing_pool
# rna_purification_method
# rin

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

# Run DESeq2 ######################################################################################
dds <- DESeq(dds)

saveRDS(dds, dds_rds)

# load dds
dds <- readRDS(dds_rds)

# Results #########################################################################################
res <- results(dds)
# shrunken results
shrunkres <- lfcShrink(dds, coef="gestation_week", type="apeglm")

#plotMA(res)
summary(res)

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


# MA Plot #########################################################################################
circColor <- "black"
labelColor <- "black"

ggplot(shrunkres.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue")) +
  geom_point(data=subset(shrunkres.df, novel == "novel" & sig == "sig"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

ggsave("~/Desktop/shrunken_l2fc_gw_novel.pdf", width = 11, height = 6, useDingbats = FALSE)

tmp <- dplyr::filter(shrunkres.df, sig == "sig", novel == "novel")


ggplot(shrunkres.df, aes(x=baseMean, y=log2FoldChange, col=sig)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change") +
    geom_hline(yintercept = 0, size=.5, color="black") +
    theme(legend.position = "none") +
    plotTheme +
    scale_x_log10(limits = c(1,10000)) +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_point(data=subset(shrunkres.df, novel == "novel" & sig == "sig"),
               aes(baseMean, log2FoldChange),
               shape=21, size=4, stroke=2, color=circColor) +
    geom_label(data=subset(shrunkres.df, mirna == "chr6_15151_mat"),
               aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5)

ggsave("~/Desktop/gw_de_novel_labels.pdf", height = 4, width = 5, useDingbats=FALSE)

####
pdfs <- list.files(here("doc/coverage_plots/pdfs/novel/sort/all/"))

df <- data.frame(files = pdfs, stringsAsFactors = FALSE)
df$trunc <- sapply(strsplit(df$files, "\\."), `[`, 1)
df$one <- sapply(strsplit(df$trunc, "_"), `[`, 1)
df$two <- sapply(strsplit(df$trunc, "_"), `[`, 2)
df$three <- sapply(strsplit(df$trunc, "_"), `[`, 3)

df$Derives_from <- NA
df$Derives_from[1:564] <- paste(df$two[1:564], df$three[1:564], sep="_")
df$Derives_from[565:620] <- paste("miRge", df$two[565:620], sep="_")

tmp <- dplyr::filter(shrunkres.df, sig == "sig", novel == "novel")
tmp$one <- sapply(strsplit(tmp$mirna, "_"), `[`, 1)
tmp$two <- sapply(strsplit(tmp$mirna, "_"), `[`, 2)
tmp$Derives_from <- paste(tmp$one, tmp$two, sep = "_")


df <- dplyr::filter(df, Derives_from %in% tmp$Derives_from)

file.copy(from=paste(here("doc/coverage_plots/pdfs/novel/sort/all/"), df$files, sep=""), to=here("doc/coverage_plots/pdfs/novel/sort/de_gw_sig/"))

####

shrunkres.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across Gestation Week") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  presTheme +
  scale_color_manual(values=c("grey70", "red")) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(shrunkres.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"diff_expression_gw_shrunklfc_labels.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Outliers Removed ################################################################################
# only keep non outliers
dds.sub <- dds[, dds$outlier == "no"]
# estimate size factors
dds.sub <- estimateSizeFactors(dds.sub)

# run DESeq2
#dds.sub <- DESeq(dds.sub)
#saveRDS(dds.sub, dds_outliers_removed_rds)

# load dds
dds.sub <- readRDS(dds_outliers_removed_rds)
# shrunken results
shrunkres.sub <- lfcShrink(dds.sub, coef="gestation_week", type="apeglm")

# data frame from shrunken results
shrunkres.sub.df <- as.data.frame(shrunkres.sub)
shrunkres.sub.df$mirna <- rownames(shrunkres.sub.df)
shrunkres.sub.df$sig <- "not_sig"
shrunkres.sub.df$sig[which(shrunkres.sub.df$padj < 0.1)] <- "sig"
# order data frame for plotting
shrunkres.sub.df <- shrunkres.sub.df[order(shrunkres.sub.df$padj, decreasing = TRUE),]

# create rank of log2foldchange
shrunkres.sub.df$log2FoldChange.rank <- rank(shrunkres.sub.df$log2FoldChange)

# MA plot
circColor <- "black"
labelColor <- "black"

shrunkres.sub.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across Gestation Week (outliers removed)") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  presTheme +
  scale_color_manual(values=c("grey70", "red")) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(shrunkres.sub.df, mirna == "hsa-let-7b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(shrunkres.sub.df, mirna == "hsa-let-7b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"diff_expression_gw_shrunklfc_labels_outliers_removed_let7b.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Compare w/ and w/o Outliers #####################################################################

combined_shrunk_results_df <- left_join(shrunkres.df, shrunkres.sub.df, suffix = c(".w_outliers", ".outliers_removed"), by = "mirna")

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.w_outliers, y=log2FoldChange.outliers_removed)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log10(baseMean.w_outliers), y=log10(baseMean.outliers_removed), color=sig.outliers_removed)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.rank.w_outliers, y=log2FoldChange.rank.outliers_removed, color=sig.w_outliers)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

# Compare to GZ/CP Diff. Expression ###############################################################
dds.gzcp <- readRDS(dds_gz_cp_rds)
# shrunken results
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

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

# create rank of log2foldchange
shrunkres.gzcp.df$log2FoldChange.rank <- rank(shrunkres.gzcp.df$log2FoldChange)

combined_shrunk_results_df <- full_join(shrunkres.df, shrunkres.gzcp.df, suffix = c(".gw", ".gzcp"), by = "mirna")

#combined_shrunk_results_df <- full_join(shrunkres.sub.df, shrunkres.gzcp.df, suffix = c(".gw", ".gzcp"), by = "mirna")

combined_shrunk_results_df$sig.gw <- ifelse(is.na(combined_shrunk_results_df$log2FoldChange.gw), "not_sig", combined_shrunk_results_df$sig.gw)
combined_shrunk_results_df$sig.gzcp <- ifelse(is.na(combined_shrunk_results_df$log2FoldChange.gzcp), "not_sig", combined_shrunk_results_df$sig.gzcp)

combined_shrunk_results_df$baseMean.gw <- ifelse(is.na(combined_shrunk_results_df$baseMean.gw), 0, combined_shrunk_results_df$baseMean.gw)
combined_shrunk_results_df$log2FoldChange.gw <- ifelse(is.na(combined_shrunk_results_df$log2FoldChange.gw), 0, combined_shrunk_results_df$log2FoldChange.gw)
combined_shrunk_results_df$baseMean.gzcp <- ifelse(is.na(combined_shrunk_results_df$baseMean.gzcp), 0, combined_shrunk_results_df$baseMean.gzcp)
combined_shrunk_results_df$log2FoldChange.gzcp <- ifelse(is.na(combined_shrunk_results_df$log2FoldChange.gzcp), 0, combined_shrunk_results_df$log2FoldChange.gzcp)

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gw)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gzcp)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log10(baseMean.gw), y=log10(baseMean.gzcp))) +
  geom_point() +
  plotTheme

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.rank.gw, y=log2FoldChange.rank.gzcp, color=sig.gw)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

saveRDS(combined_shrunk_results_df, here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_novel_df.rds"))


