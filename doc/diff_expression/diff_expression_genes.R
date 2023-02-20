# diff expression of total RNA-seq expression (gene expression)
# across gz/cp and gest. week

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
#library(cqn)

source(here("src/utils/lafferty_utils.R"))
prefix <- paste(format(Sys.time(), "%Y%m%d"), "diff_gene_expression", sep="_")

# OUTPUT FILES ########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/diff_expression/pdfs/")
# write plots
write.plots <- FALSE
# intermediate rdata files (FOR SAVING)
#dds.gzcp.rds <- paste(here("doc/diff_expression/rdata/"), prefix, "_dds_gzcp.rds", sep="")
#dds.gw.rds <- paste(here("doc/diff_expression/rdata/"), prefix, "_dds_gw.rds", sep="")
# data frame for compiling DE evidence
df.rds <- paste(paste(here("doc/diff_expression/rdata/"), prefix, "_df.rds", sep=""))

# INPUT FILES #########################################################################################################
# RangedSummarizedExperiment (rse) with gene expression
rse.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")
# intermediate rdata files (FOR LOADING)
dds.gzcp.rds <- here("doc/diff_expression/rdata/20190506_diff_gene_expression_dds_gzcp.rds")
dds.gw.rds <- here("doc/diff_expression/rdata/20190506_diff_gene_expression_dds_gw.rds")

# Build DDS ###########################################################################################################
# read in rse
rse <- readRDS(rse.rds)

# count threshold
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# gzcp rse
rse.gzcp <- rse[, rse$tissue_section != "CW"]
rse.gzcp$tissue_section <- droplevels(rse.gzcp$tissue_section)
# gest week rse
rse.gw <- rse[, rse$tissue_section == "CW"]
rse.gw$tissue_section <- droplevels(rse.gw$tissue_section)

# build DESeq data set
dds.gzcp <- DESeqDataSet(rse.gzcp, design = ~1)
dds.gw <- DESeqDataSet(rse.gw, design = ~1)

rm(rse, rse.gzcp, rse.gw)

# variance stabilizing tranformed data
vsd.gzcp <- vst(dds.gzcp)
vsd.gw <- vst(dds.gw)

# PCA #################################################################################################################
# pca on uncorrected data
pca.gzcp <- prcomp(t(assay(vsd.gzcp)))
df.gzcp <- data.frame(pca.gzcp$x[,1:5], colData(vsd.gzcp))
percentVar.gzcp <- pca.gzcp$sdev^2 / sum(pca.gzcp$sdev^2) * 100

df.gzcp %>%
  ggplot(aes(PC1, PC2, col=tissue_section, shape=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar.gzcp[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.gzcp[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_manual(values = c("red", "blue"))

pca.gw <- prcomp(t(assay(vsd.gw)))
df.gw <- data.frame(pca.gw$x[,1:5], colData(vsd.gw))
percentVar.gw <- pca.gw$sdev^2 / sum(pca.gw$sdev^2) * 100

df.gw %>%
  ggplot(aes(PC1, PC2, col=gestation_week)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar.gw[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.gw[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "red"))

# batch correct
vsd.gzcp.trans <- vsd.gzcp
assay(vsd.gzcp.trans) <- limma::removeBatchEffect(assay(vsd.gzcp.trans),
                                                  batch = vsd.gzcp.trans$round,
                                                  design = model.matrix(~vsd.gzcp$tissue_section))

vsd.gw.trans <- vsd.gw
assay(vsd.gw.trans) <- limma::removeBatchEffect(assay(vsd.gw.trans),
                                                batch = vsd.gw.trans$round,
                                                batch2 = vsd.gw.trans$rna_purification_method,
                                                design = model.matrix(~vsd.gw.trans$gestation_week))

# pca on batch corrected data
pca.gzcp <- prcomp(t(assay(vsd.gzcp.trans)))
df.gzcp <- data.frame(pca.gzcp$x[,1:5], colData(vsd.gzcp.trans))
percentVar.gzcp <- pca.gzcp$sdev^2 / sum(pca.gzcp$sdev^2) * 100

df.gzcp %>%
  ggplot(aes(PC1, PC2, col=tissue_section, shape=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar.gzcp[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.gzcp[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_manual(values = c("red", "blue"))

pca.gw <- prcomp(t(assay(vsd.gw.trans)))
df.gw <- data.frame(pca.gw$x[,1:5], colData(vsd.gw.trans))
percentVar.gw <- pca.gw$sdev^2 / sum(pca.gw$sdev^2) * 100

df.gw %>%
  ggplot(aes(PC1, PC2, col=gestation_week)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar.gw[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.gw[2], 1), "%)", sep=""),
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "red"))

if(write.plots){ggsave(file.path(pdf_dir,""), width = 7, height = 6, useDingbats = FALSE)}

# Design Equation #################################################################################
# factors for design equation GZ/CP:
# rin
# donor_id
# round
# rna_concentration
# tissue_section
# gestation_week is perfectly correlated with donor

# relabel round factor labels
dds.gzcp$round[dds.gzcp$round == "round2-redo"] <- "round2_redo"
dds.gzcp$round <- factor(dds.gzcp$round)

dds.gzcp$tissue_section <- droplevels(dds.gzcp$tissue_section)

dds.gzcp$donor_id <- factor(dds.gzcp$donor_id)

design(dds.gzcp) <- formula(~ rin + donor_id + rna_concentration + tissue_section)

# factors for design equation gest. week:
# round
# rna_purification_method
# rin
# rna_concentration
# sex
# gestation_week

# relabel round factor labels
dds.gw$round[dds.gw$round == "round2-redo"] <- "round2_redo"
dds.gw$round <- factor(dds.gw$round)

# relabel rna_purification_method factor labels
dds.gw$rna_purification_method <- factor(dds.gw$rna_purification_method)
dds.gw$rna_purification_method_name <- dds.gw$rna_purification_method
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "Trizol w Glycogen"] <- "trizol"

dds.gw$sex <- factor(dds.gw$sex)

# set design equation
design(dds.gw) <- formula(~ round + rna_purification_method_name + rin + rna_concentration + sex + gestation_week)


# Run DESeq2 ##########################################################################################################
# dds.gzcp <- DESeq(dds.gzcp)
# dds.gw <- DESeq(dds.gw)
#
# saveRDS(dds.gzcp, dds.gzcp.rds)
# saveRDS(dds.gw, dds.gw.rds)

# load dds
dds.gzcp <- readRDS(dds.gzcp.rds)
dds.gw <- readRDS(dds.gw.rds)

# Results #############################################################################################################
#res.gzcp <- results(dds.gzcp)
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

#res.gw <- results(dds.gw)
shrunkres.gw <- lfcShrink(dds.gw, coef="gestation_week", type="apeglm")

# plotMA(res.gzcp)
# summary(res.gzcp)
#
# plotMA(res.gw)
# summary(res.gw)

df.res.gzcp <- as.data.frame(shrunkres.gzcp)
df.res.gzcp$gene_id <- rownames(df.res.gzcp)

df.res.gw <- as.data.frame(shrunkres.gw)
df.res.gw$gene_id <- rownames(df.res.gw)

# full join results, replace NAs in log2FoldChange with 0s
df <- full_join(df.res.gzcp, df.res.gw, by = "gene_id", suffix = c(".gzcp", ".gw"))
df$log2FoldChange.gzcp[is.na(df$log2FoldChange.gzcp)] <- 0

# significance threshold
padjusted.significance.threshold <- 0.05

df$sig.gzcp <- df$padj.gzcp < padjusted.significance.threshold
df$sig.gw <- df$padj.gw < padjusted.significance.threshold
# convert NA to FALSE
df$sig.gzcp[is.na(df$sig.gzcp)] <- FALSE
df$sig.gw[is.na(df$sig.gw)] <- FALSE

# rough order for plotting
#df <- df[order(df$padj.gw, decreasing = TRUE),]

# add in mcols data
rows <- data.frame(mcols(dds.gzcp))
rows %<>% dplyr::select(gene_id, gene_name, gene_biotype, symbol, entrezid)

df <- left_join(df, rows, by = "gene_id")
rm(rows, df.res.gw, df.res.gzcp, shrunkres.gw, shrunkres.gzcp)

# # total read counts
# rowcounts.gzcp <- rowSums(counts(dds.gzcp))
# rowcounts.gw <- rowSums(counts(dds.gw))
# all(names(rowcounts.gw) == names(rowcounts.gzcp))
# rowcounts.all <- rowcounts.gzcp + rowcounts.gw
# df.counts <- data.frame(gene_id = names(rowcounts.all), read_counts = rowcounts.all, stringsAsFactors = FALSE)
# df <- left_join(df, df.counts, by = "gene_id")
#
# rm(rowcounts.all, rowcounts.gw, rowcounts.gzcp, df.counts)

# Label Genes by Category #############################################################################################
# label significant in both analyses column
df$sig.both <- NA
df$sig.both[which(df$sig.gw & df$sig.gzcp)] <- "both"
df$sig.both[which(df$sig.gw & !df$sig.gzcp)] <- "gw only"
df$sig.both[which(!df$sig.gw & df$sig.gzcp)] <- "gzcp only"
df$sig.both[which(!df$sig.gw & !df$sig.gzcp)] <- "neither"
df$sig.both <- factor(df$sig.both)

# genes that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
genes.maturation.late <- df$gene_id[which(df$sig.both == "gw only" & df$log2FoldChange.gw > 0)]
# genes that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
genes.maturation.early <- df$gene_id[which(df$sig.both == "gw only" & df$log2FoldChange.gw < 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
genes.neurogenesis.gz <- df$gene_id[which((df$sig.both == "both" | df$sig.both == "gzcp only") & df$log2FoldChange.gzcp > 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
genes.neurogenesis.cp <- df$gene_id[which((df$sig.both == "both" | df$sig.both == "gzcp only") & df$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df$category <- "none"
df$category[which(df$gene_id %in% genes.maturation.late)] <- "MatLate"
df$category[which(df$gene_id %in% genes.maturation.early)] <- "MatEarly"
df$category[which(df$gene_id %in% genes.neurogenesis.gz)] <- "NeuroGZ"
df$category[which(df$gene_id %in% genes.neurogenesis.cp)] <- "NeuroCP"
# factor category variable
df$category <- factor(df$category)

# MA Plot #############################################################################################################

df %>%
  filter(baseMean.gzcp > 0) %>%
  arrange(sig.gzcp) %>%
  ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp)) +
  geom_point(size=1) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="GZ/CP Differential Expression",
       color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=1, color="blue") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue"))

df %>%
  filter(baseMean.gw > 0) %>%
  arrange(sig.gw) %>%
  ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw)) +
  geom_point(size=1) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Gestation Week Differential Expression",
       color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=1, color="blue") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue"))

# Save Data Frame #####################################################################################################
saveRDS(df, df.rds)

# Scratch #############################################################################################################
stop()



