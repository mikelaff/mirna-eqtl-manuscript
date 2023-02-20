# diff expression of known and novel mirnas
# across gz/cp and gest. week

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
#library(cqn)

source(here("src/utils/lafferty_utils.R"))
prefix <- paste(format(Sys.time(), "%Y%m%d"), "diff_expression", sep="_")

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
# RangedSummarizedExperiment (rse) with known and novel counts. mirbaseV22 counts quantified by mirge2.0
# (mirbaseV22, friedlander, nowakowski, mirdeep2, mirge2.0)
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# intermediate rdata files (FOR LOADING)
dds.gzcp.rds <- here("doc/diff_expression/rdata/20190506_diff_expression_dds_gzcp.rds")
dds.gw.rds <- here("doc/diff_expression/rdata/20190506_diff_expression_dds_gw.rds")

# Build DDS ###########################################################################################################
# read in rse
rse <- readRDS(rse.rds)

# count threshold
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# remove novels overlapping mirbase mirnas
rse.novels <- rse[mcols(rse)$source != "miRBase_v22",]
rse.mirbase <- rse[mcols(rse)$source == "miRBase_v22",]
hits <- findOverlaps(rse.novels, rse.mirbase,
                     ignore.strand=TRUE)
rse.novels <- rse.novels[! 1:length(rse.novels) %in% queryHits(hits),]
rse <- rbind(rse.mirbase, rse.novels)

# # non overlapping ranges
# hits <- findOverlaps(rse, ignore.strand=TRUE,
#                      drop.self=TRUE, drop.redundant=TRUE)
# rse <- rse[! 1:length(rse) %in% subjectHits(hits),]

# gzcp rse
rse.gzcp <- rse[, rse$tissue_section != "CW"]
rse.gzcp$tissue_section <- droplevels(rse.gzcp$tissue_section)
# gest week rse
rse.gw <- rse[, rse$tissue_section == "CW"]
rse.gw$tissue_section <- droplevels(rse.gw$tissue_section)

# build DESeq data set
dds.gzcp <- DESeqDataSet(rse.gzcp, design = ~1)
dds.gw <- DESeqDataSet(rse.gw, design = ~1)

rm(rse, rse.gzcp, rse.gw, rse.novels, rse.mirbase)

# variance stabilizing tranformed data
vsd.gzcp <- varianceStabilizingTransformation(dds.gzcp)
vsd.gw <- varianceStabilizingTransformation(dds.gw)

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
                                                  batch = vsd.gzcp.trans$sequencing_date,
                                                  design = model.matrix(~vsd.gzcp.trans$tissue_section))

vsd.gw.trans <- vsd.gw
assay(vsd.gw.trans) <- limma::removeBatchEffect(assay(vsd.gw.trans),
                                                batch = vsd.gw.trans$sequencing_pool,
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


# Design Equation #################################################################################
# factors for design equation GZ/CP:
# rin
# donor_id
# sequencing_date
# rna_concentration
# gestation_week is perfectly correlated with donor

levels(dds.gzcp$sequencing_date) <- c("2015_07_01", "2015_07_07")

dds.gzcp$donor_id <- factor(dds.gzcp$donor_id)

design(dds.gzcp) <- formula(~ sequencing_date + donor_id + rin + rna_concentration + tissue_section)

# factors for design equation gest. week:
# sequencing_pool
# rna_purification_method
# rin
# rna_concentration
# gestation_week

# relabel sequencing_pool factor labels
dds.gw$sequencing_pool <- factor(dds.gw$sequencing_pool)
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool1"] <- "Pool1"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool2"] <- "Pool2"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool3"] <- "Pool3"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool4"] <- "Pool4"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool5"] <- "Pool5"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool6"] <- "Pool6"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool7"] <- "Pool7"
levels(dds.gw$sequencing_pool)[levels(dds.gw$sequencing_pool) == "2015-9087Pool8"] <- "Pool8"

# relabel rna_purification_method factor labels
dds.gw$rna_purification_method <- factor(dds.gw$rna_purification_method)
dds.gw$rna_purification_method_name <- dds.gw$rna_purification_method
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
levels(dds.gw$rna_purification_method_name)[levels(dds.gw$rna_purification_method_name) == "Trizol w Glycogen"] <- "trizol"

# set design equation
design(dds.gw) <- formula(~ sequencing_pool + rna_purification_method_name + rin + rna_concentration + gestation_week)


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

#plotMA(res.gzcp)
#summary(res.gzcp)

#plotMA(res.gw)
#summary(res.gw)

df.res.gzcp <- as.data.frame(shrunkres.gzcp)
df.res.gzcp$Name <- rownames(df.res.gzcp)

df.res.gw <- as.data.frame(shrunkres.gw)
df.res.gw$Name <- rownames(df.res.gw)

# full join results, replace NAs in log2FoldChange with 0s
df <- full_join(df.res.gzcp, df.res.gw, by = "Name", suffix = c(".gzcp", ".gw"))
df$log2FoldChange.gzcp[is.na(df$log2FoldChange.gzcp)] <- 0

# significance threshold
padjusted.significance.threshold <- 0.1

df$sig.gzcp <- df$padj.gzcp < padjusted.significance.threshold
df$sig.gw <- df$padj.gw < padjusted.significance.threshold
# convert NA to FALSE
df$sig.gzcp[is.na(df$sig.gzcp)] <- FALSE
df$sig.gw[is.na(df$sig.gw)] <- FALSE

# rough order for plotting
#df <- df[order(df$padj.gw, decreasing = TRUE),]

# add in mcols data
rows <- data.frame(mcols(dds.gzcp))
rows$Name <- rownames(rows)
rows %<>% dplyr::select(Name, ID, Alias, Derives_from, source, type, score, sequence)

df <- left_join(df, rows, by = "Name")
rm(rows, df.res.gw, df.res.gzcp, shrunkres.gw, shrunkres.gzcp)

# total read counts
rowcounts.gzcp <- rowSums(counts(dds.gzcp))
rowcounts.gw <- rowSums(counts(dds.gw))
all(names(rowcounts.gw) == names(rowcounts.gzcp))
rowcounts.all <- rowcounts.gzcp + rowcounts.gw
df.counts <- data.frame(Name = names(rowcounts.all), read_counts = rowcounts.all, stringsAsFactors = FALSE)
df <- left_join(df, df.counts, by = "Name")

rm(rowcounts.all, rowcounts.gw, rowcounts.gzcp, df.counts)

# Label miRNAs by Category ############################################################################################
# label significant in both analyses column
df$sig.both <- NA
df$sig.both[which(df$sig.gw & df$sig.gzcp)] <- "both"
df$sig.both[which(df$sig.gw & !df$sig.gzcp)] <- "gw only"
df$sig.both[which(!df$sig.gw & df$sig.gzcp)] <- "gzcp only"
df$sig.both[which(!df$sig.gw & !df$sig.gzcp)] <- "neither"
df$sig.both <- factor(df$sig.both)

# mirnas that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
miRNAs.maturation.late <- df$Name[which(df$sig.both == "gw only" & df$log2FoldChange.gw > 0)]
# mirnas that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
miRNAs.maturation.early <- df$Name[which(df$sig.both == "gw only" & df$log2FoldChange.gw < 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
miRNAs.neurogenesis.gz <- df$Name[which((df$sig.both == "both" | df$sig.both == "gzcp only") & df$log2FoldChange.gzcp > 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
miRNAs.neurogenesis.cp <- df$Name[which((df$sig.both == "both" | df$sig.both == "gzcp only") & df$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df$category <- "none"
df$category[which(df$Name %in% miRNAs.maturation.late)] <- "MatLate"
df$category[which(df$Name %in% miRNAs.maturation.early)] <- "MatEarly"
df$category[which(df$Name %in% miRNAs.neurogenesis.gz)] <- "NeuroGZ"
df$category[which(df$Name %in% miRNAs.neurogenesis.cp)] <- "NeuroCP"
# factor category variable
df$category <- factor(df$category)

# MA Plot #############################################################################################################

df %>%
  filter(baseMean.gzcp > 0) %>%
  arrange(sig.gzcp) %>%
  ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="GZ/CP Differential Expression",
       color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=1, color="blue") +
  presTheme +
  scale_color_manual(values=c("grey70", "blue"))

if(write.plots){ggsave(paste(dir.pdf, prefix, "_gzcp_miRNA.pdf", sep = ""), width = 10, height = 6, useDingbats = FALSE)}

df %>%
  filter(baseMean.gw > 0) %>%
  arrange(sig.gw) %>%
  ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Gestation Week Differential Expression",
       color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=1, color="blue") +
  plotTheme +
  scale_color_manual(values=c("grey70", "blue"))

if(write.plots){ggsave(paste(dir.pdf, prefix, "_gw_miRNA.pdf", sep = ""), width = 10, height = 6, useDingbats = FALSE)}

# colored for miRNA source
df %>%
  filter(baseMean.gzcp > 0) %>%
  arrange(sig.gzcp) %>%
  ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, fill=sig.gzcp)) +
  geom_point(shape=21, stroke=0, size=2.5) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="GZ/CP Differential Expression",
       fill=paste("p.adj >", padjusted.significance.threshold, sep=""),
       color="Novel Source") +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=.7, color="black") +
  plotTheme +
  scale_fill_manual(values=c("grey80", "blue")) +
  geom_point(data = filter(df, sig.gzcp, source != "miRBase_v22"),
             mapping = aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color = source, fill=NULL),
             shape=21, size=4, stroke=2) +
  scale_color_manual(values = cbPalette[c(2,4,5,6)])

if(write.plots){ggsave(paste(dir.pdf, prefix, "_gzcp_miRNA_source.pdf", sep = ""), width = 10, height = 6, useDingbats = FALSE)}


df %>%
  filter(baseMean.gw > 0) %>%
  arrange(sig.gw) %>%
  ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, fill=sig.gw)) +
  geom_point(shape=21, stroke=0, size=2.5) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Gestation Week Differential Expression",
       fill=paste("p.adj >", padjusted.significance.threshold, sep=""),
       color="Novel Source") +
  scale_x_log10() +
  theme(legend.position = "bottom") +
  geom_hline(yintercept = 0, size=.7, color="black") +
  plotTheme +
  scale_fill_manual(values=c("grey80", "blue")) +
  geom_point(data = filter(df, sig.gw, source != "miRBase_v22"),
             mapping = aes(x=baseMean.gw, y=log2FoldChange.gw, color = source, fill=NULL),
             shape=21, size=4, stroke=2) +
  scale_color_manual(values = cbPalette[c(2,4,5,6)])

if(write.plots){ggsave(paste(dir.pdf, prefix, "_gw_miRNA_source.pdf", sep = ""), width = 10, height = 6, useDingbats = FALSE)}


# Save Data Frame #####################################################################################################
saveRDS(df, df.rds)

# Scratch #############################################################################################################
stop()

df %>%
    filter(baseMean.gzcp > 100, baseMean.gw > 100) %>%
    ggplot(aes(x = log2FoldChange.gw, y = log2FoldChange.gzcp, fill = sig.both)) +
    geom_point(shape = 21, stroke = 0, size = 2)

pdf("~/Desktop/miRNA_endog_control_candidates.pdf", width = 10, height = 10)

df %>%
    filter(baseMean.gw > 0) %>%
    arrange(sig.gw) %>%
    ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change",
         title="Gestation Week Differential Expression",
         color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_label(data=subset(df, Name == "hsa-miR-124-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-124-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-92b-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-92b-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-361-5p"), aes(label=Name), hjust=0, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-361-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-186-5p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-186-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-26a-5p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-26a-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-191-5p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-191-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-423-5p"), aes(label=Name), hjust=1, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-423-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-451a"), aes(label=Name), hjust=0, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-451a"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-320a-3p"), aes(label=Name), hjust=1, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-320a-3p"), shape=21, size=4, stroke=2, color="black")

df %>%
    filter(baseMean.gzcp > 0) %>%
    arrange(sig.gzcp) %>%
    ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change",
         title="GZ/CP Differential Expression",
         color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_label(data=subset(df, Name == "hsa-miR-124-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-124-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-92b-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-92b-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-361-5p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-361-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-186-5p"), aes(label=Name), hjust=0, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-186-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-26a-5p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-26a-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-191-5p"), aes(label=Name), hjust=0, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-191-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-423-5p"), aes(label=Name), hjust=1, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-423-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-451a"), aes(label=Name), hjust=0, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-451a"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-320a-3p"), aes(label=Name), hjust=1, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-320a-3p"), shape=21, size=4, stroke=2, color="black")

pdf("~/Desktop/miRNA_endog_control_candidates.pdf", width = 10, height = 10)

df %>%
    filter(baseMean.gw > 0) %>%
    arrange(sig.gw) %>%
    ggplot(aes(x=baseMean.gw, y=log2FoldChange.gw, color=sig.gw)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change",
         title="Gestation Week Differential Expression",
         color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_label(data=subset(df, Name == "hsa-miR-124-3p"), aes(label=Name), hjust=0, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-124-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-92b-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-92b-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-361-5p"), aes(label=Name), hjust=1, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-361-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-197-3p"), aes(label=Name), hjust=0, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-197-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-324-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-324-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-4707-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-4707-3p"), shape=21, size=4, stroke=2, color="black")

df %>%
    filter(baseMean.gzcp > 0) %>%
    arrange(sig.gzcp) %>%
    ggplot(aes(x=baseMean.gzcp, y=log2FoldChange.gzcp, color=sig.gzcp)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Shrunken Log2 Fold Change",
         title="GZ/CP Differential Expression",
         color=paste("p.adj >", padjusted.significance.threshold, sep="")) +
    scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_label(data=subset(df, Name == "hsa-miR-124-3p"), aes(label=Name), hjust=0, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-124-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-92b-3p"), aes(label=Name), hjust=1, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-92b-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-361-5p"), aes(label=Name), hjust=1, vjust=1.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-361-5p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-197-3p"), aes(label=Name), hjust=0, vjust=-0.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-197-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-324-3p"), aes(label=Name), hjust=0, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-324-3p"), shape=21, size=4, stroke=2, color="black") +
    geom_label(data=subset(df, Name == "hsa-miR-4707-3p"), aes(label=Name), hjust=1, vjust=-.5, color="black") +
    geom_point(data=subset(df, Name == "hsa-miR-4707-3p"), shape=21, size=4, stroke=2, color="black")

dev.off()

