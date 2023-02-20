# miRge quantified miRNA expression, fetal brain tissue samples

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
# summarized experiment with novel
se_rds <- here("results/mirdeep2/novel_try2.rds")

# Build DDS #######################################################################################
# read in summarized experiment
se <- readRDS(se_rds)

# label outliers by PCA (after batch effect removal of pool and rna purif. method, see below)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")
colData(se)$outlier <- ifelse(colData(se)$rnaid %in% outlier_rnaid, "yes", "no")

# build DESeq data set
dds <- DESeqDataSet(se, design = ~1)
dds <- collapseReplicates(dds, groupby = dds$rnaid)
# remove rows with only zero counts
dds <- dds[rowSums(counts(dds)) > 100, ]
# estimate size factors for normalization
dds <- estimateSizeFactors(dds)
# norm. transformed counts (DESeq size factors, log2 transformed)
#ntc <- assay(normTransform(dds))
# variance stabilizing tranformed data
vsd <- varianceStabilizingTransformation(dds)
rm(se)

# Sequencing and library depth ####################################################################
df <- as.data.frame(colData(dds))

ggplot(df, aes(x=reorder(rnaid, reads_per_sample), y=reads_per_sample/1e6, fill=tissue_section)) +
  geom_bar(stat = "identity", position = position_dodge(0),width = 1) +
  labs(x="240 Samples",
       y="Reads per Sample (million reads)",
       #title="Sequencing Depth per Sample",
       fill="Tissue\nSection") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 90)) +
  presTheme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values=c("grey70", "blue", "darkorange"))

if(write_plots){ggsave(file.path(pdf_dir, "small_rna_seq_reads_per_sample.pdf"), height = 6, width = 7)}

tmp <- data.frame(rnaid=colData(dds)$rnaid, tot_counts=colSums(counts(dds)))
df <- left_join(df, tmp, by = "rnaid")

ggplot(df, aes(x=reorder(rnaid, tot_counts), y=tot_counts/1e6, fill=tissue_section)) +
  geom_bar(stat = "identity", position = position_dodge(0),width = 1) +
  labs(x="240 Samples",
       y="miRNA Counts per Sample (million counts)",
       #title="miRNA Counts per Sample",
       fill="Tissue\nSection") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 5.5)) +
  plotTheme +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_fill_manual(values=c("grey50", "blue", "darkorange"))

ggplot(df, aes(x=reads_per_sample/1e6, y=tot_counts/1e6, col=tissue_section)) +
  geom_point() +
  labs(x="Reads per Sample (million reads)",
       y="miRNA Counts per Sample (million counts)",
       col="Tissue\nSection") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue", "darkorange"))

ggplot(df, aes(y=`sizeFactor`, x=reads_per_sample/1e6, col=tissue_section)) +
  geom_point() +
  labs(x="Reads per Sample (million reads)",
       y="DESeq2 Size Factor",
       col="Tissue\nSection") +
  plotTheme +
  scale_color_manual(values=c("grey50", "blue", "darkorange"))

# PCA #############################################################################################
pca <- prcomp(t(assay(vsd)))
df <- data.frame(pca$x[,1:2], colData(vsd))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

#rm(pca)

ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
  geom_point() +
  geom_line() +
  labs(x="PC", y="Percent Total Variance",
       title="PCA: miRNA Expression (VST Norm.)") +
  presTheme

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_scree_plot.pdf"), width = 11, height = 6, useDingbats = FALSE)}

ggplot(df, aes(PC1, PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Tissue\nSection") +
  presTheme +
  scale_color_manual(values=c("gray70", "navy", "darkorange")) +
  labs(title="PCA: miRNA Expression (VST Norm.)")

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples.pdf"), width = 11, height = 6, useDingbats = FALSE)}

ggplot(df, aes(PC1, PC2, color=sequencing_pool)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool",
       title="PCA: miRNA Expression (VST Norm.)") +
  presTheme +
  scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_by_pool.pdf"), width = 11, height = 6, useDingbats = FALSE)}

ggplot(df, aes(PC1, PC2, color=sequencing_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool") +
  plotTheme

#dev.off()

ggplot(df, aes(PC1, PC2, color=rna_extraction_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

ggplot(df, aes(y=rna_concentration, x=reorder(rnaid, rna_concentration), fill=log10(rna_concentration))) +
  geom_bar(stat="Identity", position = position_dodge(0),width = 1) +
  plotTheme +
  scale_fill_gradientn(colors = c("blue", "grey50", "grey50", "red")) +
  labs(title="RNA Concentration Across Fetal Tissue Samples",
       x="240 Samples",
       y="RNA Concentration (ng/mL)",
       fill="log10(RNA\nConc.)") +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, 2500))

ggplot(df, aes(PC1, PC2, color=log10(rna_concentration))) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)") +
  scale_color_gradientn(colors = c("blue", "grey50", "grey50", "red"))

#####################################################################
# remove seq. pool
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             design = model.matrix(~vsd$gestation_week))
# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, color=sequencing_pool)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool",
       title="PCA: miRNA Expression (VST Norm.)",
       caption="limma: remove batch effect for sequencing pool") +
  presTheme +
  scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_removed_by_pool.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=rna_extraction_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, color=rna_purification_method)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (VST Norm.)",
       caption="limma: remove batch effect for sequencing pool",
       color="RNA Purification\nMethod") +
  presTheme +
  theme(legend.text = element_text(size = 10)) +
  scale_color_manual(values=cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_removed_by_purification_method.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

# PCA analysis (scale = TRUE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, color=sequencing_pool)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool") +
  plotTheme +
  scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

df %>%
  ggplot(aes(PC1, PC2, color=rna_extraction_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

####################################################################
# remove seq. pool and rna purification method
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             batch2 = vsd$rna_purification_method,
                                             design = model.matrix(~vsd$gestation_week))

# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, color=sequencing_pool)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool") +
  plotTheme +
  scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

df %>%
  ggplot(aes(PC1, PC2, color=rna_extraction_date)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, color=rna_purification_method)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

# PCA analysis (scale = TRUE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, color=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool") +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, color=rna_purification_method)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

df %>%
  ggplot(aes(PC1, PC2, color=sequencing_pool)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Sequencing\nPool",
       title="PCA: miRNA Expression (VST Norm.)",
       caption="limma: remove batch effect for sequencing pool and RNA purification method") +
  presTheme +
  scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_pool.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=rna_purification_method)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (VST Norm.)",
       caption="          limma: remove batch effect for sequencing pool and RNA purification method",
       color="RNA Purification\nMethod") +
  presTheme +
  theme(legend.text = element_text(size = 10)) +
  scale_color_manual(values=cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_purification_method.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=rin)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="RIN") +
  presTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_rin.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=log10(reads_per_sample))) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="log10(Reads\nPer Sample)") +
  presTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_reads_per_sample.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=log10(rna_concentration))) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="log10(RNA\nConcentration)") +
  presTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_rna_concentration.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=sizeFactor)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Size Factor") +
  presTheme +
  scale_color_gradientn(colors = c("blue", "grey70", "grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_size_factor.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, PC2, color=gestation_week)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (VST Norm.)",
       caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Gestation\nWeek") +
  presTheme +
  scale_color_gradientn(colors = c("darkorange", "grey70", "grey70", "blue"))

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_gestation_week.pdf"), width = 11, height = 6, useDingbats = FALSE)}

##############################################################

df <- data.frame(df, let7b3p_exp = as.vector(assay(vsd.trans["hsa-let-7b-3p",])))

df %>%
  ggplot(aes(PC1, PC2, color=outlier)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (VST Norm.)",
       caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Outlier") +
  presTheme +
  scale_color_manual(values = cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"pca_all_samples_pool_and_pur_method_removed_by_outlier.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC1, let7b3p_exp, color=outlier)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y="let-7b-3p Expression (VST)",
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Outlier") +
  presTheme +
  scale_color_manual(values = cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"pca1_v_let7b3p_expression_by_outler.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df %>%
  ggplot(aes(PC2, let7b3p_exp, color=outlier)) +
  geom_point(size=2) +
  labs(x=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       y="let-7b-3p Expression (VST)",
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Outlier") +
  presTheme +
  scale_color_manual(values = cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"pca2_v_let7b3p_expression_by_outler.pdf"), width = 11, height = 6, useDingbats = FALSE)}

df <- df[order(df$outlier),]
df %>%
  ggplot(aes(gestation_week, let7b3p_exp, color=outlier)) +
  geom_point(size=2) +
  labs(x="Gestation Week",
       y="let-7b-3p Expression (VST)",
       #title="PCA: miRNA Expression (VST Norm.)",
       #caption="limma: remove batch effect for sequencing pool and RNA purification method",
       color="Outlier") +
  presTheme +
  scale_color_manual(values = cbPalette)

if(write_plots){ggsave(file.path(pdf_dir,"let7b3p_expression_v_gestation_week_by_outler.pdf"), width = 11, height = 6, useDingbats = FALSE)}

##################################################################

#pdf("~/Desktop/mirge_all_pca.pdf", height = 6, width = 8)

ggplot(df, aes(PC1, PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  scale_color_manual(values=c("gray50", "blue", "darkorange")) +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date")

#dev.off()
ggplot(df, aes(PC1, PC2, color=sequencing_pool)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date")

ggplot(df, aes(PC1, PC2, color=log10(rna_concentration))) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date",
       color="log10(RNA\nConc.)") +
  scale_color_gradientn(colors = c("blue", "grey50", "grey50", "red"))

ggplot(df, aes(PC1, PC2, color=rin)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date") +
  scale_color_gradientn(colors = c("blue", "grey50", "grey50", "red"))

#pdf("~/Desktop/pca_gest_week.pdf", height = 5, width = 7)

ggplot(df, aes(PC1, PC2, color=gestation_week)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       color="Gestation\nWeek") +
  plotTheme +
  scale_color_gradientn(colors = c("darkorange", "grey80", "blue"))

#dev.off()

ggplot(df, aes(PC1, PC2, color=log10(reads_per_sample))) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date") +
  scale_color_gradientn(colors = c("blue", "grey50", "grey50", "red"))

ggplot(df, aes(PC1, PC2, color=rna_extraction_date)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool")

df %>%
  ggplot(aes(PC1, PC2, color=sizeFactor)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme +
  labs(title="Mirge Expression (VST Norm. Collapsed Counts)",
       caption="limma removeBatchEffect() by sequencing_pool and rna_extraction_date") +
  scale_color_gradientn(colors = c("blue", "grey50", "grey50", "red"))



p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "mds_c1", "mds_c2", "mds_c3", "mds_c4"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, face = "bold"))
print(p)

ggpairs(df,
        columns = c("PC1", "PC2", "rin", "sequencing_pool", "sizeFactor", "sex"))

df %>%
  ggpairs(columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "mds_c1", "mds_c2", "mds_c3", "mds_c4"))



