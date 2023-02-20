# compare mirbase v22 to mirgenedb v2 annotations

library(here)
library(dplyr)
library(magrittr)
library(DESeq2)
library(ggplot2)

source(here("src/utils/lafferty_utils.R"))
#pdf("~/Desktop/mirbase_v_mirgenedb.pdf", height = 5, width = 7)

# Load SE Datasets ################################################################################
# load mirbase v22 data
se <- readRDS(here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds"))
# build DESeq data set
dds.mirbase <- DESeqDataSet(se, design = ~1)
# remove rows with only zero counts
dds.mirbase <- dds.mirbase[rowSums(counts(dds.mirbase)) > 1, ]
# estimate size factors for normalization
dds.mirbase <- estimateSizeFactors(dds.mirbase)
# variance stabilizing tranformed data
vsd.mirbase <- varianceStabilizingTransformation(dds.mirbase)
rm(se)

# load mirgenedb v2 data
se <- readRDS(here("results/rdata_files/20180906_mirgenedb_v2_mirge2.0_counts.rds"))
# build DESeq data set
dds.mirgenedb <- DESeqDataSet(se, design = ~1)
# remove rows with only zero counts
dds.mirgenedb <- dds.mirgenedb[rowSums(counts(dds.mirgenedb)) > 1, ]
# estimate size factors for normalization
dds.mirgenedb <- estimateSizeFactors(dds.mirgenedb)
# variance stabilizing tranformed data
vsd.mirgenedb <- varianceStabilizingTransformation(dds.mirgenedb)
rm(se)

# PCA analysis ####################################################################################
# mirbase v22
pca.mirbase <- prcomp(t(assay(vsd.mirbase)), center = TRUE, scale. = FALSE)
df.pca.mirbase <- data.frame(pca.mirbase$x[,1:2], colData(vsd.mirbase))
percentVar.pca.mirbase <- pca.mirbase$sdev^2 / sum(pca.mirbase$sdev^2) * 100
# mirgenedb v2
pca.mirgenedb <- prcomp(t(assay(vsd.mirgenedb)), center = TRUE, scale. = FALSE)
df.pca.mirgenedb <- data.frame(pca.mirgenedb$x[,1:2], colData(vsd.mirgenedb))
percentVar.pca.mirgenedb <- pca.mirgenedb$sdev^2 / sum(pca.mirgenedb$sdev^2) * 100

# Plot PCA ########################################################################################
df.pca.mirbase %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar.pca.mirbase[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.pca.mirbase[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (miRBase v22)",
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(3,1,2)])

df.pca.mirgenedb %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar.pca.mirgenedb[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.pca.mirgenedb[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (MirGeneDB v2)",
       caption="VST normalized, uncorrected for batch") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(3,1,2)])


# Batch Correction ################################################################################
vsd.mirbase.bc <- vsd.mirbase
assay(vsd.mirbase.bc) <- limma::removeBatchEffect(assay(vsd.mirbase),
                                                  batch = vsd.mirbase$sequencing_pool,
                                                  batch2 = vsd.mirbase$rna_extraction_date,
                                                  design = model.matrix(~vsd.mirbase$gestation_week))

vsd.mirgenedb.bc <- vsd.mirgenedb
assay(vsd.mirgenedb.bc) <- limma::removeBatchEffect(assay(vsd.mirgenedb),
                                                    batch = vsd.mirgenedb$sequencing_pool,
                                                    batch2 = vsd.mirgenedb$rna_extraction_date,
                                                    design = model.matrix(~vsd.mirgenedb$gestation_week))

# PCA analysis (BC) ###############################################################################
# mirbase v22
pca.mirbase.bc <- prcomp(t(assay(vsd.mirbase.bc)), center = TRUE, scale. = FALSE)
df.pca.mirbase.bc <- data.frame(pca.mirbase.bc$x[,1:2], colData(vsd.mirbase.bc))
percentVar.pca.mirbase.bc <- pca.mirbase.bc$sdev^2 / sum(pca.mirbase.bc$sdev^2) * 100
# mirgenedb v2
pca.mirgenedb.bc <- prcomp(t(assay(vsd.mirgenedb.bc)), center = TRUE, scale. = FALSE)
df.pca.mirgenedb.bc <- data.frame(pca.mirgenedb.bc$x[,1:2], colData(vsd.mirgenedb.bc))
percentVar.pca.mirgenedb.bc <- pca.mirgenedb.bc$sdev^2 / sum(pca.mirgenedb.bc$sdev^2) * 100

# Plot PCA (BC) ###################################################################################
df.pca.mirbase.bc %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar.pca.mirbase.bc[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.pca.mirbase.bc[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (miRBase v22)",
       caption="VST normalized, batch correct for seq. pool and rna extraction date") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(3,1,2)])

df.pca.mirgenedb.bc %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar.pca.mirgenedb.bc[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar.pca.mirgenedb.bc[2], 1), "%)", sep=""),
       title="PCA: miRNA Expression (MirGeneDB v2)",
       caption="VST normalized, batch correct for seq. pool and rna extraction date") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(3,1,2)])

# combined analysis (to see relative movement of points in PCA plot)
df <- bind_rows(df.pca.mirbase.bc, df.pca.mirgenedb.bc)
df %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point() +
  geom_line(aes(group=rnaid)) +
  labs(x="PC1",
       y="PC2",
       caption="VST normalized, batch corrected, combined mirbase and mirgenedb") +
  plotTheme +
  scale_color_manual(values = cbPalette[c(3,1,2)])

# Subset GZ/CP from CW Samples #########################################################################
# NOT IMPLIMENTED
# dds.gzcp.mirbase <- dds.mirbase[,dds.mirbase$tissue_section == "GZ" | dds.mirbase$tissue_section == "CP"]
# dds.gzcp.mirbase$tissue_section <- droplevels(dds.gzcp.mirbase$tissue_section)
# dds.gzcp.mirbase$sequencing_pool <- droplevels(dds.gzcp.mirbase$sequencing_pool)
# dds.gzcp.mirbase$rna_extraction_date <- droplevels(dds.gzcp.mirbase$rna_extraction_date)
# 
# dds.cw.mirbase <- dds.mirbase[,dds.mirbase$tissue_section == "CW"]
# dds.cw.mirbase$tissue_section <- droplevels(dds.cw.mirbase$tissue_section)
# dds.cw.mirbase$sequencing_pool <- droplevels(dds.cw.mirbase$sequencing_pool)
# dds.cw.mirbase$rna_extraction_date <- droplevels(dds.cw.mirbase$rna_extraction_date)
# 
# dds.gzcp.mirgenedb <- dds.mirgenedb[,dds.mirgenedb$tissue_section == "GZ" | dds.mirgenedb$tissue_section == "CP"]
# dds.gzcp.mirgenedb$tissue_section <- droplevels(dds.gzcp.mirgenedb$tissue_section)
# dds.gzcp.mirgenedb$sequencing_pool <- droplevels(dds.gzcp.mirgenedb$sequencing_pool)
# dds.gzcp.mirgenedb$rna_extraction_date <- droplevels(dds.gzcp.mirgenedb$rna_extraction_date)
# 
# dds.cw.mirgenedb <- dds.mirgenedb[,dds.mirgenedb$tissue_section == "CW"]
# dds.cw.mirgenedb$tissue_section <- droplevels(dds.cw.mirgenedb$tissue_section)
# dds.cw.mirgenedb$sequencing_pool <- droplevels(dds.cw.mirgenedb$sequencing_pool)
# dds.cw.mirgenedb$rna_extraction_date <- droplevels(dds.cw.mirgenedb$rna_extraction_date)


#dev.off()




