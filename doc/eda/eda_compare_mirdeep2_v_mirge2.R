# compare mirdeep2 quantification to mirge2.0 on mirbase v22

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(pheatmap)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT ##############
pdf_dir <- here("doc/eda/pdfs/")

# INPUT ###############
# summarized experiment of mirdeep2 quantified counts
mirdeep_se_rds <- here("results/rdata_files/20190204_mirbase_v22_mirdeep2_counts_se.rds")
# summarized experiment of mirge2.0 quantified counts
mirge_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# df of mirdeep2 counts
mirdeep_df_tsv <- here("results/mirdeep2/20190203_mirbase_v22_mirdeep2_counts.tsv")
# df of mirge2.0 counts
mirge_df_tsv <- here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv")


mirdeep_df <- read_tsv(mirdeep_df_tsv)
mirge_df <- read_tsv(mirge_df_tsv)


mirge_df %<>%
  mutate(key = sapply(strsplit(miRNA, "/"), `[`, 1)) %>%
  select(-miRNA) %>%
  select(key, everything())

mirdeep_df %<>%
  mutate(key = sapply(strsplit(mirna_precursor, "_"), `[`, 1), precursor = sapply(strsplit(mirna_precursor, "_"), `[`, 2)) %>%
  select(key, precursor, everything()) %>%
  group_by(key) %>%
  summarise_at(.vars = 4:243, .funs = funs(mean))

mirdeep_df[,2:241] <- lapply(mirdeep_df[,2:241], function(x) as.integer(x))

sum(duplicated(mirge_df$key))
sum(duplicated(mirdeep_df$key))

sum(mirge_df$key %in% mirdeep_df$key)
sum(mirdeep_df$key %in% mirge_df$key)



mirge_not_in_mirdeep <- mirge_df[!mirge_df$key %in% mirdeep_df$key,]
mirdeep_not_in_mirge <- mirdeep_df[!mirdeep_df$key %in% mirge_df$key,]


m.mirge <- as.matrix(mirge_df[,2:241])
rownames(m.mirge) <- mirge_df$key

m.mirdeep <- as.matrix(mirdeep_df[,2:241])
rownames(m.mirdeep) <- mirdeep_df$key

m.mirge <- m.mirge[mirge_df$key %in% mirdeep_df$key,]
m.mirdeep <- m.mirdeep[mirdeep_df$key %in% mirge_df$key,]

m.mirge <- m.mirge[rownames(m.mirdeep), colnames(m.mirdeep)]

all(rownames(m.mirge) == rownames(m.mirdeep))
all(colnames(m.mirge) == colnames(m.mirdeep))

m.mirge.trunc <- m.mirge[rowSums(m.mirge > 10) > 10 | rowSums(m.mirdeep > 10) > 10,]
m.mirdeep.trunc <- m.mirdeep[rowSums(m.mirge > 10) > 10 | rowSums(m.mirdeep > 10) > 10,]

all(rownames(m.mirge.trunc) == rownames(m.mirdeep.trunc))
all(colnames(m.mirge.trunc) == colnames(m.mirdeep.trunc))

m.diff.perc <- 100 * abs(m.mirge.trunc - m.mirdeep.trunc) / pmax(m.mirge.trunc, m.mirdeep.trunc)
m.diff.perc[is.nan(m.diff.perc)] <- 0

mean(m.diff.perc)
# percent difference from max value, on mirnas with > 10 counts in > 10 samples
pheatmap(mat = m.diff.perc,
         cluster_rows = TRUE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         show_colnames = FALSE)

df.diff <- data.frame(mirna=rownames(m.mirge.trunc),
                      mean_mirge=rowMeans(m.mirge.trunc),
                      mean_mirdeep=rowMeans(m.mirdeep.trunc),
                      perc_diff=rowMeans(m.diff.perc))

df.diff %>%
  ggplot(aes(x=mean_mirge, y=mean_mirdeep, color=perc_diff)) +
  geom_point() +
  plotTheme

df.diff %>%
  ggplot(aes(x=log2(mean_mirge+1), y=log2(mean_mirdeep+1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1,intercept = 0) +
  presTheme +
  labs(
       y="miRDeep2 log2(mean expression + 1)",
       x="miRge2.0 log2(mean expression + 1)")

ggsave(filename = paste(pdf_dir, "mirna_mean_expression_mirge_mirdeep.pdf", sep="/"), height=6, width=6, useDingbats = FALSE)

fit <- lm(log2(mean_mirdeep+1) ~ log2(mean_mirge+1), data = df.diff)
summary(fit)

cor(df.diff$mean_mirdeep, df.diff$mean_mirge)
cor.test(log2(df.diff$mean_mirdeep+1), log2(df.diff$mean_mirge+1))


# Import SEs ###############
mirdeep_se <- readRDS(mirdeep_se_rds)
mirge_se <- readRDS(mirge_se_rds)

dds.mirge <- DESeqDataSet(mirge_se, design = ~1)
dds.mirdeep <- DESeqDataSet(mirdeep_se, design = ~1)

dds.mirge <- dds.mirge[rowSums(counts(dds.mirge) > 10) > 10,]
dds.mirdeep <- dds.mirdeep[rowSums(counts(dds.mirdeep) > 10) > 10,]


vsd.mirge <- varianceStabilizingTransformation(dds.mirge)
vsd.mirdeep <- varianceStabilizingTransformation(dds.mirdeep)

#pdf(paste(pdf_dir, "mirge_v_mirdeep_pca.pdf", sep="/"))

# PCA mirge ######
# vst norm, uncorrected
pca <- prcomp(t(assay(vsd.mirge)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.mirge))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(x=PC1, y=PC2, color=sequencing_pool)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirge, vst norm, uncorrected",
       title="mirge vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

# vst norm, seq pool corrected
vsd.trans <- vsd.mirge
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.mirge),
                                             batch = vsd.mirge$sequencing_pool)

pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(x=PC1, y=PC2, color=rna_purification_method)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirge, vst norm, seq pool corrected",
       title="mirge vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

# vst norm, seq pool and purificaiton method corrected
vsd.trans <- vsd.mirge
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.mirge),
                                             batch = vsd.mirge$sequencing_pool,
                                             batch2 = vsd.mirge$rna_purification_method)

pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
df.mirge <- df

df %>%
  ggplot(aes(x=PC1, y=PC2, color=sequencing_pool)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirge, vst norm, seq pool and purif. method corrected",
       title="mirge vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=gestation_week)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirge, vst norm, seq pool and purif. method corrected",
       title="miRge2.0") +
  presTheme +
  scale_color_gradientn(colors = c("red", "grey70", "blue"))

ggsave(filename = paste(pdf_dir, "mirge_expression_gest_week.pdf", sep="/"), height=7, width=8, useDingbats = FALSE)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirge, vst norm, seq pool and purif. method corrected",
       title="miRge2.0") +
  presTheme +
  scale_color_manual(values=c("grey70", "blue", "red"))

ggsave(filename = paste(pdf_dir, "mirge_expression_tissue_section.pdf", sep="/"), height=7, width=8, useDingbats = FALSE)


# PCA mirdeep ######
# vst norm, uncorrected
pca <- prcomp(t(assay(vsd.mirdeep)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.mirdeep))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(x=PC1, y=PC2, color=rna_purification_method)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, uncorrected",
       title="mirdeep vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

# vst norm, purificaiton method corrected
vsd.trans <- vsd.mirdeep
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.mirdeep),
                                             batch = vsd.mirdeep$rna_purification_method)

pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(x=PC1, y=PC2, color=rna_purification_method)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method corrected",
       title="mirdeep vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=sequencing_pool)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method corrected",
       title="mirdeep vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

# vst norm, purificaiton method and seq pool corrected
vsd.trans <- vsd.mirdeep
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.mirdeep),
                                             batch = vsd.mirdeep$rna_purification_method,
                                             batch2 = vsd.mirdeep$sequencing_pool)

pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100
df.mirdeep <- df

df %>%
  ggplot(aes(x=PC1, y=PC2, color=rna_purification_method)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method and seq pool corrected",
       title="mirdeep vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=sequencing_pool)) +
  geom_point() +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method and seq pool corrected",
       title="mirdeep vst") +
  plotTheme +
  scale_color_manual(values=cbPalette)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=gestation_week)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method and seq pool corrected",
       title="miRDeep2") +
  presTheme +
  scale_color_gradientn(colors = c("red", "grey70", "blue"))

ggsave(filename = paste(pdf_dir, "mirdeep_expression_gest_week.pdf", sep="/"), height=7, width=8, useDingbats = FALSE)

df %>%
  ggplot(aes(x=PC1, y=PC2, color=tissue_section)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
       caption="mirdeep, vst norm, purif. method and seq pool corrected",
       title="miRDeep2") +
  presTheme +
  scale_color_manual(values=c("grey70", "blue", "red"))

ggsave(filename = paste(pdf_dir, "mirdeep_expression_tissue_section.pdf", sep="/"), height=7, width=8, useDingbats = FALSE)

df.mirge$method <- "miRge2.0"
df.mirdeep$method <- "miRDeep2"

df.combo <- bind_rows(df.mirge, df.mirdeep)

df.combo %>%
  ggplot(aes(x=PC1, y=PC2, color=method)) +
  geom_line(aes(group=rnaid), color="black") +
  geom_point(size=2) +
  presTheme +
  scale_color_manual(values=c("red", "blue"))

ggsave(filename = paste(pdf_dir, "mirge_mirdeep_expression.pdf", sep="/"), height=7, width=8, useDingbats = FALSE)


dev.off()







