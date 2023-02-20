# miRge quantified miRNA expression, fetal brain tissue samples

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
# dds of sex diff expression
dds_rds <- here("doc/eda/rdata/dds_rin.rds")
# output directory for pdf files
pdf_dir <- here("doc/eda/pdfs/")
# write plots
write_plots <- TRUE

# INPUT FILES #####################################################################################
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# male rnaids that are not labeled by genotype data
male_rnaids_file <- here("data/metadata/male_unlabeled_rnaids.txt")
# female rnaids that are not labeled by gneotype data
female_rnaids_file <- here("data/metadata/female_unlabeled_rnaids.txt")

# Build DDS #######################################################################################
# read in summarized experiment
se <- readRDS(mirbase_se_rds)

# label outliers by PCA (after batch effect removal of pool and rna purif. method)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")
colData(se)$outlier <- ifelse(colData(se)$rnaid %in% outlier_rnaid, "yes", "no")

# read male and female rnaids
males <- read_lines(male_rnaids_file)
females <- read_lines(female_rnaids_file)
# label unlabeled males and females
colData(se)$sex[which(colData(se)$rnaid %in% males)] <- "Male"
colData(se)$sex[which(colData(se)$rnaid %in% females)] <- "Female"

# build DESeq data set
dds <- DESeqDataSet(se, design = ~1)
# remove rows with only zero counts
dds <- dds[rowSums(counts(dds)) > 1, ]
# estimate size factors for normalization
dds <- estimateSizeFactors(dds)
# norm. transformed counts (DESeq size factors, log2 transformed)
#ntc <- assay(normTransform(dds))
# variance stabilizing tranformed data
vsd <- varianceStabilizingTransformation(dds)
rm(se)

# Batch Correct for Pool and RNA Purification Method ##############################################
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             batch2 = vsd$rna_purification_method,
                                             design = model.matrix(~vsd$gestation_week))

# Plot miRNA Expression PCA vs Rin ################################################################
# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = FALSE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
  ggplot(aes(PC1, PC2, color=rin)) +
  geom_point(size=2) +
  labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
       y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
  plotTheme

p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "rin"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_text(size = 18, face = "bold"))
print(p)


# Diff. Expression ################################################################################
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

dds$sex <- factor(dds$sex)

# set design equation
design(dds) <- formula(~ sequencing_pool + rna_purification_method_name + gestation_week + rin)

#dds <- DESeq(dds)

#saveRDS(dds, dds_rds)

# load dds
dds <- readRDS(dds_rds)

# Results #########################################################################################
res <- results(dds)
# shrunken results
shrunkres <- lfcShrink(dds, coef="rin", type="apeglm")

plotMA(res)
summary(res)

plotMA(shrunkres)
summary(shrunkres)

# data frame from results
res.df <- as.data.frame(res)
res.df$mirna <- rownames(res.df)
res.df$sig <- "not_sig"
res.df$sig[which(res.df$padj < 0.1)] <- "sig"
# order data frame for plotting
res.df <- res.df[order(res.df$padj, decreasing = TRUE),]

# data frame from shrunken results
shrunkres.df <- as.data.frame(shrunkres)
shrunkres.df$mirna <- rownames(shrunkres.df)
shrunkres.df$sig <- "not_sig"
shrunkres.df$sig[which(shrunkres.df$padj < 0.1)] <- "sig"
# order data frame for plotting
shrunkres.df <- shrunkres.df[order(shrunkres.df$padj, decreasing = TRUE),]

# MA Plot #########################################################################################
circColor <- "black"
labelColor <- "black"

shrunkres.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across RIN") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  presTheme +
  scale_color_manual(values=c("grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"diff_expression_rin_shrunklfc.pdf"), width = 6, height = 6, useDingbats = FALSE)}

res.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Log2 Fold Change",
       title="Diff. Expression Across RIN") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))



