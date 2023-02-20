
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(RMySQL)
library(limma)
library(pheatmap)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdf <- here("doc/mirQTL_outlier_removal/pdfs/")

# INPUT ################################################################################################################
# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment for known and novel mirna expression values
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# ranged summarized experiment for mrna expression values
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# this should be a database query, however, here is the table of all small rna-seq samples
samples.small.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")
samples.total.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_totalRNAseq_metadata.tsv")

# plink .fam file to get sexes by genotype
samples.fam <- here("results/genotypes/AllSamplesQC.fam")

# GLOBALS ##############################################################################################################
# label outliers by PCA (after batch effect removal of pool and rna purif. method, see below)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")

# Select Samples #######################################################################################################
# Select samples for fetal tissue cortical wall small-rna-seq eQTL analysis

# load small-rna-seq metadata
df.small.samples <- read_tsv(samples.small.tsv)

# subset for cortical wall samples
df.small.samples %<>%
    filter(TissueSection == "CW")

# load total-rna-seq metadata
df.total.samples <- read_tsv(samples.total.tsv)

# subset for cortical wall samples
df.total.samples %<>%
    filter(TissueSection == "CW")

# import .fam file
df.fam <- read.table(samples.fam, col.names = c("DonorID", "DNAID", "ID3", "ID4", "Sex", "Pheno"), stringsAsFactors = FALSE)
df.fam$index <- 1:length(df.fam$DonorID)

# convert sex and select columns
df.sex <- select(df.fam, DNAID, Sex, index)
df.sex$Sex <- ifelse(df.fam$Sex == 1, "M", ifelse(df.fam$Sex == 2, "F", NA))
df.sex %<>%
    select(DNAID, Sex.by.Genotype = Sex, index)

# tmp <- left_join(df.small.samples, df.total.samples, by = "RNAID", suffix = c(".small", ".total"))
# tmp %<>%
#     select(RNAID, starts_with("DonorID"), starts_with("VerifiedDNAID"))

# filter for samples with a VarifiedDNAID and no DonorID mismatches
df.small.samples %<>%
    filter(VerifiedDNAID != "MISSING" & VerifiedDNAID != "MIXTURE",
           DonorID.MISMATCH == FALSE)

df.total.samples %<>%
    filter(VerifiedDNAID != "MISSING" & VerifiedDNAID != "MIXTURE",
           DonorID.MISMATCH == FALSE)

# filter for small rna-seq samples in total rna-seq samples
df.small.samples %<>%
    filter(RNAID %in% df.total.samples$RNAID)

# join with genotype sexes
df.small.samples %<>%
    left_join(df.sex, by = c("VerifiedDNAID" = "DNAID"))

# label sex mismaches
df.small.samples$Sex.MISMATCH <- ! df.small.samples$Sex.by.XIST == df.small.samples$Sex.by.Genotype

# remove sex mismatches
df.small.samples %<>%
    filter(!Sex.MISMATCH)

# arrange by DonorID (this is how .fam files are organized)
df.small.samples %<>%
    arrange(index)

stopifnot(sum(duplicated(df.small.samples$RNAID)) == 0,
          sum(duplicated(df.small.samples$VerifiedDNAID)) == 0,
          sum(duplicated(df.small.samples$DonorID)) == 0)

rm(df.fam, df.sex, df.total.samples)

# only retain needed sample information
df.small.samples %<>%
    select(RNAID, DonorID, DNAID = VerifiedDNAID)

# Get Batch Data #######################################################################################################
# SteinLab Database ####################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")

# Import Files #########################################################################################################
# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)

# Format Sample Data ###################################################################################################
# Association design
# Expression ~ Genotype + Kinship + 10 Genotype MDS or PCA + Seq Pool + Purif Method + RIN + Sex + Gest Week

# filter samples
sample.metadata %<>%
    filter(RNAID %in% df.small.samples$RNAID)

# select only needed information
sample.metadata %<>%
    select(RNAID,
           DonorID,
           DNAID = VerifiedDNAID,
           Pool,
           PurificationMethod,
           RIN,
           Sex = Sex.by.XIST)

# join with donors for gestation week
sample.metadata %<>%
    left_join(select(donors, DonorID, GestationWeek), by = "DonorID")

# relabel sequencing_pool factor labels
sample.metadata$Pool <- factor(sample.metadata$Pool)
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool1"] <- "Pool1"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool2"] <- "Pool2"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool3"] <- "Pool3"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool4"] <- "Pool4"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool5"] <- "Pool5"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool6"] <- "Pool6"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool7"] <- "Pool7"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool8"] <- "Pool8"

# relabel rna_purification_method factor labels
sample.metadata$PurificationMethod <- factor(sample.metadata$PurificationMethod)
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "Trizol w Glycogen"] <- "trizol"

# format columns correctly
sample.metadata$Sex <- as.factor(sample.metadata$Sex)

# add outlier information
sample.metadata$Outlier <- sample.metadata$RNAID %in% outlier_rnaid

rm(con, donors)

# Import miRNA Expression Data #########################################################################################
rse <- readRDS(rse.mirna.rds)

rse.mirbase <- rse[mcols(rse)$source == "miRBase_v22",]

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]
rse.mirbase <- rse.mirbase[rowSums(assay(rse.mirbase) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)
dds.mirbase <- DESeqDataSet(rse.mirbase, design = ~1)

vsd <- varianceStabilizingTransformation(dds)
vsd.mirbase <- varianceStabilizingTransformation(dds.mirbase)

# subset samples
vsd <- vsd[,sample.metadata$RNAID]
vsd.mirbase <- vsd.mirbase[,sample.metadata$RNAID]

rm(rse, dds, rse.mirbase, dds.mirbase)

# Expression PCA  miRBase only #########################################################################################
vsd.trans.mirbase <- vsd.mirbase
assay(vsd.trans.mirbase) <- limma::removeBatchEffect(assay(vsd.mirbase),
                                                     batch = vsd.mirbase$sequencing_pool,
                                                     batch2 = vsd.mirbase$rna_purification_method,
                                                     design = model.matrix(~vsd.mirbase$gestation_week))

# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans.mirbase)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans.mirbase))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
    ggplot(aes(PC1, PC2, color=sequencing_pool)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Sequencing\nPool") +
    plotTheme() +
    scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

df %>%
    ggplot(aes(PC1, PC2, color=gestation_week)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Gestation\nWeek",
         title="miRBase Only",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_gradientn(colors = c("blue", "grey80", "red"))

df %>%
    ggplot(aes(PC1, PC2, color=Outlier)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Outlier",
         title="miRBase Only",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_manual(values = cbPalette)

# Expression PCA Known and Novel #######################################################################################
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             batch2 = vsd$rna_purification_method,
                                             design = model.matrix(~vsd$gestation_week))

# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
    ggplot(aes(PC1, PC2, color=sequencing_pool)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Sequencing\nPool") +
    plotTheme() +
    scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))

df %>%
    ggplot(aes(PC1, PC2, color=gestation_week)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Gestation\nWeek",
         title="Known and Novel",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_gradientn(colors = c("blue", "grey80", "red"))

df %>%
    ggplot(aes(PC1, PC2, color=Outlier)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Outlier",
         title="Known and Novel",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_manual(values = cbPalette)

# PCA Before Batch Correction ####################
pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
    ggplot(aes(PC1, PC2, color=Outlier)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Outlier",
         title="Known and Novel",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_manual(values = cbPalette)

source(here("src/utils/lafferty_utils.R"))

#pdf(paste0(dir.pdf, "PC1-10_by_PEER1-10.pdf"), height = 10, width = 10)
p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)
#dev.off()

# Outliers removed
vsd.or <- vsd[,!vsd$rnaid %in% outlier_rnaid]

pca <- prcomp(t(assay(vsd.or)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd.or))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

df %>%
    ggplot(aes(PC1, PC2, color=Outlier)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Outlier",
         title="Known and Novel",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme() +
    scale_color_manual(values = cbPalette)

p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)

pdf(paste0(dir.pdf, "PC1-5_with_batchVars.pdf"), height = 10, width = 10)
p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "gestation_week", "Sex", "RIN", "Pool"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)
dev.off()

# Heatmap PCA and Batch ####################

# Outliers removed mirbase only
vsd.mirbase.or <- vsd.mirbase[,!vsd.mirbase$rnaid %in% outlier_rnaid]

pca <- prcomp(t(assay(vsd.mirbase.or)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd.mirbase.or))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# df %>%
#     ggplot(aes(PC1, PC2)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Outlier",
#          title="miRBase Only",
#          caption="Batch Corrected: Seq. Pool and Purif. Method") +
#     plotTheme() +
#     scale_color_manual(values = cbPalette)


# outliers removed
df.or <- data.frame(pca$x[,1:20])
df.or$RNAID <- rownames(df.or)
df.or <- left_join(df.or, sample.metadata, by = "RNAID")

batchVars <- c("Pool", "PurificationMethod", "Sex", "RIN", "GestationWeek")

# create model matrix with i pca components
exprPCA <- paste0("PC", 1:20)

formula.string <- paste("~",
                        paste(paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.or))

cormat <- cor(modelmat[,-1], use = "pair")

pdf(paste0(dir.pdf, "outliers_removed_mirbase_only_PC20_heatmap.pdf"))
pheatmap(cormat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Outliers removed all known and novel
vsd.or <- vsd[,!vsd$rnaid %in% outlier_rnaid]

pca <- prcomp(t(assay(vsd.or)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd.or))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# df %>%
#     ggplot(aes(PC1, PC2)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Outlier",
#          title="miRBase Only",
#          caption="Batch Corrected: Seq. Pool and Purif. Method") +
#     plotTheme() +
#     scale_color_manual(values = cbPalette)


# outliers removed
df.or <- data.frame(pca$x[,1:20])
df.or$RNAID <- rownames(df.or)
df.or <- left_join(df.or, sample.metadata, by = "RNAID")

batchVars <- c("Pool", "PurificationMethod", "Sex", "RIN", "GestationWeek")

# create model matrix with i pca components
exprPCA <- paste0("PC", 1:20)

formula.string <- paste("~",
                        paste(paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.or))

cormat <- cor(modelmat[,-1], use = "pair")

pdf(paste0(dir.pdf, "outliers_removed_known_and_novel_PC20_heatmap.pdf"))
pheatmap(cormat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Remove extra oulier ###############
# remove extra outlier
or.plus <- "RNAID1546"

# Outliers removed mirbase only
vsd.mirbase.or <- vsd.mirbase[,!vsd.mirbase$rnaid %in% c(outlier_rnaid, or.plus)]

pca <- prcomp(t(assay(vsd.mirbase.or)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd.mirbase.or))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# df %>%
#     ggplot(aes(PC1, PC2)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Outlier",
#          title="miRBase Only",
#          caption="Batch Corrected: Seq. Pool and Purif. Method") +
#     plotTheme() +
#     scale_color_manual(values = cbPalette)


# outliers removed
df.or <- data.frame(pca$x[,1:20])
df.or$RNAID <- rownames(df.or)
df.or <- left_join(df.or, sample.metadata, by = "RNAID")
df.or$PurificationMethod <- droplevels(df.or$PurificationMethod)

batchVars <- c("Pool", "PurificationMethod", "Sex", "RIN", "GestationWeek")

# create model matrix with i pca components
exprPCA <- paste0("PC", 1:20)

formula.string <- paste("~",
                        paste(paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.or))

cormat <- cor(modelmat[,-1], use = "pair")

pdf(paste0(dir.pdf, "outliers_removedPLUS_mirbase_only_PC20_heatmap.pdf"))
pheatmap(cormat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()

# Outliers removed all known and novel
vsd.or <- vsd[,!vsd$rnaid %in% c(outlier_rnaid, or.plus)]

pca <- prcomp(t(assay(vsd.or)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:20], colData(vsd.or))
df <- left_join(df, sample.metadata, by = c("rnaid" = "RNAID"))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# df %>%
#     ggplot(aes(PC1, PC2)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Outlier",
#          title="miRBase Only",
#          caption="Batch Corrected: Seq. Pool and Purif. Method") +
#     plotTheme() +
#     scale_color_manual(values = cbPalette)


# outliers removed
df.or <- data.frame(pca$x[,1:20])
df.or$RNAID <- rownames(df.or)
df.or <- left_join(df.or, sample.metadata, by = "RNAID")
df.or$PurificationMethod <- droplevels(df.or$PurificationMethod)

batchVars <- c("Pool", "PurificationMethod", "Sex", "RIN", "GestationWeek")

# create model matrix with i pca components
exprPCA <- paste0("PC", 1:20)

formula.string <- paste("~",
                        paste(paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.or))

cormat <- cor(modelmat[,-1], use = "pair")

pdf(paste0(dir.pdf, "outliers_removedPLUS_known_and_novel_PC20_heatmap.pdf"))
pheatmap(cormat,
         cluster_rows = FALSE,
         cluster_cols = FALSE)
dev.off()
