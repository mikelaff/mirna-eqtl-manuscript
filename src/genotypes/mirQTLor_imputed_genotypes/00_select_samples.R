# select mirQTLor samples
# outliers removed by PC1/PC2 after batch correction for sequencing pool and rna purification method

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(RMySQL)
library(limma)
library(mikelaffr)

# OUTPUT ###############################################################################################################
# list of mirQTLor samples for plink filtering
output.samples.txt <- here("results/genotypes/mirQTLor_imputed_genotypes/mirQTLor_samples.txt")

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

if (FALSE) {
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

    # Import Gene Expression Data ##########################################################################################
    rse.gene <- readRDS(rse.gene.rds)

    # threshold for expression
    rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]

    # build DESeq data set
    dds.gene <- DESeqDataSet(rse.gene, design = ~1)

    # vst normalization
    vsd.gene <- vst(dds.gene)

    # filter for samples within mirna dataset
    vsd.gene <- vsd.gene[,sample.metadata$RNAID]

    rm(rse.gene, dds.gene)

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

    # Expression PCA miRBase Only Outliers Removed #########################################################################
    vsd.trans.mirbase.or <- vsd.mirbase[,!vsd.mirbase$rnaid %in% outlier_rnaid]
    assay(vsd.trans.mirbase.or) <- limma::removeBatchEffect(assay(vsd.trans.mirbase.or),
                                                            batch = vsd.trans.mirbase.or$sequencing_pool,
                                                            batch2 = vsd.trans.mirbase.or$rna_purification_method,
                                                            design = model.matrix(~vsd.trans.mirbase.or$gestation_week))

    # PCA analysis (scale = FALSE)
    pca <- prcomp(t(assay(vsd.trans.mirbase.or)), center = TRUE, scale. = TRUE)
    df <- data.frame(pca$x[,1:5], colData(vsd.trans.mirbase.or))
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
    vsd.trans.or <- vsd[,!vsd$rnaid %in% outlier_rnaid]
    assay(vsd.trans.or) <- limma::removeBatchEffect(assay(vsd.trans.or),
                                                    batch = vsd.trans.or$sequencing_pool,
                                                    batch2 = vsd.trans.or$rna_purification_method,
                                                    design = model.matrix(~vsd.trans.or$gestation_week))

    # PCA analysis (scale = FALSE)
    pca <- prcomp(t(assay(vsd.trans.or)), center = TRUE, scale. = TRUE)
    df <- data.frame(pca$x[,1:5], colData(vsd.trans.or))
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

}
# Write Samples to File ################################################################################################

#--keep accepts a space/tab-delimited text file with family IDs (DonorID) in the first column and within-family
# IDs (DNAID) in the second column, and removes all unlisted samples from the current analysis. --remove does the
# same for all listed samples.

df.small.samples %>%
    filter(!RNAID %in% outlier_rnaid) -> df.output

write_lines(paste(df.output$DonorID, df.output$DNAID), output.samples.txt)


