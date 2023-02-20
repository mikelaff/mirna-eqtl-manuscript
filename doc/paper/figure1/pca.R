
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(limma)
library(RMySQL)
library(mikelaffr)

# FIGURE ###############################################################################################################
# figure 1: batch corrected PCA on mirQTL samples (212), using only miRBase miRNA expression
output.pdf <- paste0(here("doc/paper/figure1/pdfs/"), "figure1_pca_by_gest_week.pdf")

# OUTPUT FILES #########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/paper/figure1/pdfs/")


# INPUT FILES ##########################################################################################################
# mirQTL samples used in EMMAX association analysis
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# SteinLab Database ###################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")
rm(con)

# Load Samples #########################################################################################################
# mirQTL samples
df.samples <- read_table2(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"))

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)
sample.metadata %<>%
    dplyr::filter(RNAID %in% df.samples$RNAID)

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
rm(donors)

# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################
sample.metadata$PurificationMethod[which(sample.metadata$PurificationMethod == "Trizol w Glycogen")] <- "*Trizol w Glycogen(colum purified)"
# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################

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
#levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "Trizol w Glycogen"] <- "trizol"

# format columns correctly
sample.metadata$Sex <- as.factor(sample.metadata$Sex)


# Load Expression Data #################################################################################################
rse <- readRDS(rse.rds)

# use only miRBase expression
rse <- rse[mcols(rse)$source == "miRBase_v22",]
# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd <- varianceStabilizingTransformation(dds)

# subset samples
vsd <- vsd[,df.samples$RNAID]

rm(rse, dds)

# Batch Effect Correction ##############################################################################################
# remove seq. pool and rna purification method
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.trans),
                                             batch = vsd.trans$sequencing_pool,
                                             batch2 = vsd.trans$rna_purification_method,
                                             design = model.matrix(~vsd.trans$gestation_week))


# Expression PCA #######################################################################################################
expr.pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100

df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))


# Plot #################################################################################################################


# df.pca %>%
#     ggplot(aes(PC1, PC2, color=Pool)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Sequencing\nPool") +
#     plotTheme() +
#     scale_color_manual(values = cbPalette, labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))
#
#
# df.pca %>%
#     ggplot(aes(PC1, PC2, color=PurificationMethod)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
#     plotTheme()
#
# df.pca %>%
#     ggplot(aes(PC1, PC2, color=Sex)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
#     plotTheme()
#
# df.pca %>%
#     ggplot(aes(PC1, PC2, color=RIN)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
#     plotTheme()

df.pca %>%
    ggplot(aes(PC1, PC2, color=GestationWeek)) +
    geom_point(size=1) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme(theme = "figure") +
    theme(legend.position = "none") +
    scale_color_gradientn(colors = paperRedToBlue)

ggsave(output.pdf, width = 2, height = 2, useDingbats = FALSE)

stop()
# # Expression PCA #######################################################################################################
# expr.pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
# df.expr.pca <- data.frame(expr.pca$x[,1:10])
# df.expr.pca$RNAID <- rownames(df.expr.pca)
#
# percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100
#
# df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))
#
#
# # Plot #################################################################################################################
#
#
# df.pca %>%
#     ggplot(aes(PC1, PC2, color=Pool)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Seq.\nPool") +
#     plotTheme(theme = "presentation") +
#     scale_color_manual(values = paperCatEight[c(2,1,3,4,5,6,7,8)], labels = c("Pool 1", "Pool 2", "Pool 3", "Pool 4", "Pool 5", "Pool 6", "Pool 7", "Pool 8"))
#
# ggsave("~/Desktop/pca_pool.pdf", width = 6, height = 5, useDingbats = FALSE)
#
# vsd.trans <- vsd
# assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.trans),
#                                              batch = vsd.trans$sequencing_pool,
#                                              design = model.matrix(~vsd.trans$gestation_week))
#
# expr.pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
# df.expr.pca <- data.frame(expr.pca$x[,1:10])
# df.expr.pca$RNAID <- rownames(df.expr.pca)
#
# percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100
#
# df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))
#
# df.pca %>%
#     ggplot(aes(PC1, PC2, color=PurificationMethod)) +
#     geom_point(size=2) +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
#          color="Purif.\nMethod") +
#     plotTheme(theme = "presentation") +
#     scale_color_manual(values = paperCatEight[c(2,1,3,4,5,6,7,8)], labels = c("Trizol+Col","miRNeasy","miRNeasy_mini"))
#
# ggsave("~/Desktop/pca_purif.pdf", width = 7, height = 5, useDingbats = FALSE)
#
#



# Suppliments ##############
# mirQTL samples
#df.samples <- read_table2(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"))

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)
#sample.metadata %<>%
#    dplyr::filter(RNAID %in% df.samples$RNAID)

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
#sample.metadata %<>%
#    left_join(select(donors, DonorID, GestationWeek), by = "DonorID")
#rm(donors)

# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################
sample.metadata$PurificationMethod[which(sample.metadata$PurificationMethod == "Trizol w Glycogen")] <- "*Trizol w Glycogen(colum purified)"
# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################

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
#levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "Trizol w Glycogen"] <- "trizol"

# format columns correctly
sample.metadata$Sex <- as.factor(sample.metadata$Sex)

rse <- readRDS(rse.rds)

# use only miRBase expression
rse <- rse[mcols(rse)$source == "miRBase_v22",]
# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd <- varianceStabilizingTransformation(dds)

# subset samples
vsd <- vsd[,colData(vsd)$tissue_section == "CW"]

rm(rse, dds)


expr.pca <- prcomp(t(assay(vsd)), center = TRUE, scale. = TRUE)
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100

df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))


df.pca %>%
    ggplot(aes(PC1, PC2, color=Pool)) +
    geom_point(size=0.7) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme(theme = "figure") +
    theme(legend.position = "none") +
    scale_color_manual(values = paperCatEight)

ggsave("~/Desktop/pca_pool.pdf", width = 2, height = 2, useDingbats = FALSE)


vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.trans),
                                             batch = vsd.trans$sequencing_pool,
                                             design = model.matrix(~vsd.trans$gestation_week))


expr.pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100

df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))


df.pca %>%
    ggplot(aes(PC1, PC2, color=PurificationMethod)) +
    geom_point(size=0.7) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme(theme = "figure") +
    theme(legend.position = "none") +
    scale_color_manual(values = paperCatEight)

ggsave("~/Desktop/pca_purMethod.pdf", width = 2, height = 2, useDingbats = FALSE)

vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd.trans),
                                             batch = vsd.trans$sequencing_pool,
                                             batch2 = vsd.trans$rna_purification_method,
                                             design = model.matrix(~vsd.trans$gestation_week))

expr.pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
df.expr.pca <- data.frame(expr.pca$x[,1:10])
df.expr.pca$RNAID <- rownames(df.expr.pca)

percentVar <- expr.pca$sdev^2 / sum(expr.pca$sdev^2) * 100

df.pca <- as_tibble(left_join(df.expr.pca, sample.metadata, by = "RNAID"))

# label outliers by PCA (after batch effect removal of pool and rna purif. method, see below)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")

df.pca$outlier <- FALSE
df.pca$outlier[df.pca$RNAID %in% outlier_rnaid] <- TRUE

df.pca %>%
    ggplot(aes(PC1, PC2, color=outlier)) +
    geom_point(size=0.7) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme(theme = "figure") +
    theme(legend.position = "right") +
    scale_color_manual(values = paperCatEight)

ggsave("~/Desktop/pca_outlier_legend.pdf", width = 2, height = 2, useDingbats = FALSE)

