# look and the correlations between PEER factors and known batch effects and PCs

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(GGally)
library(RMySQL)
library(limma)

# OUTPUT ###############################################################################################################
dir.pdf <- here("doc/peer/pdfs/")

# INPUT ################################################################################################################
# small rna-seq samples used
small.qtl.samples.txt <- here("results/emmax/samples/20191101_mirQTL_RNAID_DonorID_DNAID.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# pre-computed PEER factors
df.peer.tsv <- here("results/peer/peer_factors_mirQTL.tsv")

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
# sample in analysis
samples <- read_delim(small.qtl.samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)

# Format Sample Data ###################################################################################################
# Association design
# Expression ~ Genotype + Kinship + 10 Genotype MDS or PCA + Seq Pool + Purif Method + RIN + Sex + Gest Week

# filter samples
sample.metadata %<>%
    filter(RNAID %in% samples$RNAID)

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

rm(con, donors, samples)

# Import Expression Data ###############################################################################################
rse <- readRDS(rse.rds)

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

# Expression PCA #######################################################################################################
pca.expr <- prcomp(t(assay(vsd)))
df.pca.expr <- data.frame(pca.expr$x[,1:20])
df.pca.expr$RNAID <- rownames(df.pca.expr)
percentVar.expr <- pca.expr$sdev^2 / sum(pca.expr$sdev^2) * 100

# PEER Factors #########################################################################################################
df.peer <- read_tsv(df.peer.tsv)
colnames(df.peer) <- c(paste0("PEER", seq(1,20,1)), "RNAID")

# Join PCA, PEER, and Metadata #########################################################################################

df <- left_join(sample.metadata, df.pca.expr, by = "RNAID")
df <- left_join(df, df.peer, by = "RNAID")

# Plot #################################################################################################################

df %>%
    ggplot(aes(x = PC1, y = PC2, color = GestationWeek)) +
    geom_point()

ggpairs(data = df,
        mapping = c("PC1", "PC2", "PC3", "PC4"))

source(here("src/utils/lafferty_utils.R"))

pdf(paste0(dir.pdf, "PC1-10_by_PEER1-10.pdf"), height = 10, width = 10)
p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10",
                         "PEER1", "PEER2", "PEER3", "PEER4", "PEER5", "PEER6", "PEER7", "PEER8", "PEER9", "PEER10"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)
dev.off()

pdf(paste0(dir.pdf, "PC1-5_by_batch.pdf"), height = 10, width = 10)
p <- ggpairs(df,
             columns = c("PC1", "PC2", "PC3", "PC4", "PC5", "Pool", "PurificationMethod", "RIN", "Sex", "GestationWeek"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)
dev.off()

pdf(paste0(dir.pdf, "PEER1-5_by_batch.pdf"), height = 10, width = 10)
p <- ggpairs(df,
             columns = c("PEER1", "PEER2", "PEER3", "PEER4", "PEER5", "Pool", "PurificationMethod", "RIN", "Sex", "GestationWeek"),
             upper = list(continuous = wrap(my_custom_cor, sizeRange = c(4,4))),
             lower = list(continuous = wrap(my_custom_smooth)),
             diag = list(continuous = wrap(my_custom_density))) +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_text(size = 10, face = "bold"))
print(p)
dev.off()

# Remove Batch Effects ##############
# remove seq. pool
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             design = model.matrix(~vsd$gestation_week))
# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans)), center = TRUE, scale. = TRUE)
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

# remove seq. pool and rna purification method
vsd.trans <- vsd
assay(vsd.trans) <- limma::removeBatchEffect(assay(vsd),
                                             batch = vsd$sequencing_pool,
                                             batch2 = vsd$rna_purification_method,
                                             design = model.matrix(~vsd$gestation_week))

# PCA analysis (scale = FALSE)
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
    ggplot(aes(PC1, PC2, color=gestation_week)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Gestation\nWeek",
         title="miRBase + Novel",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme +
    scale_color_gradientn(colors = c("blue", "grey80", "red"))

ggsave(paste0(dir.pdf, "mirbase_novel_batch_corrected_pca.pdf"))

# remove seq. pool and rna purification method
vsd.trans.mirbase <- vsd.mirbase
assay(vsd.trans.mirbase) <- limma::removeBatchEffect(assay(vsd.mirbase),
                                             batch = vsd.mirbase$sequencing_pool,
                                             batch2 = vsd.mirbase$rna_purification_method,
                                             design = model.matrix(~vsd.mirbase$gestation_week))

# PCA analysis (scale = FALSE)
pca <- prcomp(t(assay(vsd.trans.mirbase)), center = TRUE, scale. = TRUE)
df <- data.frame(pca$x[,1:5], colData(vsd.trans.mirbase))
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
    ggplot(aes(PC1, PC2, color=gestation_week)) +
    geom_point(size=2) +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep=""),
         color="Gestation\nWeek",
         title="miRBase Only",
         caption="Batch Corrected: Seq. Pool and Purif. Method") +
    plotTheme +
    scale_color_gradientn(colors = c("blue", "grey80", "red"))

ggsave(paste0(dir.pdf, "mirbase_only_batch_corrected_pca.pdf"))
