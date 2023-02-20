# looking at different kinship matrices

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)
library(pheatmap)
library(RColorBrewer)
library(ggplot2)

# OUTPUT ###############################################################################################################
dir.pdf <- here("doc/emmax/pdfs/")

# INPUT ################################################################################################################
dir.kinship.mats <- "~/Projects/BACKUP_EXCLUDED/"

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# small rna-seq samples used in mirQTL analysis
small.qtl.samples.txt <- here("results/emmax/samples/20191101_mirQTL_RNAID_DonorID_DNAID.txt")

# GLOBALS ##############################################################################################################

# Load Sample Metadata #################################################################################################

# samples used in mirQTL analysis
mirQTL.samples <- read_delim(small.qtl.samples.txt, delim = " ", col_names = c("RNAID", "DonorID", "DNAID"))

# sample metadata
samples <- as.data.frame(read_tsv(samples.metadata.tsv))
# filter for mirQTLor samples
samples %<>%
    filter(RNAID %in% mirQTL.samples$RNAID)




# Load Kinship Matrices ################################################################################################

tfam <- read_delim("~/Projects/BACKUP_EXCLUDED/chrAll.mirQTL.tfam",
                   delim = " ",
                   col_names = c("DonorID", "DNAID", "ID3", "ID4", "Sex", "Phenotype"))

kmat.chrALL <- as.matrix(read_tsv("~/Projects/BACKUP_EXCLUDED/chrAll.mirQTL.BN.kinf", col_names = FALSE))
rownames(kmat.chrALL) <- tfam$DonorID
colnames(kmat.chrALL) <- tfam$DonorID


rownames(samples) <- samples$DonorID

# Plot Heatmaps ###############

cols <- rev(colorRampPalette((brewer.pal(11,"RdBu")))(250))

df.anno <- data.frame(Sex = samples$Sex.by.XIST, row.names = samples$DonorID)
cols.anno <- list(Sex = c(`F` = "red", `M` = "blue"))


pheatmap(mat = kmat.chrALL,
         color = cols,
         breaks = seq(0,1,length.out = 250),
         annotation_row = df.anno,
         annotation_col = df.anno,
         annotation_colors = cols.anno,
         show_colnames = FALSE,
         show_rownames = FALSE,
         main = "chrAll")



kmat.not23.26 <- as.matrix(read_tsv("~/Projects/BACKUP_EXCLUDED/chrAll.not23-26.mirQTL.BN.kinf", col_names = FALSE))
rownames(kmat.not23.26) <- tfam$DonorID
colnames(kmat.not23.26) <- tfam$DonorID

pheatmap(mat = kmat.not23.26,
         color = cols,
         breaks = seq(0,1,length.out = 250),
         annotation_row = df.anno,
         annotation_col = df.anno,
         annotation_colors = cols.anno,
         show_colnames = FALSE,
         show_rownames = FALSE,
         main = "Not23-26")

kmat.not1 <- as.matrix(read_tsv("~/Projects/BACKUP_EXCLUDED/chrAll.not1.mirQTL.BN.kinf", col_names = FALSE))
rownames(kmat.not1) <- tfam$DonorID
colnames(kmat.not1) <- tfam$DonorID

pheatmap(mat = kmat.not1,
         color = cols,
         breaks = seq(0,1,length.out = 250),
         annotation_row = df.anno,
         annotation_col = df.anno,
         annotation_colors = cols.anno,
         show_colnames = FALSE,
         show_rownames = FALSE,
         main = "Not1")

kmat.not1.not23.26 <- as.matrix(read_tsv("~/Projects/BACKUP_EXCLUDED/chrAll.not1.not23-26.mirQTL.BN.kinf", col_names = FALSE))
rownames(kmat.not1.not23.26) <- tfam$DonorID
colnames(kmat.not1.not23.26) <- tfam$DonorID

pheatmap(mat = kmat.not1.not23.26,
         color = cols,
         breaks = seq(0,1,length.out = 250),
         annotation_row = df.anno,
         annotation_col = df.anno,
         annotation_colors = cols.anno,
         show_colnames = FALSE,
         show_rownames = FALSE,
         main = "Not1 Not23-26")




