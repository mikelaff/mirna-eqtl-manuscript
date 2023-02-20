# Which cell types for each miRNA?

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)


df <- readRDS(here("doc/diff_expression/rdata/20190514_diff_expression_evidence_df.rds"))

colnames(df)

df.mirbase <- dplyr::filter(df, source == "miRBase_v22")

mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))

rownames(mat.odds) <- df.mirbase$Name
rownames(mat.enriched) <- df.mirbase$Name

oddsNames <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 3)
oddsNames[1] <- "GZ>CP"
oddsNames[2] <- "CP>GZ"
oddsNames[3] <- "Early"
oddsNames[4] <- "Late"
colnames(mat.odds) <- oddsNames

enrichNames <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
enrichNames[1] <- "GZ>CP"
enrichNames[2] <- "CP>GZ"
enrichNames[3] <- "Early"
enrichNames[4] <- "Late"
colnames(mat.enriched) <- enrichNames

# reorder cols
mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
ClusterType <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
                                          "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
                                          "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.cluster <- data.frame(ClusterType = ClusterType[5:20])
rownames(df.cluster) <- oddsNames[5:20]

#####

mat <- mat.odds
mat[!mat.enriched] <- 1
mat[is.na(mat)] <- 1

mat <- mat[,5:20]
mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gzcp), c(14,15,16,12,9,1,6,5,4,3,7,8,2,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp,
                      MeanExp_GZ_CP = df.mirbase$baseMean.gzcp)
rownames(df.mirs) <- df.mirbase$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50))

pdf("~/Desktop/heatmap_all_mirbase_mirs_by_l2fc.pdf", useDingbats = FALSE, height=10, width=15)
# row cluster order: progenitors, exNeuron, inneuron, other
# miRNA with targets enriched specific in oRG and not vRG
pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_col = df.cluster,
         annotation_row = df.mirs,
         annotation_color = ann_colors) -> p

print(p)
dev.off()


# YEP ##############
df.mirbase <- dplyr::filter(df, source == "miRBase_v22", category == "neurogenesis_cp" | category == "neurogenesis_gz", baseMean.gzcp > 200)

mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))

rownames(mat.odds) <- df.mirbase$Name
rownames(mat.enriched) <- df.mirbase$Name

oddsNames <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 3)
oddsNames[1] <- "GZ>CP"
oddsNames[2] <- "CP>GZ"
oddsNames[3] <- "Early"
oddsNames[4] <- "Late"
colnames(mat.odds) <- oddsNames

enrichNames <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
enrichNames[1] <- "GZ>CP"
enrichNames[2] <- "CP>GZ"
enrichNames[3] <- "Early"
enrichNames[4] <- "Late"
colnames(mat.enriched) <- enrichNames

# reorder cols
mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
ClusterType <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
                 "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
                 "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.cluster <- data.frame(ClusterType = ClusterType[5:20])
rownames(df.cluster) <- oddsNames[5:20]

#####

mat <- mat.odds
mat[!mat.enriched] <- 1
mat[is.na(mat)] <- 1

mat <- mat[,5:20]
mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gzcp), c(14,15,16,12,9,1,6,5,4,3,7,8,2,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp,
                      MeanExp_GZ_CP = df.mirbase$baseMean.gzcp)
rownames(df.mirs) <- df.mirbase$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50))

pdf("~/Desktop/heatmap_all_mirbase_mirs_by_l2fc_highExp.pdf", useDingbats = FALSE, height=10, width=15)
# row cluster order: progenitors, exNeuron, inneuron, other
# miRNA with targets enriched specific in oRG and not vRG
pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = df.cluster,
         annotation_row = df.mirs,
         annotation_color = ann_colors) -> p

print(p)
dev.off()

# YEP2 ####################################
df.mirbase <- dplyr::filter(df, source == "miRBase_v22", category == "maturation_early" | category == "maturation_late", baseMean.gw > 200)

mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))

rownames(mat.odds) <- df.mirbase$Name
rownames(mat.enriched) <- df.mirbase$Name

oddsNames <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 3)
oddsNames[1] <- "GZ>CP"
oddsNames[2] <- "CP>GZ"
oddsNames[3] <- "Early"
oddsNames[4] <- "Late"
colnames(mat.odds) <- oddsNames

enrichNames <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
enrichNames[1] <- "GZ>CP"
enrichNames[2] <- "CP>GZ"
enrichNames[3] <- "Early"
enrichNames[4] <- "Late"
colnames(mat.enriched) <- enrichNames

# reorder cols
mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
ClusterType <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
                 "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
                 "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.cluster <- data.frame(ClusterType = ClusterType[5:20])
rownames(df.cluster) <- oddsNames[5:20]

#####

mat <- mat.odds
mat[!mat.enriched] <- 1
mat[is.na(mat)] <- 1

mat <- mat[,5:20]
mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gw), c(14,15,16,12,9,1,6,5,4,3,7,8,2,10,11,13)]

df.mirs <- data.frame(L2FC_GW = df.mirbase$log2FoldChange.gw,
                      MeanExp_GW = df.mirbase$baseMean.gw)
rownames(df.mirs) <- df.mirbase$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50))

pdf("~/Desktop/heatmap_all_mirbase_mirs_by_l2fc_highExp_Maturation.pdf", useDingbats = FALSE, height=10, width=15)
# row cluster order: progenitors, exNeuron, inneuron, other
# miRNA with targets enriched specific in oRG and not vRG
pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = df.cluster,
         annotation_row = df.mirs,
         annotation_color = ann_colors) -> p

print(p)
dev.off()


###################

df.mirbase <- dplyr::filter(df, source == "miRBase_v22", category == "neurogenesis_cp" | category == "neurogenesis_gz", baseMean.gzcp > 200)

mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))

rownames(mat.odds) <- df.mirbase$Name
rownames(mat.enriched) <- df.mirbase$Name

oddsNames <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 3)
oddsNames[1] <- "GZ>CP"
oddsNames[2] <- "CP>GZ"
oddsNames[3] <- "Early"
oddsNames[4] <- "Late"
colnames(mat.odds) <- oddsNames

enrichNames <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
enrichNames[1] <- "GZ>CP"
enrichNames[2] <- "CP>GZ"
enrichNames[3] <- "Early"
enrichNames[4] <- "Late"
colnames(mat.enriched) <- enrichNames

# reorder cols
mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
ClusterType <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
                 "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
                 "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.cluster <- data.frame(ClusterType = ClusterType[5:20])
rownames(df.cluster) <- oddsNames[5:20]

#####

mat <- mat.odds
mat[!mat.enriched] <- 1
mat[is.na(mat)] <- 1

mat <- mat[,5:20]
mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gzcp), c(14,15,16,12,9,1,6,5,4,3,7,8,2,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp,
                      MeanExp_GZ_CP = df.mirbase$baseMean.gzcp)
rownames(df.mirs) <- df.mirbase$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(8,"RdBu")))(50))

pdf("~/Desktop/heatmap_all_mirbase_mirs_by_l2fc_highExp.pdf", useDingbats = FALSE, height=10, width=15)
# row cluster order: progenitors, exNeuron, inneuron, other
# miRNA with targets enriched specific in oRG and not vRG
pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = df.cluster,
         annotation_row = df.mirs,
         annotation_color = ann_colors) -> p

print(p)
dev.off()


stop()



























































#################################################
#####

df.clust <- df.clusters[5:20,]

mat <- mat.enriched[,df.clust$ClusterType]

rownames(mat.odds) <- df.mirbase$Name
mat.odds[is.na(mat.odds)] <- 1

colNames <- colnames(mat.odds)
colNames <- sapply(strsplit(colNames, "_"), `[`, 3)
colNames[1] <- "GZ>CP"
colNames[2] <- "CP>GZ"
colNames[3] <- "Early"
colNames[4] <- "Late"

colTypes <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
              "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
              "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.clusters <- data.frame(colTypes)
rownames(df.clusters) <- colNames

colnames(mat.odds) <- colNames

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp)
rownames(df.mirs) <- df.mirbase$Name

#mat <- (sweep(sweep(mat.odds, 2, colMins(mat.odds), "-"), 2, (colMaxs(mat.odds) - colMins(mat.odds)), "/"))

mat <- mat.odds
mat[mat == 0] <- 1.0

mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gzcp),]

cols <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)

pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_col = df.clusters,
         annotation_row = df.mirs
         )


# binary enrichment
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))
rownames(mat.enriched) <- df.mirbase$Name
mat.enriched[is.na(mat.enriched)] <- FALSE

colNames <- colnames(mat.enriched)
colNames <- sapply(strsplit(colNames, "_"), `[`, 2)
colNames[1] <- "GZ>CP"
colNames[2] <- "CP>GZ"
colNames[3] <- "Early"
colNames[4] <- "Late"

colNames

colTypes <- c("Tissue", "Tissue", "Tissue", "Tissue", "Other", "ExNeuron", "ExNeuron",
              "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
              "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.clusters <- data.frame(colTypes)
rownames(df.clusters) <- colNames

colnames(mat.enriched) <- colNames

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp)
rownames(df.mirs) <- df.mirbase$Name

mat <- mat.enriched
storage.mode(mat) <- "numeric"

mat <- mat[order(df.mirbase$log2FoldChange.gzcp),]

cols <- colorRampPalette(rev(brewer.pal(3,"Reds")))(2)

pheatmap(mat = mat,
         color = c("white", "red"),
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE,
         annotation_col = df.clusters,
         annotation_row = df.mirs
)

# filter for category
df.mirbase <- dplyr::filter(df, source == "miRBase_v22", category == "neurogenesis_cp" | category == "neurogenesis_gz", baseMean.gzcp > 100)
mat.enriched <- as.matrix(dplyr::select(df.mirbase, starts_with("enriched")))
mat.enriched[is.na(mat.enriched)] <- FALSE

mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
rownames(mat.odds) <- df.mirbase$Name
mat.odds[is.na(mat.odds)] <- 1
mat.odds[!mat.enriched] <- 1

colNames <- colnames(mat.odds)
colNames <- sapply(strsplit(colNames, "_"), `[`, 3)
colNames[1] <- "GZ>CP"
colNames[2] <- "CP>GZ"
colNames[3] <- "Early"
colNames[4] <- "Late"

colTypes <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
              "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
              "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.clusters <- data.frame(colTypes[5:20])
rownames(df.clusters) <- colNames[5:20]

colnames(mat.odds) <- colNames

df.mirs <- data.frame(L2FC_GZ_CP = df.mirbase$log2FoldChange.gzcp,
                      MeanExp_GZCP = df.mirbase$baseMean.gzcp)
rownames(df.mirs) <- df.mirbase$Name

#mat <- (sweep(sweep(mat.odds, 2, colMins(mat.odds), "-"), 2, (colMaxs(mat.odds) - colMins(mat.odds)), "/"))

mat <- mat.odds[,5:20]
mat[mat == 0] <- 1.0

mat <- log10(mat)
mat[mat == 0] <- NA

mat <- mat[order(df.mirbase$log2FoldChange.gzcp),]

cols <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(50))

# row cluster order: progenitors, exNeuron, inneuron, other
# miRNA with targets enriched specific in oRG and not vRG
pheatmap(mat = mat,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = df.clusters,
         annotation_row = df.mirs,
         annotation_color = ann_colors,
         na_col = "white"
)

# filter for category: maturation
df.mirbase <- dplyr::filter(df, source == "miRBase_v22", category == "maturation_early" | category == "maturation_late", baseMean.gw > 500)
mat.odds <- as.matrix(dplyr::select(df.mirbase, starts_with("FT_odds")))
rownames(mat.odds) <- df.mirbase$Name
mat.odds[is.na(mat.odds)] <- 1

colNames <- colnames(mat.odds)
colNames <- sapply(strsplit(colNames, "_"), `[`, 3)
colNames[1] <- "GZ>CP"
colNames[2] <- "CP>GZ"
colNames[3] <- "Early"
colNames[4] <- "Late"

colTypes <- c("Tissue", "Tissue", "Tissue", "Tissue", "ExNeuron", "Other", "ExNeuron",
              "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron", "Progenitor",
              "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")

df.clusters <- data.frame(colTypes)
rownames(df.clusters) <- colNames

colnames(mat.odds) <- colNames

df.mirs <- data.frame(L2FC_GW = df.mirbase$log2FoldChange.gw,
                      MeanExp_GW = df.mirbase$baseMean.gw)
rownames(df.mirs) <- df.mirbase$Name

#mat <- (sweep(sweep(mat.odds, 2, colMins(mat.odds), "-"), 2, (colMaxs(mat.odds) - colMins(mat.odds)), "/"))

mat <- mat.odds
mat[mat == 0] <- 1.0

mat <- log10(mat)

mat <- mat[order(df.mirbase$log2FoldChange.gw),]

cols <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)


ann_colors <- list(L2FC_GW = colorRampPalette(rev(brewer.pal(5,"RdBu")))(50))
pheatmap(mat = mat,
         color = cols,
         cluster_rows = FALSE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = TRUE,
         annotation_col = df.clusters,
         annotation_row = df.mirs,
         annotation_colors = ann_colors
)





