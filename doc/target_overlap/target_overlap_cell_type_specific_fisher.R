# look at miRNA targets and their overlap with differentially expressed genes, and genes defined
# as being enriched in a cell type by scRNA-seq (from Luis at UCLA)

# Input:
# differentially expressed miRNAs and genes from fetal tissue
# predicted targets from miRNAtap databases
# cell type specific genes from scRNA-seq

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(readxl)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

source(here("src/utils/lafferty_utils.R"))
date.prefix <- format(Sys.time(), "%Y%m%d")

date.prefix <- "20190519"
# OUTPUT FILES ########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/target_overlap/pdfs/")
# write plots
write.plots <- FALSE
# data frame with fisher test results 2,3,4 sources miRNAtap
df.enrichments.2sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_2sources.rds", sep=""))
df.enrichments.3sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_3sources.rds", sep=""))
df.enrichments.4sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_4sources.rds", sep=""))
# data frame with fisher test results for miRDB 2019 predicted targets
df.enrichments.mirdb.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRDB_enrichments.rds", sep=""))

# INPUT FILES #########################################################################################################
# total RNA-seq (gene expression) DE data
df.gene.diff.expression.rds <- here("doc/diff_expression/rdata/20190519_diff_gene_expression_df.rds")
# small RNA-seq (miRNA expression) DE data
df.mirna.diff.expression.rds <- here("doc/diff_expression/rdata/20190519_diff_expression_df.rds")
# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.xlsx <- here("data/ucla_single_cell/TableS4 Cluster analysis.xlsx")
# target predictions by miRNAtap databases
predictions.2sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_2sources.csv")
predictions.3sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_3sources.csv")
predictions.4sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_4sources.csv")
# miRDB 2019 predictions
predictions.mirdb.txt.gz <- here("data/target_predictions/miRDB_v6.0_prediction_result.txt.gz")

# adjusted p-value cuttoff to be considered significantly enriched or depleted
P.VALUE.CUTTOFF.ENRICHMENTS <- 0.05
# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Import DE Data ######################################################################################################
df.de.gene <- readRDS(df.gene.diff.expression.rds)
df.de.mirna <- readRDS(df.mirna.diff.expression.rds)

# all expressed genes
EXPRESSED.GENES <- df.de.gene$gene_id

# filter out non-protein-coding genes
df.de.gene$gene_biotype <- factor(df.de.gene$gene_biotype)
df.de.gene %<>%
  dplyr::filter(gene_biotype == "protein_coding")

# all expressed protein coding
EXPRESSED.GENES.PROT.CODING <- df.de.gene$gene_id

# filter for only significant DE genes and mirnas
df.de.gene %<>%
  dplyr::filter(sig.gzcp | sig.gw)
df.de.mirna %<>%
  dplyr::filter(sig.gzcp | sig.gw)

# filter for only miRBase annotated miRNAs
df.de.mirna.mirbase <- dplyr::filter(df.de.mirna, source == "miRBase_v22")

# Import Single Cell Gene Data ########################################################################################

df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

# Import Predicted Targets ############################################################################################

# miRDB 2019 predictions
df.mirdb.predictions <- read_tsv(predictions.mirdb.txt.gz, col_names = c("Name", "REFSEQ", "TargetScore", "ENTREZID"))
df.mirdb.predictions %<>% dplyr::filter(Name %in% df.de.mirna.mirbase$Name)

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.predictions$REFSEQ),
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL"))

# combine with predictions df
df.mirdb.predictions %<>%
  dplyr::left_join(rowdata, by = "REFSEQ")

df.mirdb.predictions %<>%
  dplyr::select(Name, ENSEMBL) %>%
  dplyr::filter(!is.na(ENSEMBL))

# miRNAtap compiled predictions
df.2sources.predictions <- read_csv(predictions.2sources.csv)
df.3sources.predictions <- read_csv(predictions.3sources.csv)
df.4sources.predictions <- read_csv(predictions.4sources.csv)

df.2sources.predictions %<>%
  dplyr::select(Name = miRNA, ENSEMBL) %>%
  dplyr::filter(!is.na(ENSEMBL))

df.3sources.predictions %<>%
  dplyr::select(Name = miRNA, ENSEMBL) %>%
  dplyr::filter(!is.na(ENSEMBL))

df.4sources.predictions %<>%
  dplyr::select(Name = miRNA, ENSEMBL) %>%
  dplyr::filter(!is.na(ENSEMBL))

# Gene Universes ######################################################################################################
# all protein-coding genes
edb <- EnsDb.Hsapiens.v86
GENES.PROT.CODING <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

# Find Target Overlaps ######################################

# Things to Set ###############
# which miRNAs: neurogenesis implicated mirnas (DE in GZ/CP or both datasets)
# which predictions: 2sources in miRNATap
# which gene universe: expressed genes, protein coding

# df.mirna.enrichments
#df.de.mirna.mirbase %>%
#  dplyr::filter(category == "NeuroGZ" | category == "NeuroCP") -> df.mirna.enrichments
df.mirna.enrichments <- df.de.mirna.mirbase

# df.targets
df.targets <- df.4sources.predictions

# GENE.UNIVERSE
GENE.UNIVERSE <- EXPRESSED.GENES.PROT.CODING

# fisher test statistics
columns.ft <- c("FTodds", "FTpval", "FTpvalAdj", "FTciLow", "FTciHigh", "FTenriched", "FTdepleted")

# loop over each cluster
for (cluster in clusters) {
  print(paste("############## Cluster: ", cluster, " ###############", sep=""))

  # initialize data frame for FT stats
  column.names <- c("Name", columns.ft)
  df.ft <- data.frame(matrix(nrow = length(df.mirna.enrichments$Name), ncol = length(column.names)))
  colnames(df.ft) <- column.names

  # genes for this cluster out of the gene universe
  genes.in.cluster <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == cluster)]
  # remove duplicates
  genes.in.cluster <- genes.in.cluster[!duplicated(genes.in.cluster)]
  # make sure genes are in gene universe
  genes.in.cluster <- genes.in.cluster[genes.in.cluster %in% GENE.UNIVERSE]
  # genes in universe that are not in cluster
  genes.not.in.cluster <- GENE.UNIVERSE[! GENE.UNIVERSE %in% genes.in.cluster]

  # check all genes are accounted for
  if ( length(GENE.UNIVERSE) != length(genes.in.cluster) + length(genes.not.in.cluster) ) {
    stop("GENE.UNIVERSE not same length as combination of genes in and out of cluster!")
  }

  # loop over each miRNA
  for (i in 1:length(df.mirna.enrichments$Name)) {
    mir <- df.mirna.enrichments$Name[i]
    df.ft$Name[i] <- mir

    #print(paste("miR: ", mir, sep=""))

    # if no targets exist, skip
    if (! mir %in% df.targets$Name) {
      next
    }

    # targets for this mirna
    targets <- dplyr::filter(df.targets, Name == mir)$ENSEMBL
    # remove duplicate targets
    targets <- targets[!duplicated(targets)]
    # make sure targets are in the gene universe
    targets <- targets[targets %in% GENE.UNIVERSE]
    # genes in universe that are not targets
    nonTargets <- GENE.UNIVERSE[! GENE.UNIVERSE %in% targets]

    # check all genes are accounted for
    if ( length(GENE.UNIVERSE) != length(targets) + length(nonTargets) ) {
      stop("GENE.UNIVERSE not same length as combination of targets and non targets!")
    }

    # contingency table
    #                   not.in.cluster                  /   in.cluster
    # nonTargets        num.nonTargets.not.in.cluster   /   num.nonTargets.in.cluster
    # targets           num.targets.not.in.cluster      /   num.targets.in.cluster

    num.nonTargets.not.in.cluster <- sum(nonTargets %in% genes.not.in.cluster)
    num.nonTargets.in.cluster <- sum(nonTargets %in% genes.in.cluster)
    num.targets.not.in.cluster <- sum(targets %in% genes.not.in.cluster)
    num.targets.in.cluster <- sum(targets %in% genes.in.cluster)

    # check all genes are accounted for
    if (length(GENE.UNIVERSE) != num.nonTargets.not.in.cluster + num.nonTargets.in.cluster +
        num.targets.not.in.cluster + num.targets.in.cluster) {
      stop("Contingency table incorrect!")
    }

    # Fisher Test for enrichment
    ft <- NULL
    ft <- fisher.test(matrix(c(num.nonTargets.not.in.cluster, num.nonTargets.in.cluster,
                               num.targets.not.in.cluster, num.targets.in.cluster),
                             byrow = TRUE,
                             nrow = 2))

    # Fill in df.ft of results
    df.ft$FTodds[i] <- ft$estimate
    df.ft$FTpval[i] <- ft$p.value
    df.ft$FTciLow[i] <- ft$conf.int[1]
    df.ft$FTciHigh[i] <- ft$conf.int[2]
  }

  # adjusted p-values
  df.ft$FTpvalAdj <- p.adjust(df.ft$FTpval, method = "fdr")

  # call enrichment or depletion
  df.ft$FTenriched <- (df.ft$FTpvalAdj <= P.VALUE.CUTTOFF.ENRICHMENTS) & (df.ft$FTodds > 1) & (df.ft$FTciLow > 1)
  df.ft$FTdepleted <- (df.ft$FTpvalAdj <= P.VALUE.CUTTOFF.ENRICHMENTS) & (df.ft$FTodds < 1) & (df.ft$FTciHigh < 1)

  # rename columns
  column.names <- paste(columns.ft, cluster, sep = "_")
  column.names <- c("Name", column.names)
  colnames(df.ft) <- column.names

  # join with df.mirna.enrichments
  df.mirna.enrichments %<>%
    dplyr::left_join(df.ft, by = "Name")
}

# Save Results ########################################################################################################
saveRDS(df.mirna.enrichments, df.enrichments.4sources.rds)

# Heatmap #############################################################################################################

df2 <- readRDS(df.enrichments.2sources.rds)
df3 <- readRDS(df.enrichments.3sources.rds)
df4 <- readRDS(df.enrichments.4sources.rds)
dfM <- readRDS(df.enrichments.mirdb.rds)

# Maturation Heatmap: miRDB ##########
df <- dfM
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "MatEarly" | category == "MatLate")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
                 "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gw), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gw))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                   ExNeuron = "#D95F02",
                                   InNeuron = "#7570B3",
                                   Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_mat_mirs_mirdb_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Maturation Associated miRNAs (miRDB2019, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Maturation Heatmap: 4sources ##########
df <- df4
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "MatEarly" | category == "MatLate")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gw), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gw))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_mat_mirs_4sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Maturation Associated miRNAs (4sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Maturation Heatmap: 3sources ##########
df <- df3
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "MatEarly" | category == "MatLate")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gw), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gw))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_mat_mirs_3sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Maturation Associated miRNAs (3sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Maturation Heatmap: 2sources ##########
df <- df2
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "MatEarly" | category == "MatLate")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gzcp <- ifelse(df$sig.gzcp, df$log2FoldChange.gzcp, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gw), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gw))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_mat_mirs_2sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Maturation Associated miRNAs (2sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Neurogenesis Heatmap: miRDB ##########
df <- dfM
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_neuro_mirs_mirdb_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Neuro. Associated miRNAs (miRDB2019, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Neurogenesis Heatmap: 4sources ##########
df <- df4
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_neuro_mirs_4sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Neuro. Associated miRNAs (4sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Neurogenesis Heatmap: 3sources ##########
df <- df3
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_neuro_mirs_3sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Neuro. Associated miRNAs (3sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()

# Neurogenesis Heatmap: 2sources ##########
df <- df2
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.odds <- as.matrix(dplyr::select(df, starts_with("FTodds")))
mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.odds) <- sapply(strsplit(colnames(mat.odds), "_"), `[`, 2)
colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)

all(colnames(mat.odds) == colnames(mat.enriched))
all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.odds)

# plotting matrix
mat <- mat.odds
# if not enriched, set odds to 1
mat[!mat.enriched] <- 1
#mat[is.na(mat)] <- 1

# transform odds ratios
mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp))
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(10)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

pdf(paste(dir.pdf, "heatmap_neuro_mirs_2sources_10counts.pdf", sep = "/"), useDingbats = FALSE, height=8, width=10)
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
         na_col = "grey50",
         breaks = seq(from=0, to=1, by=0.1),
         annotation_colors = ann_colors,
         main = "Target Enrichment of Neuro. Associated miRNAs (2sources, counts>10)",
         fontsize_row = 6)

#print(p)
dev.off()







