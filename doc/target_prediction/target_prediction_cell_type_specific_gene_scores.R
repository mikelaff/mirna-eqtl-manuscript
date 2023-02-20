# Investigate target prediction scores by cell-type specific genes (do prediction of neuron specific
# genes have higher scores than progenitors? possible reason for enrichment bias?)
# 29 May 2019

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(reshape2)
library(AnnotationHub)
library(beepr)
library(pheatmap)
library(RColorBrewer)

library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ########################################################################################################
# output directory for png files
dir.png <- here("doc/target_prediction/pngs/")

# INPUT FILES #########################################################################################################
# total RNA-seq (gene expression) DE data
df.gene.diff.expression.rds <- here("doc/diff_expression/rdata/20190519_diff_gene_expression_df.rds")
# small RNA-seq (miRNA expression) DE data
df.mirna.diff.expression.rds <- here("doc/diff_expression/rdata/20190519_diff_expression_df.rds")

# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.xlsx <- here("data/ucla_single_cell/TableS4 Cluster analysis.xlsx")

# miRDB 2019 predictions, known mirnas
predictions.mirdb.txt.gz <- here("data/target_predictions/miRDB_v6.0_prediction_result.txt.gz")
# miRDB 2019 predictions, novel mirnas
predictions.mirdb.novel.csv <- here("data/target_predictions/miRDB_v6.0_novel_prediction_results.csv")

# GLOBALS #############################################################################################################
# write plots
WRITE.PLOTS <- TRUE
# top n number of genes by log2fold change to define as cluster specific
NUM.GENES.PER.CLUSTER <- 99999

# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Import Single Cell Gene Data ########################################################################################

df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

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

# Import Predicted Targets ############################################################################################

# miRDB 2019 predictions for novel miRNAs
df.mirdb.predictions.novel <- read_csv(predictions.mirdb.novel.csv)

# remove NA symbols
df.mirdb.predictions.novel %<>% dplyr::filter(!is.na(GeneSymbol))
df.mirdb.predictions.novel %<>% dplyr::rename(SYMBOL = GeneSymbol)
df.mirdb.predictions.novel %<>% dplyr::rename(Name = miRNAName)

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.predictions.novel$SYMBOL),
                                 keytype = "SYMBOL",
                                 columns = c("ENSEMBL"))

# join with predictions
df.mirdb.predictions.novel %<>%
    dplyr::left_join(rowdata, by = "SYMBOL")

df.mirdb.predictions.novel %<>%
    dplyr::select(Name, ENSEMBL, TargetScore) %>%
    dplyr::filter(!is.na(ENSEMBL))

# miRDB 2019 predictions
df.mirdb.predictions <- read_tsv(predictions.mirdb.txt.gz, col_names = c("Name", "REFSEQ", "TargetScore", "ENTREZID"))
df.mirdb.predictions %<>% dplyr::filter(Name %in% df.de.mirna$Name)

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.predictions$REFSEQ),
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL"))

# combine with predictions df
df.mirdb.predictions %<>%
    dplyr::left_join(rowdata, by = "REFSEQ")

df.mirdb.predictions %<>%
    dplyr::select(Name, ENSEMBL, TargetScore) %>%
    dplyr::filter(!is.na(ENSEMBL))

# combine known and novel
df.mirdb.predictions %<>%
    dplyr::bind_rows(df.mirdb.predictions.novel)

rm(df.mirdb.predictions.novel)

# Find Mean TargetScore ###############################################################################################
# data frame to fill with mean scores
df.mirna.scores <- df.de.mirna

# df.targets
df.targets <- df.mirdb.predictions

# score columns
columns.ft <- c("mean", "median", "max", "min", "var")

# loop over each cluster
for (cluster in clusters) {
    print(paste("############## Cluster: ", cluster, " ###############", sep=""))

    # initialize data frame for stats
    column.names <- c("Name", columns.ft)
    df.ft <- data.frame(matrix(nrow = length(df.mirna.scores$Name), ncol = length(column.names)))
    colnames(df.ft) <- column.names

    # genes for this cluster
    genes.in.cluster <- top_n(dplyr::filter(df.genes.by.cluster, Cluster == cluster), NUM.GENES.PER.CLUSTER, wt = Log2_fold_change)$Ensembl
    # remove duplicates
    genes.in.cluster <- genes.in.cluster[!duplicated(genes.in.cluster)]

    # loop over each miRNA
    for (i in 1:length(df.mirna.scores$Name)) {
        mir <- df.mirna.scores$Name[i]
        df.ft$Name[i] <- mir

        #print(paste("miR: ", mir, sep=""))

        # if no targets exist, skip
        if (! mir %in% df.targets$Name) {
            next
        }

        targets <- NULL
        # targets for this mirna
        targets <- dplyr::filter(df.targets, Name == mir)
        # target overlaping with cluster specific genes
        targets <- dplyr::filter(targets, ENSEMBL %in% genes.in.cluster)

        # if no overlaps exist, skip
        if (!length(targets$Name) > 0) {
            # Fill in df.ft of results
            df.ft$mean[i] <- 0
            df.ft$median[i] <- 0
            df.ft$max[i] <- 0
            df.ft$min[i] <- 0
            df.ft$var[i] <- 0
        } else {
            # Fill in df.ft of results
            df.ft$mean[i] <- mean(targets$TargetScore)
            df.ft$median[i] <- median(targets$TargetScore)
            df.ft$max[i] <- max(targets$TargetScore)
            df.ft$min[i] <- min(targets$TargetScore)
            df.ft$var[i] <- var(targets$TargetScore)
        }
    }

    # rename columns
    column.names <- paste(columns.ft, cluster, sep = "_")
    column.names <- c("Name", column.names)
    colnames(df.ft) <- column.names

    # join with df.mirna.enrichments
    df.mirna.scores %<>%
        dplyr::left_join(df.ft, by = "Name")
}

beep(sound = 8)

# Heatmap ##########
df <- df.mirna.scores
# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
#df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
#df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.mean <- as.matrix(dplyr::select(df, starts_with("mean")))
mat.median <- as.matrix(dplyr::select(df, starts_with("median")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.mean) <- sapply(strsplit(colnames(mat.mean), "_"), `[`, 2)
colnames(mat.median) <- sapply(strsplit(colnames(mat.median), "_"), `[`, 2)

#all(colnames(mat.odds) == colnames(mat.enriched))
#all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.median)

# plotting matrix
mat <- mat.median
# if not enriched, set odds to 1
#mat[!mat.enriched] <- 1
mat[is.na(mat)] <- 0

# transform odds ratios
#mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp),
                      Source = df$source)
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(100)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"),
                   Source = c(miRBase_v22 = "#66A61E",
                              Friedlander2014 = "#E6AB02",
                              miRDeep2 = "#BF5B17",
                              miRge = "#386CB0",
                              Nowakowski2018 = "grey70"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

#rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

png(paste(dir.png, "heatmap_median_scores_neuro_mirs_mirdb_10counts_known_novel.png", sep = "/"), height=800, width=800)
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
         breaks = seq(from=1, to=100, by=1),
         annotation_colors = ann_colors,
         main = "Median Target Scores of Neurogenesis\nAssociated miRNAs (miRDB2019, counts>10)",
         fontsize_row = 4,
         gaps_col = c(5,10,12))

#print(p)
dev.off()

dfE <- readRDS(here("doc/target_overlap/rdata/20190521_miRDB_known_novel_enrichments.rds"))
df <- df.mirna.scores

df %<>%
    dplyr::left_join(dplyr::select(dfE, Name, starts_with("FT")), by = "Name")

# filter for only "category" mirnas
df <- dplyr::filter(df, category == "NeuroGZ" | category == "NeuroCP")
# filter for high expressers
df <- dplyr::filter(df, baseMean.gw > 10)
# remove mirnas with no enrichments
df <- dplyr::filter(df, rowSums(dplyr::select(df, starts_with("FTenriched"))) > 0)
# set l2fold change to 0 if not significant
#df$log2FoldChange.gw <- ifelse(df$sig.gw, df$log2FoldChange.gw, 0)

df <- as.data.frame(df)
rownames(df) <- df$Name
colnames(df)

mat.enriched <- as.matrix(dplyr::select(df, starts_with("FTenriched")))
mat.median <- as.matrix(dplyr::select(df, starts_with("median")))

#rownames(mat.odds) <- df$Name
#rownames(mat.enriched) <- df$Name

colnames(mat.enriched) <- sapply(strsplit(colnames(mat.enriched), "_"), `[`, 2)
colnames(mat.median) <- sapply(strsplit(colnames(mat.median), "_"), `[`, 2)

#all(colnames(mat.odds) == colnames(mat.enriched))
#all(rownames(mat.odds) == rownames(mat.enriched))

# reorder cols
#mat.enriched <- mat.enriched[,colnames(mat.odds)]

# df for cluster labels
CellType <- c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor")
df.cluster <- data.frame(CellType = CellType)
df.cluster$CellType <- as.factor(df.cluster$CellType)
rownames(df.cluster) <- colnames(mat.median)

# plotting matrix
mat <- mat.median
# if not enriched, set to 0
mat[!mat.enriched] <- 0
mat[is.na(mat)] <- 0

# transform odds ratios
#mat <- log10(mat)

# order matrix
mat <- mat[order(df$log2FoldChange.gzcp), c(14,15,16,12,9,6,5,3,4,2,7,8,1,10,11,13)]

df.mirs <- data.frame(L2FC_GZ_CP = df$log2FoldChange.gzcp,
                      L2FC_GW = df$log2FoldChange.gw,
                      L2Expression = log2(df$baseMean.gzcp),
                      Source = df$source)
rownames(df.mirs) <- df$Name

cols <- colorRampPalette((brewer.pal(9,"Reds")))(10)
cols <- colorRampPalette(c("white", "red3"))(100)

ann_colors <- list(L2FC_GZ_CP = colorRampPalette(rev(brewer.pal(5,"RdBu")))(10),
                   L2FC_GW = colorRampPalette(brewer.pal(5,"PRGn"))(10),
                   L2Expression = colorRampPalette(brewer.pal(7,"Greys"), bias=1)(10),
                   CellType = c(Other = "#1B9E77",
                                ExNeuron = "#D95F02",
                                InNeuron = "#7570B3",
                                Progenitor = "#E7298A"),
                   Source = c(miRBase_v22 = "#66A61E",
                              Friedlander2014 = "#E6AB02",
                              miRDeep2 = "#BF5B17",
                              miRge = "#386CB0",
                              Nowakowski2018 = "grey70"))

#colorRampPalette(brewer.pal(4, "Dark2"))(4)

#rowNames <- ifelse(rownames(mat) == "hsa-miR-92b-3p" | rownames(mat) == "hsa-miR-124-3p", rownames(mat), "")

png(paste(dir.png, "heatmap_median_scores_whenEnriched_neuro_mirs_mirdb_10counts_known_novel.png", sep = "/"), height=800, width=800)
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
         breaks = seq(from=1, to=100, by=1),
         annotation_colors = ann_colors,
         main = "Median Target Scores of Neurogenesis\nAssociated miRNAs (miRDB2019, counts>10)",
         fontsize_row = 4,
         gaps_col = c(5,10,12))

#print(p)
dev.off()
