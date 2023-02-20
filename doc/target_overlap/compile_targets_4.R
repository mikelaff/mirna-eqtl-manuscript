# compile predictions


library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(AnnotationHub)
library(miRNAtap)
library(EnsDb.Hsapiens.v86)


# INPUT FILES #########################################################################################################
# total RNA-seq (gene expression) DE data
df.gene.diff.expression.rds <- here("doc/diff_expression/rdata/20190506_diff_gene_expression_df.rds")
# small RNA-seq (miRNA expression) DE data
df.mirna.diff.expression.rds <- here("doc/diff_expression/rdata/20190506_diff_expression_df.rds")
# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.csv <- here("data/ucla_single_cell/TableS6 Cluster enriched genes.csv")

# Import Data #########################################################################################################
df.gene <- readRDS(df.gene.diff.expression.rds)
df.mirna <- readRDS(df.mirna.diff.expression.rds)

df.gene$gene_biotype <- factor(df.gene$gene_biotype)

df.genes.by.cluster <- read_csv(genes.by.cluster.csv)
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Filter Genes ########################################################################################################
# edb <- EnsDb.Hsapiens.v86
# prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

df.gene %<>%
    dplyr::filter(gene_biotype == "protein_coding")

df.mirna %<>%
    dplyr::filter(sig.gzcp | sig.gw)

edb <- EnsDb.Hsapiens.v86
prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

df.genes.by.cluster %<>%
    dplyr::filter(Ensembl %in% prot_genes_ensg)

# Label miRNAs and Genes by Category ##################################################################################
# label significant in both analyses column
df.mirna$sig.both[which(df.mirna$sig.gw & df.mirna$sig.gzcp)] <- "both"
df.mirna$sig.both[which(df.mirna$sig.gw & !df.mirna$sig.gzcp)] <- "gw only"
df.mirna$sig.both[which(!df.mirna$sig.gw & df.mirna$sig.gzcp)] <- "gzcp only"
df.mirna$sig.both[which(!df.mirna$sig.gw & !df.mirna$sig.gzcp)] <- "neither"
df.mirna$sig.both <- factor(df.mirna$sig.both)

# mirnas that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
miRNAs.maturation.late <- df.mirna$Name[which(df.mirna$sig.both == "gw only" & df.mirna$log2FoldChange.gw > 0)]
# mirnas that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
miRNAs.maturation.early <- df.mirna$Name[which(df.mirna$sig.both == "gw only" & df.mirna$log2FoldChange.gw < 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
miRNAs.neurogenesis.gz <- df.mirna$Name[which((df.mirna$sig.both == "both" | df.mirna$sig.both == "gzcp only") & df.mirna$log2FoldChange.gzcp > 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
miRNAs.neurogenesis.cp <- df.mirna$Name[which((df.mirna$sig.both == "both" | df.mirna$sig.both == "gzcp only") & df.mirna$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df.mirna$category <- "none"
df.mirna$category[which(df.mirna$Name %in% miRNAs.maturation.late)] <- "maturation_late"
df.mirna$category[which(df.mirna$Name %in% miRNAs.maturation.early)] <- "maturation_early"
df.mirna$category[which(df.mirna$Name %in% miRNAs.neurogenesis.gz)] <- "neurogenesis_gz"
df.mirna$category[which(df.mirna$Name %in% miRNAs.neurogenesis.cp)] <- "neurogenesis_cp"
# factor category variable
df.mirna$category <- factor(df.mirna$category)

# label significant in both analyses column
df.gene$sig.both <- NA
df.gene$sig.both[which(df.gene$sig.gw & df.gene$sig.gzcp)] <- "both"
df.gene$sig.both[which(df.gene$sig.gw & !df.gene$sig.gzcp)] <- "gw only"
df.gene$sig.both[which(!df.gene$sig.gw & df.gene$sig.gzcp)] <- "gzcp only"
df.gene$sig.both[which(!df.gene$sig.gw & !df.gene$sig.gzcp)] <- "neither"
df.gene$sig.both <- factor(df.gene$sig.both)

# genes that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
genes.maturation.late <- df.gene$gene_id[which(df.gene$sig.both == "gw only" & df.gene$log2FoldChange.gw > 0)]
# genes that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
genes.maturation.early <- df.gene$gene_id[which(df.gene$sig.both == "gw only" & df.gene$log2FoldChange.gw < 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
genes.neurogenesis.gz <- df.gene$gene_id[which((df.gene$sig.both == "both" | df.gene$sig.both == "gzcp only") & df.gene$log2FoldChange.gzcp > 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
genes.neurogenesis.cp <- df.gene$gene_id[which((df.gene$sig.both == "both" | df.gene$sig.both == "gzcp only") & df.gene$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df.gene$category <- "none"
df.gene$category[which(df.gene$gene_id %in% genes.maturation.late)] <- "maturation_late"
df.gene$category[which(df.gene$gene_id %in% genes.maturation.early)] <- "maturation_early"
df.gene$category[which(df.gene$gene_id %in% genes.neurogenesis.gz)] <- "neurogenesis_gz"
df.gene$category[which(df.gene$gene_id %in% genes.neurogenesis.cp)] <- "neurogenesis_cp"
# factor category variable
df.gene$category <- factor(df.gene$category)

# Find Targets ######################################
# minimum number of sources required for a target to be considered
min_sources <- 4

df.mirna$mir <- NA

df.mirna <- dplyr::filter(df.mirna, source == "miRBase_v22")

df.predictions <- data.frame()

for (i in 1:length(df.mirna$Name)) {

    # if novel, skip
    if (df.mirna$source[i] != "miRBase_v22") {
        next
    }


    predictions <- NULL
    rowdata <- NULL

    # short name used for database lookup
    df.mirna$mir[i] <- strsplit(strsplit(df.mirna$Name[i], '/')[[1]][1], "hsa-")[[1]][2]
    print(paste(i, df.mirna$mir[i], sep=": "))
    # get predicted targets
    predictions <- data.frame(getPredictedTargets(df.mirna$mir[i], species = "hsa", method = "geom", min_src = min_sources))
    if (length(predictions) == 0) {
        print(paste("No Targets:", df.mirna$mirna[i]))
        next
    }
    predictions$ENTREZID <- rownames(predictions)
    # get ensembl ids
    rowdata <- AnnotationDbi::select(orgdb,
                                     keys = predictions$ENTREZID,
                                     keytype = "ENTREZID",
                                     columns = c("ENSEMBL"))
    # join with predictions
    predictions <- left_join(predictions, rowdata, by = "ENTREZID")
    predictions$miRNA <- df.mirna$Name[i]

    df.predictions <- bind_rows(predictions, df.predictions)
}

df.predictions <- dplyr::rename(df.predictions, pictar = source_1, diana = source_2, targetscan = source_3, miranda = source_4, mirdb = source_5)

saveRDS(df.predictions, here("doc/target_overlap/rdata/20190514_miRNA_predictions_4sources.rds"))
write_csv(df.predictions, here("doc/target_overlap/20190514_miRNA_predictions_4sources.csv"))
