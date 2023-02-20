# look at miRNA targets and their overlap with differentially expressed genes, and genes defined
# as being enriched in a cell type by scRNA-seq (from Luis at UCLA)

# Subset target genes to ones that are diff. expressed and negatively correlated with miRNA expression

# Also look at global enrichments against DE genes (mRNAs) from tissue RNA-seq dataset

# High Threshold for cell type specific genes (top 20 in each cluster)
# Novel miRNA predicted targets by miRDB 2019

# Input:
# differentially expressed miRNAs and genes from fetal tissue
# predicted targets from miRNAtap databases
# cell type specific genes from scRNA-seq
# miRNA and gene expression

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(readxl)
library(AnnotationHub)
library(EnsDb.Hsapiens.v86)
library(DESeq2)
library(psych)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(mikelaff)

date.prefix <- format(Sys.time(), "%Y%m%d")

# OUTPUT FILES ########################################################################################################
# output directory for graphs
dir.pdf <- here("doc/target_overlap/pdfs/")
dir.png <- here("doc/target_overlap/pngs/")

# data frame with fisher test results for miRDB 2019 predicted targets with novel miRNAs neg. correlated with target expression
df.known.novel.enrichments.mirdb.rds <- paste(paste(here("doc/target_overlap/rdata/"),
                                                    date.prefix, "_miRDB_known_novel_enrichments_pos_corr.rds", sep=""))
# data frame with fisher test results 2,3,4 sources miRNAtap, neg. correlated with target expression
df.enrichments.2sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_2sources_pos_corr.rds", sep=""))
df.enrichments.3sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_3sources_pos_corr.rds", sep=""))
df.enrichments.4sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_4sources_pos_corr.rds", sep=""))

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

# target predictions by miRNAtap databases
predictions.2sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_2sources.csv")
predictions.3sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_3sources.csv")
predictions.4sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_4sources.csv")

# RangedSummarizedExperiment (rse) with known and novel counts. mirbaseV22 counts quantified by mirge2.0
# (mirbaseV22, friedlander, nowakowski, mirdeep2, mirge2.0)
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# RangedSummarizedExperiment (rse) with gene expression
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# GLOBALS #############################################################################################################
# write plots
WRITE.PLOTS <- TRUE

# adjusted p-value cuttoff to be considered significantly enriched or depleted
P.VALUE.CUTTOFF.ENRICHMENTS <- 0.05

# number of cluster genes by highest log2fold change
NUM.GENES.PER.CLUSTER <- 9999999

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

# modify DE genes for combining with cluster genes below
df.de.gene.gzcp <- dplyr::filter(df.de.gene, category  %in% c("NeuroCP", "NeuroGZ"))
df.de.gene.gzcp %<>%
  dplyr::select(Ensembl = gene_id, Gene = symbol, Cluster = category, Log2_fold_change = log2FoldChange.gzcp,
                `P-value` = pvalue.gzcp, FDR = padj.gzcp)
df.de.gene.gw <- dplyr::filter(df.de.gene, category  %in% c("MatEarly", "MatLate"))
df.de.gene.gw %<>%
  dplyr::select(Ensembl = gene_id, Gene = symbol, Cluster = category, Log2_fold_change = log2FoldChange.gw,
                `P-value` = pvalue.gw, FDR = padj.gw)

# import ranged summarized experiment for plotting expression
rse.gene <- readRDS(rse.gene.rds)
rse.mirna <- readRDS(rse.mirna.rds)

# count threshold
rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]
rse.mirna <- rse.mirna[rowSums(assay(rse.mirna) >= 10) >= 10, ]

dds.gene <- DESeqDataSet(rse.gene, design = ~1)
dds.mirna <- DESeqDataSet(rse.mirna, design = ~1)

vsd.gene <- vst(dds.gene)
vsd.mirna <- varianceStabilizingTransformation(dds.mirna)

# Get Expression Correlations #########################################################################################

# Subset for only DE miRNAs
mat.de.mirs <- assay(vsd.mirna[rownames(vsd.mirna) %in% df.de.mirna$Name,])
# Subset for only DE genes
mat.de.genes <- assay(vsd.gene)
# make sure both tables have samples in the same order
if (all(colnames(mat.de.mirs) == colnames(mat.de.genes))) {
  # get correlations and significance values
  corr.de <- corr.test(x = t(mat.de.mirs),
                       y = t(mat.de.genes),
                       method = "pearson",
                       adjust = "fdr",
                       ci=FALSE)
} else {
  stop("Correlations not completed.")
}

# Import Single Cell Gene Data ########################################################################################

df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)
df.genes.by.cluster$Log2_fold_change <- as.numeric(df.genes.by.cluster$Log2_fold_change)
df.genes.by.cluster$`P-value` <- as.numeric(df.genes.by.cluster$`P-value`)
df.genes.by.cluster$FDR <- as.numeric(df.genes.by.cluster$FDR)
df.genes.by.cluster$Percent_expressed_cluster <- as.numeric(df.genes.by.cluster$Percent_expressed_cluster)
df.genes.by.cluster$Percent_expressed_all_cells <- as.numeric(df.genes.by.cluster$Percent_expressed_all_cells)

# append DE genes to cluster genes for enrichment tests
df.genes.by.cluster %<>%
  dplyr::bind_rows(df.de.gene.gzcp, df.de.gene.gw)

df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)
clusters <- levels(df.genes.by.cluster$Cluster)

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
  dplyr::select(Name, ENSEMBL) %>%
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
  dplyr::select(Name, ENSEMBL) %>%
  dplyr::filter(!is.na(ENSEMBL))

# combine known and novel
df.mirdb.predictions %<>%
  dplyr::bind_rows(df.mirdb.predictions.novel)

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
# edb <- EnsDb.Hsapiens.v86
# GENES.PROT.CODING <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

# Find Target Overlaps ################################################################################################

# GENE.UNIVERSE
GENE.UNIVERSE <- EXPRESSED.GENES.PROT.CODING

# loop over the 4 prediction lists and save results
pred.list <- list(df.mirdb.predictions,
                  df.2sources.predictions,
                  df.3sources.predictions,
                  df.4sources.predictions)
save.list <- list(df.known.novel.enrichments.mirdb.rds,
                  df.enrichments.2sources.rds,
                  df.enrichments.3sources.rds,
                  df.enrichments.4sources.rds)

for (p in 1:length(pred.list)) {
  printMessage(fillChar = "%")
  printMessage(paste("Prediction list:", p, "of", length(pred.list)), fillChar = "%")
  printMessage(fillChar = "%")

  # set prediction list
  df.targets <- pred.list[[p]]

  geneSets <- c("all_prot_coding", "all_expressed", "all_expressed_prot_coding", "all_diff_expressed", "all_")
  # set enrichments dataframe
  df.geneSet.enrichments <- data.frame(geneSet = c())

  # fisher test statistics column names
  columns.ft <- c("FTodds", "FTpval", "FTpvalAdj", "FTciLow", "FTciHigh", "FTenriched", "FTdepleted")

  # loop over each cluster
  for (cluster in clusters) {
    printMessage(paste("Cluster: ", cluster, sep=""))

    # initialize data frame for FT stats
    column.names <- c("Name", columns.ft)
    df.ft <- data.frame(matrix(nrow = length(df.mirna.enrichments$Name), ncol = length(column.names)))
    colnames(df.ft) <- column.names

    # genes for this cluster out of the gene universe
    genes.in.cluster <- top_n(dplyr::filter(df.genes.by.cluster, Cluster == cluster), NUM.GENES.PER.CLUSTER, wt = Log2_fold_change)$Ensembl
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

      # targets significantly negatively correlated with miRNA expression
      corrs <- corr.de$r[mir, targets]
      pvals <- corr.de$p[mir, targets]

      if (all(names(corrs) == names(pvals))) {
        # set correlations to 0 if pvalue not significant
        corrs[pvals > 0.05] <- 0
        # select as targets only the genes that are negatively correlated
        targets <- names(corrs[corrs < 0])
      } else {
        stop("Corrs and Pvals don't have same names.")
      }


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

  # save results
  print("Saving results...")
  saveRDS(df.mirna.enrichments, save.list[[p]])
}


stop("Stop before plotting")

# Heatmaps ############################################################################################################

dfM <- readRDS(here("doc/target_overlap/rdata/20190614_miRDB_known_novel_enrichments_neg_corr.rds"))

source(here("doc/target_overlap/utils/target_overlap_heatmap.R"))


png(paste(dir.png, "heatmap_enrich_deplet_matur_mirs_mirdb_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(dfM, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, miRDB2019")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_neuro_mirs_mirdb_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(dfM, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, miRDB2019")
dev.off()

# pdfs
pdf(paste(dir.pdf, "heatmap_enrich_deplet_matur_mirs_mirdb_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(dfM, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, miRDB2019")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_neuro_mirs_mirdb_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(dfM, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, miRDB2019")
dev.off()

# miRNAtap
df2 <- readRDS(here("doc/target_overlap/rdata/20190614_miRNAtap_enrichments_2sources_neg_corr.rds"))
df3 <- readRDS(here("doc/target_overlap/rdata/20190614_miRNAtap_enrichments_3sources_neg_corr.rds"))
df4 <- readRDS(here("doc/target_overlap/rdata/20190614_miRNAtap_enrichments_4sources_neg_corr.rds"))

source(here("doc/target_overlap/utils/target_overlap_heatmap.R"))

# pngs
png(paste(dir.png, "heatmap_enrich_deplet_matur_mirs_2sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df2, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 2sources")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_neuro_mirs_2sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df2, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 2sources")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_matur_mirs_3sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df3, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 3sources")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_neuro_mirs_3sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df3, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 3sources")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_matur_mirs_4sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df4, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 4sources")
dev.off()

png(paste(dir.png, "heatmap_enrich_deplet_neuro_mirs_4sources_10counts_known_novel_neg_corr.png", sep="/"),
    height = 800,
    width = 800)
plot.enrichments(df4, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 4sources")
dev.off()

# pdfs
pdf(paste(dir.pdf, "heatmap_enrich_deplet_matur_mirs_2sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df2, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 2sources")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_neuro_mirs_2sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df2, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 2sources")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_matur_mirs_3sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df3, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 3sources")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_neuro_mirs_3sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df3, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 3sources")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_matur_mirs_4sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df4, categories = c("MatEarly", "MatLate"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Maturation\nAssociated miRNAs, 4sources")
dev.off()

pdf(paste(dir.pdf, "heatmap_enrich_deplet_neuro_mirs_4sources_10counts_known_novel_neg_corr.pdf", sep="/"),
    height = 10,
    width = 8,
    useDingbats = FALSE)
plot.enrichments(df4, categories = c("NeuroGZ", "NeuroCP"),
                 title = "Neg. Correlated Target Enrichment/Depletion of Neurogenesis\nAssociated miRNAs, 4sources")
dev.off()

