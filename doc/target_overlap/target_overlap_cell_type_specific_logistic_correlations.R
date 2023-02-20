# look at miRNA targets and their overlap with differentially expressed genes, and genes defined
# as being enriched in a cell type by scRNA-seq (from Luis at UCLA) using logistic regression
# to correct for covariates (transcript length, gc bias, etc.)

# Novel miRNA predicted targets by miRDB 2019

# Input:
# differentially expressed miRNAs and genes from fetal tissue
# predicted targets from miRNAtap databases
# cell type specific genes from scRNA-seq
# miRNA and gene expression

time.start <- Sys.time()

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
library(beepr)

date.prefix <- format(Sys.time(), "%Y%m%d")

# OUTPUT FILES ########################################################################################################
# output directory for graphs
dir.pdf <- here("doc/target_overlap/pdfs/")
dir.png <- here("doc/target_overlap/pngs/")

# data frame with logistic regression results from miRDB2019 predictions
df.logistic.reg.mirdb.rds <- paste(paste(here("doc/target_overlap/rdata/"),
                                                    date.prefix, "_miRDB_logistic_regression_corr.rds", sep=""))
# # data frame with logistic regression results 2,3,4 sources miRNAtap
# df.logistic.reg.2sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix,
#                                             "_miRNAtap_2sources_logistic_regression_tx.rds", sep=""))
# df.logistic.reg.3sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix,
#                                             "_miRNAtap_3sources_logistic_regression_tx.rds", sep=""))
# df.logistic.reg.4sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix,
#                                             "_miRNAtap_4sources_logistic_regression_tx.rds", sep=""))

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
# all expressed genes
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

# Find Target Overlaps: Logistic Regression ###########################################################################

# GENE.UNIVERSE
GENE.UNIVERSE <- EXPRESSED.GENES.PROT.CODING

pred.list <- list(df.mirdb.predictions)
save.list <- list(df.logistic.reg.mirdb.rds)

# # loop over the 4 prediction lists and save results
# pred.list <- list(df.mirdb.predictions,
#                   df.2sources.predictions,
#                   df.3sources.predictions,
#                   df.4sources.predictions)
# save.list <- list(df.logistic.reg.mirdb.rds,
#                   df.logistic.reg.2sources.rds,
#                   df.logistic.reg.3sources.rds,
#                   df.logistic.reg.4sources.rds)

# # get transcript lengths, cds, 5'utr, and 3'utr lenghts
# df.transcripts <- transcriptLengths(EnsDb.Hsapiens.v86,
#                                     with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
# # summarize lenghts by gene
# df.transcripts %<>%
#   dplyr::filter(gene_id %in% GENE.UNIVERSE) %>%
#   dplyr::group_by(gene_id) %>%
#   dplyr::summarise_at(.vars = c("nexon", "tx_len", "cds_len", "utr5_len", "utr3_len"),
#                       .funs = c(mean=mean, median=median, min=min, max=max, var=var)) %>%
#   dplyr::rename(Ensembl = gene_id)
#
# df.transcripts$log10tx <- log10(df.transcripts$tx_len_mean)
#
# saveRDS(df.transcripts, here("doc/target_overlap/rdata/tmp.rds"))
df.transcripts <- readRDS(here("doc/target_overlap/rdata/tmp.rds"))

for (p in 1:length(pred.list)) {
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")
  print(paste("Prediction list:", p, "of", length(pred.list)))
  print(Sys.time() - time.start)
  print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")

  # set prediction list
  df.targets <- pred.list[[p]]

  # set enrichments dataframe
  df.mirna.enrichments <- df.de.mirna

  # logistic regression statistics column names
  columns.ft <- c("LGlogit", "LGoddsRatio", "LGprob", "LGpval", "LGpvalAdj", "LGsig", "LGenriched", "LGdepleted")

  # loop over each cluster
  for (cluster in clusters) {
    print(paste("############## Cluster: ", cluster, " ###############", sep=""))
    print(Sys.time() - time.start)

    # initialize data frame for LG stats
    column.names <- c("Name", columns.ft)
    df.lg <- data.frame(matrix(nrow = length(df.mirna.enrichments$Name), ncol = length(column.names)))
    colnames(df.lg) <- column.names

    # genes for this cluster out of the gene universe
    genes.in.cluster <- top_n(dplyr::filter(df.genes.by.cluster, Cluster == cluster), NUM.GENES.PER.CLUSTER, wt = Log2_fold_change)$Ensembl
    # remove duplicates
    genes.in.cluster <- genes.in.cluster[!duplicated(genes.in.cluster)]
    # make sure genes are in gene universe
    genes.in.cluster <- genes.in.cluster[genes.in.cluster %in% GENE.UNIVERSE]

    # loop over each miRNA
    for (i in 1:length(df.mirna.enrichments$Name)) {
      mir <- df.mirna.enrichments$Name[i]
      df.lg$Name[i] <- mir

      #print(paste("miR ", i, " of ", length(df.mirna.enrichments$Name), ": ", mir, sep=""))

      # if mir not in target list, skip
      if (! mir %in% df.targets$Name) {
        next
      }

      # targets for this mirna
      targets <- dplyr::filter(df.targets, Name == mir)$ENSEMBL
      # remove duplicate targets
      targets <- targets[!duplicated(targets)]
      # make sure targets are in the gene universe
      targets <- targets[targets %in% GENE.UNIVERSE]

      # if no targets, skip
      if (!length(targets) > 0) {
        next
      }

      # dataframe, one row per gene, for logistic regression fitting
      df.genes <- data.frame(Ensembl = GENE.UNIVERSE, stringsAsFactors = FALSE)

      # cluster specific genes get a 1
      df.genes$cluster_specific <- ifelse(df.genes$Ensembl %in% genes.in.cluster, 1, 0)
      # target genes get a 1
      df.genes$target <- ifelse(df.genes$Ensembl %in% targets, 1, 0)

      # convert to factor
      df.genes$cluster_specific <- as.factor(df.genes$cluster_specific)
      df.genes$target <- as.factor(df.genes$target)

      # add correlations
      df.corrs <- data.frame(corr = t(corr.de$r)[,mir])
      df.corrs$Ensembl <- rownames(df.corrs)
      df.genes %<>% dplyr::left_join(df.corrs, by = "Ensembl")

      # add gene length information
      df.genes %<>% dplyr::left_join(df.transcripts, by = "Ensembl")

      lg <- glm(cluster_specific ~ target + log10tx + corr, data = df.genes, family = binomial)

      coef.lg <- coef(lg)

      # Fill in df.lg: coefficient for "target1"
      df.lg$LGlogit[i] <- coef.lg["target1"]
      df.lg$LGoddsRatio[i] <- exp(coef.lg)["target1"]
      df.lg$LGprob[i] <- (exp(coef.lg) / (1 + exp(coef.lg)))["target1"]
      df.lg$LGpval[i] <- coef(summary(lg))["target1", 4]
    }

    # adjusted p-values
    df.lg$LGpvalAdj <- p.adjust(df.lg$LGpval, method = "fdr")

    # significant target beta value
    df.lg$LGsig <- df.lg$LGpvalAdj <= P.VALUE.CUTTOFF.ENRICHMENTS

    # call enrichment or depletion
    df.lg$LGenriched <- (df.lg$LGsig) & (df.lg$LGoddsRatio > 1)
    df.lg$LGdepleted <- (df.lg$LGsig) & (df.lg$LGoddsRatio < 1)

    # rename columns
    column.names <- paste(columns.ft, cluster, sep = "_")
    column.names <- c("Name", column.names)
    colnames(df.lg) <- column.names

    # join with df.mirna.enrichments
    df.mirna.enrichments %<>%
      dplyr::left_join(df.lg, by = "Name")
  }

  # save results
  print("Saving results...")
  saveRDS(df.mirna.enrichments, save.list[[p]])
}


stop("Stop before plotting")

# Heatmaps ############################################################################################################

source(here("doc/target_overlap/utils/target_overlap_heatmap.R"))

dfM.lg <- readRDS(here("doc/target_overlap/rdata/20190625_miRDB_logistic_regression_corr.rds"))
dfM.ft <- readRDS(here("doc/target_overlap/rdata/20190530_miRDB_known_novel_enrichments.rds"))

pdf(paste(dir.pdf, "heatmap_miRDB2019_ft_lg_comparison_log10tx_corr.pdf"), height=10, width=10, useDingbats = FALSE)
plot.enrichments(dfM.ft, categories = c("NeuroGZ", "NeuroCP"), title = "Fisher, miRDB2019")
plot.lg.enrichments(dfM.lg, categories = c("NeuroGZ", "NeuroCP"), title = "Log. Regression (+ log10 mean tx length + corr), miRDB2019")
dev.off()



