# Predicted targets of diff. expressed miRNAs:
# How many of the miRNA predicted targets are expressed?
# How many targets are differentially expressed?
# How many targets have expression levels negatively correlated with miRNA expression?
# Of the negatively correlated targets, are they cell-type specific?
# 3 June 2019

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
library(DESeq2)
#library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(reshape2)
library(corrplot)
library(psych)
#library(mikelaff)

source(here("src/utils/lafferty_utils.R"))
date.prefix <- format(Sys.time(), "%Y%m%d")

# OUTPUT FILES ########################################################################################################
# output directory for graphs
dir.pdf <- here("doc/target_expression/pdfs/")
dir.png <- here("doc/target_expression/pngs/")

# # data frame with fisher test results 2,3,4 sources miRNAtap
# df.enrichments.2sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_2sources.rds", sep=""))
# df.enrichments.3sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_3sources.rds", sep=""))
# df.enrichments.4sources.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRNAtap_enrichments_4sources.rds", sep=""))
# # data frame with fisher test results for miRDB 2019 predicted targets
# df.enrichments.mirdb.rds <- paste(paste(here("doc/target_overlap/rdata/"), date.prefix, "_miRDB_enrichments.rds", sep=""))

# INPUT FILES #########################################################################################################
# RangedSummarizedExperiment (rse) with known and novel counts. mirbaseV22 counts quantified by mirge2.0
# (mirbaseV22, friedlander, nowakowski, mirdeep2, mirge2.0)
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# RangedSummarizedExperiment (rse) with gene expression
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

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
# miRDB 2019 predictions, novel mirnas
predictions.mirdb.novel.csv <- here("data/target_predictions/miRDB_v6.0_novel_prediction_results.csv")

# GLOBALS #############################################################################################################
# write plots
WRITE.PLOTS <- TRUE

# Annotation Hub ######################################################################################################
# used to annotate prediction databases with ENSG number
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Import DE Data ######################################################################################################
df.de.gene <- readRDS(df.gene.diff.expression.rds)
df.de.mirna <- readRDS(df.mirna.diff.expression.rds)
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

# modify DE mirnas for correlation matrix
df.de.mirna.gzcp <- dplyr::filter(df.de.mirna, category  %in% c("NeuroCP", "NeuroGZ"))
df.de.mirna.gzcp %<>%
    dplyr::select(Name, category, Log2_fold_change = log2FoldChange.gzcp,
                  `P-value` = pvalue.gzcp, FDR = padj.gzcp)
df.de.mirna.gw <- dplyr::filter(df.de.mirna, category  %in% c("MatEarly", "MatLate"))
df.de.mirna.gw %<>%
    dplyr::select(Name, category, Log2_fold_change = log2FoldChange.gw,
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

# Import Single Cell Gene Data ########################################################################################

df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)
df.genes.by.cluster$Log2_fold_change <- as.numeric(df.genes.by.cluster$Log2_fold_change)
df.genes.by.cluster$`P-value` <- as.numeric(df.genes.by.cluster$`P-value`)
df.genes.by.cluster$FDR <- as.numeric(df.genes.by.cluster$FDR)
df.genes.by.cluster$Percent_expressed_cluster <- as.numeric(df.genes.by.cluster$Percent_expressed_cluster)
df.genes.by.cluster$Percent_expressed_all_cells <- as.numeric(df.genes.by.cluster$Percent_expressed_all_cells)

# append DE genes to cluster genes
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

# Plot Specific miRs ##################################################################################################
mir <- "hsa-miR-124-3p"

# predicted targets: miRDB2019
mirdb.targets <- dplyr::filter(df.mirdb.predictions, Name == mir)$ENSEMBL
mirdb.targets <- mirdb.targets[!duplicated(mirdb.targets)]
# filter for expressed targets
mirdb.targets <- mirdb.targets[mirdb.targets %in% rownames(vsd.gene)]
# filter for DE targets
mirdb.targets <- mirdb.targets[mirdb.targets %in% df.de.gene$gene_id]

# expression of this miR across samples
mir.expression <- assay(vsd.mirna)[mir,]

# data frame for plotting, rnaid per sample, miR expression per sample
df <- data.frame(rnaid = names(mir.expression), mir_expression = mir.expression, stringsAsFactors = FALSE)

# create data frame for target expression, rownames as rnaid
df.targets.expression <- as.data.frame(t(assay(vsd.gene[mirdb.targets,])))
df.targets.expression$rnaid <- rownames(df.targets.expression)

# combine with plotting data frame by rnaid
df %<>%
    dplyr::left_join(df.targets.expression, by = "rnaid")

# plot
df %>%
    dplyr::filter(rnaid %in% colData(rse.mirna)$rnaid[colData(rse.mirna)$tissue_section != "CW"]) %>%
    melt(id.vars = c("rnaid", "mir_expression")) %>%
    ggplot(aes(x=mir_expression, y=value, group=variable)) +
    geom_point(size=.1) +
    geom_smooth(method = "lm", se=FALSE, size=.2) +
    theme(legend.position = "none") +
    labs(y="Target Expression",
         x="miR-124-3p Expression",
         title="miRDB 2019 Predicted Targets: CW Samples Only") +
    plotTheme

ggsave(paste(dir.png, "mir124_expression_by_target_expression_cw_samples.png", sep="/"), height=8, width=8)


# Correlation Heatmaps ################################################################################################


mat.mirs <- assay(vsd.mirna)

mat.genes <- assay(vsd.gene)

all(colnames(mat.mirs) == colnames(mat.genes))

# correlation between all expressed mirnas and genes
# corr <- cor(x=t(mat.mirs), y=t(mat.genes), method = "pearson")
#
# cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(255)
#
# pheatmap(mat = corr,
#          color = cols)

# just gz/cp samples
mat.de.mirs.gzcp <- assay(vsd.mirna[rownames(vsd.mirna) %in% df.de.mirna.gzcp$Name, colData(vsd.mirna)$tissue_section != "CW"])

mat.de.genes.gzcp <- assay(vsd.gene[rownames(vsd.gene) %in% df.de.gene.gzcp$Ensembl, colData(vsd.gene)$tissue_section != "CW"])

all(colnames(mat.de.mirs.gzcp) == colnames(mat.de.genes.gzcp))

mat.corr.de.gzcp <- cor(x=t(mat.de.mirs.gzcp), y=t(mat.de.genes.gzcp), method="pearson")

# get correlations and significance values
corr.de.gzcp <- corr.test(x = t(mat.de.mirs.gzcp),
                          y = t(mat.de.genes.gzcp),
                          method = "pearson",
                          adjust = "fdr",
                          ci=FALSE)
# extract correlation
mat.corr.de.gzcp <- corr.de.gzcp$r
# extract corrected p-values
mat.pval.de.gzcp <- corr.de.gzcp$p

#mat.corr.de.gzcp[mat.pval.de.gzcp > 0.05] <- 0

cols <- colorRampPalette(c("darkorange", "white", "darkgreen"))(255)

cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(255)

png(paste(dir.png, "heatmap_pears_corr_de_mirna_genes_sig_only_gzcp.png", sep="/"),
    height = 800,
    width = 1200)
pheatmap(mat = mat.corr.de.gzcp,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Pearson Corr. of Diff. Expressed miRNAs and mRNAs (only significant values)")
dev.off()

png(paste(dir.png, "heatmap_pears_corr_de_mirna_genes_gzcp.png", sep="/"),
    height = 800,
    width = 1200)
pheatmap(mat = corr.de.gzcp$r,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Pearson Corr. of Diff. Expressed miRNAs and mRNAs")
dev.off()

# across gestation week
mat.de.mirs.gw <- assay(vsd.mirna[rownames(vsd.mirna) %in% df.de.mirna.gw$Name, colData(vsd.mirna)$tissue_section == "CW"])

mat.de.genes.gw <- assay(vsd.gene[rownames(vsd.gene) %in% df.de.gene.gw$Ensembl, colData(vsd.gene)$tissue_section == "CW"])

all(colnames(mat.de.mirs.gw) == colnames(mat.de.genes.gw))

# get correlations and significance values
corr.de.gw <- corr.test(x = t(mat.de.mirs.gw),
                          y = t(mat.de.genes.gw),
                          method = "pearson",
                          adjust = "fdr",
                          ci=FALSE)
# extract correlation
mat.corr.de.gw <- corr.de.gw$r
# extract corrected p-values
mat.pval.de.gw <- corr.de.gw$p

mat.corr.de.gw[mat.pval.de.gw > 0.05] <- 0

cols <- colorRampPalette(rev(brewer.pal(11,"RdBu")))(255)

png(paste(dir.png, "heatmap_pears_corr_de_mirna_genes_sig_only_gw.png", sep="/"),
    height = 800,
    width = 1200)
pheatmap(mat = mat.corr.de.gw,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Cortical Wall Samples:\nPearson Corr. of Diff. Expressed miRNAs and mRNAs (only significant values)")
dev.off()

png(paste(dir.png, "heatmap_pears_corr_de_mirna_genes_gw.png", sep="/"),
    height = 800,
    width = 1200)
pheatmap(mat = corr.de.gw$r,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE,
         main = "Cortical Wall Samples:\nPearson Corr. of Diff. Expressed miRNAs and mRNAs")
dev.off()

# Individual miRNA Expression Correlations ############################################################################

# just gz/cp samples, DE mirnas
mat.de.mirs.gzcp <- assay(vsd.mirna[rownames(vsd.mirna) %in% df.de.mirna.gzcp$Name, colData(vsd.mirna)$tissue_section != "CW"])
# gz/cp samples, DE genes
mat.de.genes.gzcp <- assay(vsd.gene[rownames(vsd.gene) %in% df.de.gene.gzcp$Ensembl, colData(vsd.gene)$tissue_section != "CW"])

all(colnames(mat.de.mirs.gzcp) == colnames(mat.de.genes.gzcp))

# get correlations and significance values
corr.de.gzcp <- corr.test(x = t(mat.de.mirs.gzcp),
                          y = t(mat.de.genes.gzcp),
                          method = "pearson",
                          adjust = "fdr",
                          ci=FALSE)

mir <- "hsa-miR-124-3p"

# predicted targets: miRDB2019
mirdb.targets <- dplyr::filter(df.mirdb.predictions, Name == mir)$ENSEMBL
mirdb.targets <- mirdb.targets[!duplicated(mirdb.targets)]
# filter for DE targets
mirdb.targets <- mirdb.targets[mirdb.targets %in% df.de.gene.gzcp$Ensembl]

mir.target.corrs <- corr.de.gzcp$r[mir,]
mir.target.pvals <- corr.de.gzcp$p[mir,]

# check both lists are in the same order (by ensg#), combine into dataframe
if (all(names(mir.target.corrs) == names(mir.target.pvals))) {
    df <- data.frame(ensg = names(mir.target.corrs), mir.target.corrs, mir.target.pvals)
} else {
    stop("Named lists not in same order")
}


df %>%
    ggplot(aes(x=mir.target.corrs, y=-log10(mir.target.pvals))) +
    geom_point() +
    plotTheme +
    labs(title="miR-124-3p Expression Correlation\nto DE mRNA Expression (GZ/CP)",
         x="Pearson's Correlation",
         y="-Log10(P-value) (FDR adjusted)") +
    geom_hline(yintercept = -log10(0.05))

ggsave(paste(dir.png, "mir124_expression_pearsons_corr_pvalue_gzcp.png", sep="/"), height=6, width=6)

# only DE mRNA targets
mir.target.corrs <- corr.de.gzcp$r[mir,mirdb.targets]
mir.target.pvals <- corr.de.gzcp$p[mir,mirdb.targets]

# check both lists are in the same order (by ensg#), combine into dataframe
if (all(names(mir.target.corrs) == names(mir.target.pvals))) {
    df <- data.frame(ensg = names(mir.target.corrs), mir.target.corrs, mir.target.pvals)
} else {
    stop("Named lists not in same order")
}


df %>%
    ggplot(aes(x=mir.target.corrs, y=-log10(mir.target.pvals))) +
    geom_point() +
    plotTheme +
    labs(title="miR-124-3p Expression Correlation\nto DE mRNA Expression (targets only)",
         x="Pearson's Correlation",
         y="P-value (FDR adjusted)") +
    geom_hline(yintercept = -log10(0.05))

ggsave(paste(dir.png, "mir124_expression_pearsons_corr_pvalue_gzcp_targets_only.png", sep="/"), height=6, width=6)

# Enrichments of Negatively Correlated Targets ########################################################################



# Linear Regression: miRNA expression vs target expression ############################################################

# just gz/cp samples, DE mirnas
mat.de.mirs.gzcp <- assay(vsd.mirna[rownames(vsd.mirna) %in% df.de.mirna.gzcp$Name, colData(vsd.mirna)$tissue_section != "CW"])
# gz/cp samples, DE genes
mat.de.genes.gzcp <- assay(vsd.gene[rownames(vsd.gene) %in% df.de.gene.gzcp$Ensembl, colData(vsd.gene)$tissue_section != "CW"])

all(colnames(mat.de.mirs.gzcp) == colnames(mat.de.genes.gzcp))

lm.de.gzcp <- lm(t(mat.de.genes.gzcp) ~ t(mat.de.mirs.gzcp))

mir1 <- mat.de.mirs.gzcp["hsa-miR-92b-3p",]
gene1 <- mat.de.genes.gzcp["ENSG00000174145",]

all(names(mir1) == names(gene1))

lm1 <- lm(gene1 ~ mir1)

plot(gene1 ~ mir1)
abline(lm1)



# Scratch ####

cols <- colorRampPalette(rev(brewer.pal(11,"RdGy")))(255)

png(paste("~/Desktop/", "heatmap_background_rdgy.png", sep="/"),
    height = 900,
    width = 1400)
pheatmap(mat = mat.corr.de.gzcp,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         show_colnames = FALSE)
dev.off()




