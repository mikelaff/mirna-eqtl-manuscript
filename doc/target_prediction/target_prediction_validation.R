# how well does target prediction work?
# use linear regression or correlation?
# does target prediction predict correlation?

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
# # output directory for graphs
dir.pdf <- here("doc/target_prediction/pdfs/")
dir.png <- here("doc/target_prediction/pngs/")
#
# # data frame with results for miRDB 2019 predicted targets with novel miRNAs
# df.known.novel.enrichments.mirdb.rds <- paste(paste(here("doc/target_prediction/rdata/"),
#                                                     date.prefix, "_miRDB_results.rds", sep=""))
# # data frame with results 2,3,4 sources miRNAtap
# df.enrichments.2sources.rds <- paste(paste(here("doc/target_prediction/rdata/"),
#                                            date.prefix, "_miRNAtap_2sources_results.rds", sep=""))
# df.enrichments.3sources.rds <- paste(paste(here("doc/target_prediction/rdata/"),
#                                            date.prefix, "_miRNAtap_3sources_results.rds", sep=""))
# df.enrichments.4sources.rds <- paste(paste(here("doc/target_prediction/rdata/"),
#                                            date.prefix, "_miRNAtap_4sources_results.rds", sep=""))

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
#WRITE.PLOTS <- TRUE

# adjusted p-value cuttoff
#P.VALUE.CUTTOFF <- 0.05

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
#
# df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
# df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)
# df.genes.by.cluster$Log2_fold_change <- as.numeric(df.genes.by.cluster$Log2_fold_change)
# df.genes.by.cluster$`P-value` <- as.numeric(df.genes.by.cluster$`P-value`)
# df.genes.by.cluster$FDR <- as.numeric(df.genes.by.cluster$FDR)
# df.genes.by.cluster$Percent_expressed_cluster <- as.numeric(df.genes.by.cluster$Percent_expressed_cluster)
# df.genes.by.cluster$Percent_expressed_all_cells <- as.numeric(df.genes.by.cluster$Percent_expressed_all_cells)
#
# # append DE genes to cluster genes for enrichment tests
# df.genes.by.cluster %<>%
#     dplyr::bind_rows(df.de.gene.gzcp, df.de.gene.gw)
#
# df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)
# clusters <- levels(df.genes.by.cluster$Cluster)

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

# Scratch #############################################################################################################

mir <- "hsa-miR-92b-3p"

# targets for this mirna
targets <- dplyr::filter(df.mirdb.predictions, Name == mir)$ENSEMBL
# remove duplicate targets
targets <- targets[!duplicated(targets)]

# dataframe of targets and correlations to this mir
df.targets <- data.frame(Ensembl = targets, stringsAsFactors = FALSE)

# add correlations
df.corrs <- data.frame(corr = t(corr.de$r)[,mir])
df.corrs$Ensembl <- rownames(df.corrs)
df.targets %<>% dplyr::left_join(df.corrs, by = "Ensembl")

# plot most negatively correlated
target <- "ENSG00000087586"

df.target.expr <- data.frame(target.expr = mat.de.genes[target,])
df.target.expr$rnaid <- rownames(df.target.expr)

df.mir.expr <- data.frame(mir.expr = mat.de.mirs[mir,])
df.mir.expr$rnaid <- rownames(df.mir.expr)

df <- dplyr::full_join(df.mir.expr, df.target.expr, by = "rnaid")

df %>%
    ggplot(aes(mir.expr, target.expr)) +
    geom_point() +
    geom_smooth(method = "lm") +
    plotTheme()

corr <- corr.de$r[mir, target]
corr.pval <- corr.de$p[mir, target]

lm.test <- lm(target.expr ~ mir.expr, data = df)

df %>%
    ggplot(aes(mir.expr, target.expr)) +
    geom_point() +
    geom_abline(intercept = coef(lm.test)[1], slope = coef(lm.test)[2]) +
    plotTheme(theme = "presentation") +
    labs(x=paste(mir, "Expression"),
         y=paste("Target Expression:", target),
         title="VST Normalized Expression") +
    geom_smooth(method = "lm")

ggsave("~/Desktop/92b_pos_corr_expression.pdf", height=7, width=7)

 # Scratch ############

# for each miR, does correlation predict target prediction
# predicted target ~ corr + len

# GENE.UNIVERSE
GENE.UNIVERSE <- EXPRESSED.GENES.PROT.CODING

df.transcripts <- readRDS(here("doc/target_overlap/rdata/tmp.rds"))

# set prediction list
df.targets <- df.4sources.predictions

# set dataframe
df.mirna.pred.validation <- df.de.mirna

# logistic regression statistics column names
columns.ft <- c("LGlogit", "LGoddsRatio", "LGprob", "LGpval", "LGpvalAdj", "LGsig", "LGenriched", "LGdepleted")

# initialize data frame for LG stats
column.names <- c("Name", columns.ft)
df.lg <- data.frame(matrix(nrow = length(df.mirna.pred.validation$Name), ncol = length(column.names)))
colnames(df.lg) <- column.names

# loop over each miRNA
for (i in 1:length(df.mirna.pred.validation$Name)) {
    mir <- df.mirna.pred.validation$Name[i]
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

    # target genes get a 1
    df.genes$target <- ifelse(df.genes$Ensembl %in% targets, 1, 0)

    # convert to factor
    df.genes$target <- as.factor(df.genes$target)

    # add correlations
    df.corrs <- data.frame(corr = t(corr.de$r)[,mir])
    df.corrs$Ensembl <- rownames(df.corrs)
    df.genes %<>% dplyr::left_join(df.corrs, by = "Ensembl")

    # add gene length information
    df.genes %<>% dplyr::left_join(df.transcripts, by = "Ensembl")

    lg <- glm(target ~ corr + log10tx, data = df.genes, family = binomial)

    coef.lg <- coef(lg)

    # Fill in df.lg: coefficient for "corr"
    df.lg$LGlogit[i] <- coef.lg["corr"]
    df.lg$LGoddsRatio[i] <- exp(coef.lg)["corr"]
    df.lg$LGprob[i] <- (exp(coef.lg) / (1 + exp(coef.lg)))["corr"]
    df.lg$LGpval[i] <- coef(summary(lg))["corr", 4]
}

binomial_smooth <- function(...) {
    geom_smooth(method = "glm", method.args = list(family = "binomial"), ...)
}

df.genes %>%
    ggplot(aes(corr, target)) +
    geom_point() +
    labs(title=mir) +
    plotTheme()

ggsave(paste(dir.png, mir, "_target_corr.png", sep=""), width = 6, height=6)

lm.test <- lm(corr ~ target, data = df.genes)

df.genes %>%
    ggplot(aes(target, corr)) +
    geom_point() +
    plotTheme() +
    labs(title=mir) +
    geom_abline(intercept = coef(lm.test)[1], slope = coef(lm.test)[2])

ggsave(paste(dir.png, mir, "_corr_target.png", sep=""), width = 6, height=6)


# Scratch ###########3

# targets for this mirna
targets <- dplyr::filter(df.2sources.predictions, Name == "hsa-miR-124-3p")$ENSEMBL
# remove duplicate targets
targets <- targets[!duplicated(targets)]
#1720

# targets for this mirna
targets <- dplyr::filter(df.3sources.predictions, Name == "hsa-miR-124-3p")$ENSEMBL
# remove duplicate targets
targets <- targets[!duplicated(targets)]
#1029

# targets for this mirna
targets <- dplyr::filter(df.4sources.predictions, Name == "hsa-miR-124-3p")$ENSEMBL
# remove duplicate targets
targets <- targets[!duplicated(targets)]
#413

# targets for this mirna
targets <- dplyr::filter(df.mirdb.predictions, Name == "chr19_38431_star")$ENSEMBL
# remove duplicate targets
targets <- targets[!duplicated(targets)]
#1749

df.mirdb.predictions %>%
    dplyr::group_by(Name) %>%
    dplyr::summarise(count = n()) -> tmp

