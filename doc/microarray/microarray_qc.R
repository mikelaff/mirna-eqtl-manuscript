
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
library(ggplot2)

library(ArrayExpress)
library(Biobase)
library(pd.clariom.s.human.ht)
library(clariomshumantranscriptcluster.db)
#library(pd.hugene.1.0.st.v1)
#library(hugene10sttranscriptcluster.db)
#library(AffyCompatible)
library(oligo)
library(arrayQualityMetrics)

library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)

library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)

library(mikelaffr)


# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/microarray/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

dir.quality.metrics <- here("doc/microarray/quality_metrics/")
dir.create(dir.quality.metrics, recursive = TRUE, showWarnings = FALSE)

# background removed, normalized, and filtered expression set
output.expressionSet.rds <- here("results/rdata_files/20220829_es_HNP_Differential_Microarray_SST-RMA.rds")

# INPUT ################################################################################################################
# Microarray ExpressionSet for HNP Differentiataion Experiment
expressionSet.rds <- here("results/rdata_files/20220829_es_HNP_Differentiation_Microarray.rds")

# yep
#tmp.txt <- here("results/microarray/Stein_Mike_Clariom_S_Human_24_08172022/log2.SST-RMA-GENE-FULL - Group 1.TXT")

# GLOBALS ##############################################################################################################

# Import ExpressionSet #################################################################################################
raw.es <- readRDS(expressionSet.rds)

# # yep
# temp.log2.expression <- read_tsv(tmp.txt)
#
# log2.expression %<>%
#     select(1:25)
#
# samples <- colnames(log2.expression)[2:25]
# colnames(log2.expression)[2:25] <- paste0(sapply(strsplit(samples, "\\."), `[`, 1), ".CEL")
#
# df.expression <- as.data.frame(log2.expression)
# rownames(df.expression) <- df.expression[,1]
#
# df.expression <- df.expression[,2:25]

# QC  Before Correction #################
# simple log2 of expression
log2.expr <- log2(exprs(raw.es))

pca <- prcomp(t(log2.expr), center = TRUE, scale. = FALSE)

# combine with phenotype data
df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(raw.es), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


df.pca %<>%
    mutate(Condition = paste(Expression, Timepoint, sep = "_"))

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
    geom_point() +
    geom_line() +
    labs(x="PC", y="Percent Total Variance",
         title="PCA: Microarray Expression (log2 Norm.)") +
    scale_x_continuous(breaks = seq(1:10)) +
    plotTheme()

#colnames(df.pca)

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

#ggsave(paste0(dir.pdfs, "pca_log2_norm.pdf"), height = 6, width = 8)

#pdf(paste0(dir.pdfs, "pca_by_technical.pdf"))

# df.pca %>%
#     ggplot(aes(`CRNA Yield (ug)`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`Extraction Date`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`Extraction Date`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(Expression, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(Replicate, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(Timepoint, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`ng/ul`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(A260, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(A280, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`260/280`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`260/230`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`340 raw`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`Total RNA (ng)`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`Core Sample ID`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(RIN, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`28S/18S (Area)`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`TapStation Conc. ng/ul`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(`Qubit Conc. ng/ul`, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")
#
# df.pca %>%
#     ggplot(aes(Condition, PC1)) +
#     geom_point(size=3) +
#     labs(y=paste("PC1 (", round(percentVar[1], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Gene Expression (log2 Norm.)")


#dev.off()

rm(df.pca, log2.expr, pca, percentVar)

# sort(abs(pca$rotation[,1]), decreasing = TRUE)[1:20]
#
# sort(abs(pca$rotation[,2]), decreasing = TRUE)[1:20]

# log2 intensities
oligo::boxplot(raw.es)

# quality metrics
# arrayQualityMetrics(expressionset = raw.es,
#                     outdir = dir.quality.metrics,
#                     force = TRUE,
#                     do.logtransform = TRUE,
#                     intgroup = c("Expression", "Timepoint"))

# RLE #################
# Relative Log Expression

# rle.es <- rma(raw.es, normalize = FALSE)
#
# row_medians_assayData <- Biobase::rowMedians(as.matrix(exprs(rle.es)))
#
# RLE_data <- sweep(exprs(rle.es), 1, row_medians_assayData)
#
# RLE_data <- as.data.frame(RLE_data)
# RLE_data_gathered <- gather(RLE_data, sample, log2_expression_deviation)
#
# ggplot(RLE_data_gathered, aes(sample,log2_expression_deviation)) +
#     geom_boxplot(outlier.shape = NA) +
#     ylim(c(-2, 2)) +
#     theme(axis.text.x = element_text(colour = "aquamarine4",
#                                      angle = 60, size = 6.5, hjust = 1 ,
#                                      face = "bold"))

# rm(rle.es, RLE_data, RLE_data_gathered, row_medians_assayData)
# RMA Calibration ##################

norm.es <- rma(raw.es)

pca <- prcomp(t(exprs(norm.es)), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(norm.es), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


df.pca %<>%
    mutate(Condition = paste(Expression, Timepoint, sep = "_"))

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
    geom_point() +
    geom_line() +
    labs(x="PC", y="Percent Total Variance",
         title="PCA: Microarray Expression (log2 Norm.)") +
    scale_x_continuous(breaks = seq(1:10)) +
    plotTheme()

#colnames(df.pca)

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

rm(pca, df.pca, percentVar)
# Intensity Filtering ##############

# median expression per probe
norm.medians <- rowMedians(exprs(norm.es))

hist(norm.medians, 100, freq = FALSE)

# manual intensity threshold
man_threshold <- 4

abline(v = man_threshold, col = "red", lwd = 2)

# number of samples per group
samples_cutoff <- 6

# indexes of probes with greater than threshold in at least cutoff number of samples
idx_man_threshold <- apply(exprs(norm.es), 1,
                           function(x){
                               sum(x > man_threshold) >= samples_cutoff})
table(idx_man_threshold)

# filter for expressed probes
filtered.es <- subset(norm.es, idx_man_threshold)

rm(raw.es, norm.es, idx_man_threshold, norm.medians, man_threshold, samples_cutoff)

# Cluster Annotation #############

annotations <- AnnotationDbi::select(clariomshumantranscriptcluster.db,
                                     keys = featureNames(filtered.es),
                                     columns = c("SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID"),
                                     keytype = "PROBEID")

# remove probes without a gene symbol or other ID
annotations %<>%
    dplyr::filter(!(is.na(SYMBOL) & is.na(GENENAME) & is.na(ENSEMBL) & is.na(ENTREZID)))

# probes with more than one annotation, created duplicate probeids
annotations$duplicate.annotation <- duplicated(annotations$PROBEID) | duplicated(annotations$PROBEID, fromLast = TRUE)
# remove all but one duplicate probeid
annotations %<>%
    dplyr::filter(!(duplicated(PROBEID)))

# filter expression set for probes with an annotation
final.es <- subset(filtered.es, featureNames(filtered.es) %in% annotations$PROBEID)

# add annotation to expression set
fData(final.es)$PROBEID <- rownames(fData(final.es))
fData(final.es) <- left_join(fData(final.es), annotations, by = "PROBEID")

rownames(fData(final.es)) <- fData(final.es)$PROBEID

validObject(final.es)

rm(annotations, filtered.es)


pca <- prcomp(t(exprs(final.es)), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(final.es), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


df.pca %<>%
    mutate(Condition = paste(Expression, Timepoint, sep = "_"))

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
    geom_point() +
    geom_line() +
    labs(x="PC", y="Percent Total Variance",
         title="PCA: Microarray Expression (log2 Norm.)") +
    scale_x_continuous(breaks = seq(1:10)) +
    plotTheme()

#colnames(df.pca)

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

rm(pca, df.pca, percentVar)


# Remove batch effect #################

corrected.expression <- limma::removeBatchEffect(final.es,
                                                  covariates = final.es$`CRNA Yield (ug)`,
                                                 design = model.matrix(~final.es$Timepoint + final.es$Expression))

pca <- prcomp(t(corrected.expression), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(final.es), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


df.pca %<>%
    mutate(Condition = paste(Expression, Timepoint, sep = "_"))

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
    geom_point() +
    geom_line() +
    labs(x="PC", y="Percent Total Variance",
         title="PCA: Microarray Expression (log2 Norm.)") +
    scale_x_continuous(breaks = seq(1:10)) +
    plotTheme()

#colnames(df.pca)

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

rm(pca, df.pca, percentVar)

# possible outliers
pData(final.es)$outlier <- FALSE
pData(final.es)$outlier[pData(final.es)$CEL %in% c("MA13.CEL", "MA21.CEL")] <- TRUE

filt.es <- final.es[,!pData(final.es)$outlier]

corrected.expression <- limma::removeBatchEffect(filt.es,
                                                 covariates = filt.es$`CRNA Yield (ug)`,
                                                 design = model.matrix(~filt.es$Timepoint + filt.es$Expression))

pca <- prcomp(t(corrected.expression), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(filt.es), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100


df.pca %<>%
    mutate(Condition = paste(Expression, Timepoint, sep = "_"))

# scree plot
ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
    geom_point() +
    geom_line() +
    labs(x="PC", y="Percent Total Variance",
         title="PCA: Microarray Expression (log2 Norm.)") +
    scale_x_continuous(breaks = seq(1:10)) +
    plotTheme()

#colnames(df.pca)

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

df.pca %>%
    ggplot(aes(PC1, PC2, color = Qubit.Conc..ng.ul)) +
    geom_point(size = 3) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme() +
    #scale_color_manual(values=cbPalette) +
    labs(title="Microarray Expression (log2 Norm.)")

# Output ##########
# final es
saveRDS(final.es, output.expressionSet.rds)


