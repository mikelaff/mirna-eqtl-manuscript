
library(here)
library(readr)
library(dplyr)
library(tidyr)
library(readxl)
library(magrittr)
library(ggplot2)
library(ggrepel)

#library(ArrayExpress)
library(Biobase)
library(pd.clariom.s.human.ht)
library(clariomshumantranscriptcluster.db)
#library(pd.hugene.1.0.st.v1)
#library(hugene10sttranscriptcluster.db)
#library(AffyCompatible)
library(oligo)
#library(arrayQualityMetrics)

library(limma)
library(topGO)
library(org.Hs.eg.db)
#library(ReactomePA)
#library(clusterProfiler)

#library(gplots)
#library(ggplot2)
#library(geneplotter)
#library(RColorBrewer)
#library(pheatmap)
#library(enrichplot)

#library(stringr)
#library(matrixStats)
library(genefilter)
#library(openxlsx)

library(multiMiR)
library(miRNAtap)

library(ComplexHeatmap)

library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/microarray/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

dir.topgo <- here("doc/microarray/topgo/")
dir.create(dir.topgo, recursive = TRUE, showWarnings = FALSE)
# INPUT ################################################################################################################
# Microarray ExpressionSet for HNP Differentiataion Experiment, normalized, probe filtered
expressionSet.rds <- here("results/rdata_files/20220829_es_HNP_Differential_Microarray_SST-RMA.rds")

# TargetScan 8.0 Predicted Targets for miR-4707-3p
targetScan8.0.mir4707.targets.txt <- here("data/target_predictions/TargetScan8.0_miR-4707-3p.predicted_targets.txt")

# miRDB v6.0 Target Predictions
miRDB6.0.target.predictions.txt <- here("data/target_predictions/miRDB_v6.0_prediction_result.txt.gz")

# miRTarBase v7 Validated Target Interactions
miRTarBase7.target.interactions.xlsx <- here("data/mirtarbase_v7/hsa_MTI.xlsx")

# miRNAtap predictions, 2 sources
miRNAtap.predictions.2sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_2sources.csv")
# GLOBALS ##############################################################################################################

# Import ExpressionSet #################################################################################################
eSet <- readRDS(expressionSet.rds)

eSet$Condition <- paste(eSet$Timepoint, eSet$Expression, sep = "_")

# PCA ###############

# Uncorrected PCA #################3

pca <- prcomp(t(exprs(eSet)), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(eSet), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# scree plot
# ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
#     geom_point() +
#     geom_line() +
#     labs(x="PC", y="Percent Total Variance",
#          title="PCA: Microarray Expression (log2 Norm.)") +
#     scale_x_continuous(breaks = seq(1:10)) +
#     plotTheme()

df.pca %>%
    ggplot(aes(PC1, PC2, color = cDNA_Yield_ug)) +
    geom_point(size = 1) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme("figure") +
    scale_color_gradientn(colors = paperRedToBlue)

ggsave("~/Desktop/uncorrected_pca_legend.pdf", height = 2, width = 2)



rm(pca, df.pca, percentVar)


# Batch Corrected PCA ##################

corrected.expression <- limma::removeBatchEffect(eSet,
                                                 covariates = eSet$cDNA_Yield_ug,
                                                 design = model.matrix(~eSet$Timepoint + eSet$Expression))

pca <- prcomp(t(corrected.expression), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(eSet), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# scree plot
# ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
#     geom_point() +
#     geom_line() +
#     labs(x="PC", y="Percent Total Variance",
#          title="PCA: Microarray Expression (log2 Norm.)") +
#     scale_x_continuous(breaks = seq(1:10)) +
#     plotTheme()

# df.pca %>%
#     ggplot(aes(PC1, PC2, color = Condition)) +
#     geom_point(size = 3) +
#     #geom_label(hjust="inward", vjust="inward") +
#     labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
#          y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
#     plotTheme() +
#     scale_color_manual(values=cbPalette) +
#     labs(title="Microarray Expression (log2 Norm.)")

df.pca %>%
    ggplot(aes(PC1, PC2, color = outlier)) +
    geom_point(size = 1) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme("figure") +
    scale_color_manual(values=paperCatEight)

ggsave("~/Desktop/corrected_pca_legend.pdf", height = 2, width = 2)


rm(pca, df.pca, percentVar, corrected.expression)


# Outliers Removed ############

filt.eSet <- eSet[,!eSet$outlier]

corrected.expression <- limma::removeBatchEffect(filt.eSet,
                                                 covariates = filt.eSet$cDNA_Yield_ug,
                                                 design = model.matrix(~filt.eSet$Timepoint + filt.eSet$Expression))

pca <- prcomp(t(corrected.expression), center = TRUE, scale. = FALSE)

df.pca <- as_tibble(data.frame(pca$x[,1:20], pData(filt.eSet), stringsAsFactors = FALSE))
percentVar <- pca$sdev^2 / sum(pca$sdev^2) * 100

# scree plot
# ggplot(data.frame(perc_var=percentVar, pc=1:length(percentVar))[1:10,], aes(x=pc, y=perc_var)) +
#     geom_point() +
#     geom_line() +
#     labs(x="PC", y="Percent Total Variance",
#          title="PCA: Microarray Expression (log2 Norm.)") +
#     scale_x_continuous(breaks = seq(1:10)) +
#     plotTheme()

df.pca %>%
    ggplot(aes(PC1, PC2, color = Condition)) +
    geom_point(size = 1) +
    #geom_label(hjust="inward", vjust="inward") +
    labs(x=paste("PC1 (", round(percentVar[1], 1), "%)", sep=""),
         y=paste("PC2 (", round(percentVar[2], 1), "%)", sep="")) +
    plotTheme("figure") +
    scale_color_manual(values=paperCatEight)

ggsave("~/Desktop/outlier_removed_pca_legend.pdf", height = 2, width = 2)





rm(pca, df.pca, percentVar)


# Diff Expression #############

# conduct differential expression primarily based on the 4 condition groups
condition <- factor(filt.eSet$Condition)

# include cDNA_yield in the model matrix
design <- model.matrix(~ 0 + filt.eSet$cDNA_Yield_ug + condition)
colnames(design) <- c("cDNA_yield", levels(condition))

# fit to the expression data
fit <- limma::lmFit(filt.eSet, design)

# contrasts based on expression (contrl v 4707) or differentiation (week1 v week2)
cont.matrix <- makeContrasts(Expression_in_Week1 = Week1_Control-Week1_4707,
                             Expression_in_Week2 = Week2_Control-Week2_4707,
                             Diff_in_Control = Week2_Control-Week1_Control,
                             Diff_in_4707 = Week2_4707-Week1_4707,
                             Interaction = (Week2_Control-Week2_4707)-(Week1_Control-Week1_4707),
                             levels = design)

# fit for contrasts
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# results tables
df.results.expression_in_week1 <- as_tibble(topTable(fit2, coef = "Expression_in_Week1", number=Inf, adjust.method = "BH"))
df.results.expression_in_week2 <- as_tibble(topTable(fit2, coef = "Expression_in_Week2", number=Inf, adjust.method = "BH"))
df.results.diff_in_control <- as_tibble(topTable(fit2, coef = "Diff_in_Control", number=Inf, adjust.method = "BH"))
df.results.diff_in_4707 <- as_tibble(topTable(fit2, coef = "Diff_in_4707", number=Inf, adjust.method = "BH"))
df.results.interaction <- as_tibble(topTable(fit2, coef = "Interaction", number=Inf, adjust.method = "BH"))

# Export DE genes for paper
library(writexl)

write_xlsx(list(expression_in_week1 = df.results.expression_in_week1,
                expression_in_week2 = df.results.expression_in_week2,
                diff_in_control = df.results.diff_in_control,
                diff_in_4707 = df.results.diff_in_4707),
           path = "~/Desktop/supplementaryFile6_diff_expressed_genes_in_phNPCs.xlsx")

# MA Plots ####################

P.ADJ.SIG.THRESH <- 0.05

interestingGenes <- c("NES", "TUBB3", "MKI67", "CCND1", "PAX6", "SOX2", "DCX", "HAUS4")

neuronGenes <- c("TUBB3", "DCX", "NRXN1", "MAPT", "ROR1", "STMN4", "NSG1", "NTRK3")

progGenes <- c("NES", "MKI67", "CCND1", "PAX6", "SOX2", "ASCL1", "CEP135", "CENPF", "PCNA", "TOP2A")

all.genes <- fData(filt.eSet)$SYMBOL

df.results.diff_in_control %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, size=1, color="black") +
    geom_label(data = function(x) subset(x, SYMBOL %in% c(neuronGenes, progGenes)), mapping = aes(label = SYMBOL), size = 2) +
    geom_point(data = function(x) subset(x, SYMBOL %in% c(neuronGenes, progGenes)), size = 2, shape = 21, color = "black") +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    annotate(geom = "text", label = "Up in Week2", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
    annotate(geom = "text", label = "Up in Week1", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6) +
    theme(legend.position = "none") +
    plotTheme("figure") +
    scale_color_manual(values=c("grey70", paperCatEight[1]))

ggsave("~/Desktop/ma_plot_differentiation_control_moreGenes.pdf", height = 2.5, width = 4)


# df.results.diff_in_4707 %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     #filter() %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
#     geom_point() +
#     labs(x="Mean of Normalized Counts",
#          y="Log2 Fold Change",
#          title="Differential Expression: Week1 v Week2 in 4707",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     #scale_x_log10() +
#     theme(legend.position = "bottom") +
#     geom_hline(yintercept = 0, size=1, color="blue") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", "blue"))
#
# df.results.expression_in_week1 %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     #filter() %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
#     geom_point() +
#     labs(x="Mean of Normalized Counts",
#          y="Log2 Fold Change",
#          title="Differential Expression: Control v 4707 in Week1",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     #scale_x_log10() +
#     theme(legend.position = "bottom") +
#     geom_hline(yintercept = 0, size=1, color="blue") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", "blue"))
#
# df.results.expression_in_week2 %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     #filter() %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
#     geom_point() +
#     labs(x="Mean of Normalized Counts",
#          y="Log2 Fold Change",
#          title="Differential Expression: Control v 4707 in Week2",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     #scale_x_log10() +
#     theme(legend.position = "bottom") +
#     geom_hline(yintercept = 0, size=1, color="blue") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", "blue"))
#
# df.results.interaction %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     #filter() %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
#     geom_point() +
#     labs(x="Mean of Normalized Counts",
#          y="Log2 Fold Change",
#          title="Differential Expression: Interaction",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     #scale_x_log10() +
#     theme(legend.position = "bottom") +
#     geom_hline(yintercept = 0, size=1, color="blue") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", "blue"))
#
# df.results.interaction %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = SIG)) +
#     geom_point() +
#     labs(x="Log2 Fold Change",
#          y="-Log10(adj. P)",
#          title="Differential Expression: Interaction",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     #scale_x_log10() +
#     theme(legend.position = "bottom") +
#     #geom_vline(xintercept = 0, size=1, color="black", linetype = "dashed") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", "blue")) +
#     geom_label_repel(data = function(x) subset(x, SIG), mapping = aes(label = SYMBOL), size = 3)

# topGO Analysis ######################3

# topGo data frame file for this analysis
# df.topgo.upreg.week1.control.rds <- paste0(dir.topgo, "topgo_results_upreg_week1_control_elim_fisher_df.rds")
#
# if (file.exists(df.topgo.upreg.week1.control.rds)) {
#     df.topgo.upreg.week1.control.elim.fisher <- readRDS(df.topgo.upreg.week1.control.rds)
#     rm(df.topgo.upreg.week1.control.rds)
# } else {
#     df.results.diff_in_control %>%
#         filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
#         pull(SYMBOL) -> genes.diff.up.control
#
#     geneList <- factor(as.integer(all.genes %in% genes.diff.up.control))
#     names(geneList) <- all.genes
#
#     myGOdata <- new("topGOdata",
#                     ontology = "BP",
#                     description = "my data",
#                     allGenes = geneList,
#                     nodeSize = 10,
#                     annotationFun = annFUN.org,
#                     mapping = "org.Hs.eg.db",
#                     ID = "symbol")
#
#     resultsFisher <- runTest(myGOdata, algorithm = "elim", statistic = "fisher")
#
#     df.topgo.upreg.week1.control.elim.fisher <- GenTable(myGOdata, elim_fisher = resultsFisher, topNodes = 80, orderBy = "elim_fisher")
#
#     saveRDS(df.topgo.upreg.week1.control.elim.fisher, df.topgo.upreg.week1.control.rds)
#
#     rm(genes.diff.up.control, geneList, myGOdata, resultsFisher, df.topgo.upreg.week1.control.rds)
# }
#
# df.topgo.upreg.week1.control.elim.fisher %>%
#     mutate(elim_fisher = as.numeric(elim_fisher),
#            label = paste(GO.ID, Term, sep = "_")) %>%
#     top_n(10, wt = -elim_fisher) %>%
#     ggplot(aes(x = reorder(label, elim_fisher), y = -log10(elim_fisher))) +
#     geom_bar(stat = "identity") +
#     geom_hline(yintercept = -log10(0.05)) +
#     theme(axis.text.x = element_text(angle = 90)) +
#     plotTheme("figure")
#
# ggsave("~/Desktop/go_up_week1.pdf", height = 3, width = 3)

df.topgo.upreg.week1.control.rds <- paste0(dir.topgo, "topgo_results_upreg_week1_control_classic_fisher_df.rds")

if (file.exists(df.topgo.upreg.week1.control.rds)) {
    df.topgo.upreg.week1.control.classic.fisher <- readRDS(df.topgo.upreg.week1.control.rds)
    rm(df.topgo.upreg.week1.control.rds)
} else {
    df.results.diff_in_control %>%
        filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
        pull(SYMBOL) -> genes.diff.up.control

    geneList <- factor(as.integer(all.genes %in% genes.diff.up.control))
    names(geneList) <- all.genes

    myGOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "my data",
                    allGenes = geneList,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "symbol")

    resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

    df.topgo.upreg.week1.control.classic.fisher <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")

    saveRDS(df.topgo.upreg.week1.control.classic.fisher, df.topgo.upreg.week1.control.rds)

    rm(genes.diff.up.control, geneList, myGOdata, resultsFisher, df.topgo.upreg.week1.control.rds)
}

df.topgo.upreg.week1.control.classic.fisher %>%
    mutate(classic_fisher = as.numeric(classic_fisher),
           label = paste(GO.ID, Term, sep = "_")) %>%
    top_n(15, wt = -classic_fisher) %>%
    ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    theme(axis.text.x = element_text(angle = 90)) +
    plotTheme("figure")

ggsave("~/Desktop/go_up_week1.pdf", height = 3, width = 4)

# topGo data frame file for this analysis
# df.topgo.upreg.week2.control.rds <- paste0(dir.topgo, "topgo_results_upreg_week2_control_elim_fisher_df.rds")
#
# if (file.exists(df.topgo.upreg.week2.control.rds)) {
#     df.topgo.upreg.week2.control.elim.fisher <- readRDS(df.topgo.upreg.week2.control.rds)
#     rm(df.topgo.upreg.week2.control.rds)
# } else {
#     df.results.diff_in_control %>%
#         filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
#         pull(SYMBOL) -> genes.diff.down.control
#
#     geneList <- factor(as.integer(all.genes %in% genes.diff.down.control))
#     names(geneList) <- all.genes
#
#     myGOdata <- new("topGOdata",
#                     ontology = "BP",
#                     description = "my data",
#                     allGenes = geneList,
#                     nodeSize = 10,
#                     annotationFun = annFUN.org,
#                     mapping = "org.Hs.eg.db",
#                     ID = "symbol")
#
#     resultsFisher <- runTest(myGOdata, algorithm = "elim", statistic = "fisher")
#
#     df.topgo.upreg.week2.control.elim.fisher <- GenTable(myGOdata, elim_fisher = resultsFisher, topNodes = 80, orderBy = "elim_fisher")
#
#     saveRDS(df.topgo.upreg.week2.control.elim.fisher, df.topgo.upreg.week2.control.rds)
#
#     rm(genes.diff.down.control, myGOdata, resultsFisher, df.topgo.upreg.week2.control.rds)
# }
#
#
# df.topgo.upreg.week2.control.elim.fisher %>%
#     mutate(elim_fisher = as.numeric(elim_fisher)) %>%
#     arrange(elim_fisher) %>%
#     top_n(10, wt = -elim_fisher) %>%
#     ggplot(aes(x = reorder(Term, elim_fisher), y = -log10(elim_fisher))) +
#     geom_bar(stat = "identity") +
#     theme(axis.text.x = element_text(angle = 90, size = 16))
#
#
# ggsave("~/Desktop/go_up_week2.pdf")

df.topgo.upreg.week2.control.rds <- paste0(dir.topgo, "topgo_results_upreg_week2_control_classic_fisher_df.rds")

if (file.exists(df.topgo.upreg.week2.control.rds)) {
    df.topgo.upreg.week2.control.classic.fisher <- readRDS(df.topgo.upreg.week2.control.rds)
    rm(df.topgo.upreg.week2.control.rds)
} else {
    df.results.diff_in_control %>%
        filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
        pull(SYMBOL) -> genes.diff.down.control

    geneList <- factor(as.integer(all.genes %in% genes.diff.down.control))
    names(geneList) <- all.genes

    myGOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "my data",
                    allGenes = geneList,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "symbol")

    resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

    df.topgo.upreg.week2.control.classic.fisher <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")

    saveRDS(df.topgo.upreg.week2.control.classic.fisher, df.topgo.upreg.week2.control.rds)

    rm(genes.diff.down.control, myGOdata, resultsFisher, df.topgo.upreg.week2.control.rds)
}


df.topgo.upreg.week2.control.classic.fisher %>%
    mutate(classic_fisher = as.numeric(classic_fisher),
           label = paste(GO.ID, Term, sep = "_")) %>%
    top_n(15, wt = -classic_fisher) %>%
    ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    theme(axis.text.x = element_text(angle = 90)) +
    plotTheme("figure")

ggsave("~/Desktop/go_up_week2.pdf", height = 3, width = 4)




# # topGo data frame file for this analysis
# df.topgo.upreg.control.week2.rds <- paste0(dir.topgo, "topgo_results_upreg_control_week2_classic_fisher_df.rds")
#
# if (file.exists(df.topgo.upreg.control.week2.rds)) {
#     df.topgo.upreg.control.week2.classic.fisher <- readRDS(df.topgo.upreg.control.week2.rds)
#     rm(df.topgo.upreg.control.week2.rds)
# } else {
#     df.results.expression_in_week2 %>%
#         filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
#         pull(SYMBOL) -> genes.exp.up.week2
#
#     geneList <- factor(as.integer(all.genes %in% genes.exp.up.week2))
#     names(geneList) <- all.genes
#
#     myGOdata <- new("topGOdata",
#                     ontology = "BP",
#                     description = "my data",
#                     allGenes = geneList,
#                     nodeSize = 10,
#                     annotationFun = annFUN.org,
#                     mapping = "org.Hs.eg.db",
#                     ID = "symbol")
#
#     resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")
#
#     df.topgo.upreg.control.week2.classic.fisher <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")
#
#     saveRDS(df.topgo.upreg.control.week2.classic.fisher, df.topgo.upreg.control.week2.rds)
#
#     rm(genes.exp.up.week2, geneList, myGOdata, resultsFisher, df.topgo.upreg.control.week2.rds)
# }
#
# df.topgo.upreg.control.week2.classic.fisher %>%
#     mutate(classic_fisher = as.numeric(classic_fisher ),
#            label = paste(GO.ID, Term, sep = "_")) %>%
#     arrange(classic_fisher) %>%
#     top_n(20, wt = -classic_fisher) %>%
#     ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
#     geom_bar(stat = "identity") +
#     geom_hline(yintercept = -log10(0.05)) +
#     plotTheme("figure") +
#     theme(axis.text.x = element_text(angle = 90, size = 6))
#
#
# ggsave("~/Desktop/go_up_control.pdf", height = 4, width = 3)
#
#
# # topGo data frame file for this analysis
# df.topgo.upreg.4707.week2.rds <- paste0(dir.topgo, "topgo_results_upreg_4707_week2_classic_fisher_df.rds")
#
# if (file.exists(df.topgo.upreg.4707.week2.rds)) {
#     df.topgo.upreg.4707.week2.classic.fisher <- readRDS(df.topgo.upreg.4707.week2.rds)
#     rm(df.topgo.upreg.4707.week2.rds)
# } else {
#     df.results.expression_in_week2 %>%
#         filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
#         pull(SYMBOL) -> genes.exp.down.week2
#
#     geneList <- factor(as.integer(all.genes %in% genes.exp.down.week2))
#     names(geneList) <- all.genes
#
#     myGOdata <- new("topGOdata",
#                     ontology = "BP",
#                     description = "my data",
#                     allGenes = geneList,
#                     nodeSize = 10,
#                     annotationFun = annFUN.org,
#                     mapping = "org.Hs.eg.db",
#                     ID = "symbol")
#
#     resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")
#
#     df.topgo.upreg.4707.week2.classic.fisher <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")
#
#     saveRDS(df.topgo.upreg.4707.week2.classic.fisher, df.topgo.upreg.4707.week2.rds)
#
#     rm(genes.exp.down.week2, geneList, myGOdata, resultsFisher, df.topgo.upreg.4707.week2.rds)
# }
#
# df.topgo.upreg.4707.week2.classic.fisher %>%
#     mutate(classic_fisher = as.numeric(classic_fisher),
#            label = paste(GO.ID, Term, sep = "_")) %>%
#     arrange(classic_fisher) %>%
#     top_n(20, wt = -classic_fisher) %>%
#     ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
#     geom_bar(stat = "identity") +
#     geom_hline(yintercept = -log10(0.05)) +
#     plotTheme("figure") +
#     theme(axis.text.x = element_text(angle = 90, size = 6))
#
#
# ggsave("~/Desktop/go_up_4707.pdf", width = 3, height = 4)

# topGo data frame file for this analysis
df.topgo.week2.expression.rds <- paste0(dir.topgo, "topgo_results_week2_expression_classic_fisher_df.rds")

if (file.exists(df.topgo.week2.expression.rds)) {
    df.topgo.week2.expression.classic.fisher <- readRDS(df.topgo.week2.expression.rds)
    rm(df.topgo.week2.expression.rds)
} else {
    df.results.expression_in_week2 %>%
        filter(adj.P.Val < P.ADJ.SIG.THRESH) %>%
        pull(SYMBOL) -> genes.week2.expression

    geneList <- factor(as.integer(all.genes %in% genes.week2.expression))
    names(geneList) <- all.genes

    myGOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "my data",
                    allGenes = geneList,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "symbol")

    resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

    df.topgo.week2.expression.classic.fisher <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")

    saveRDS(df.topgo.week2.expression.classic.fisher, df.topgo.week2.expression.rds)

    rm(genes.week2.expression, geneList, myGOdata, resultsFisher, df.topgo.week2.expression.rds)
}

df.topgo.week2.expression.classic.fisher %>%
    mutate(classic_fisher = as.numeric(classic_fisher),
           label = paste(GO.ID, Term, sep = "_")) %>%
    top_n(15, wt = -classic_fisher) %>%
    ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    theme(axis.text.x = element_text(angle = 90)) +
    plotTheme("figure")

# ggsave("~/Desktop/go_up_week2.pdf", height = 3, width = 4)

# topGo data frame file for this analysis
df.topgo.week2.expression.rds <- paste0(dir.topgo, "topgo_results_week2_expression_elim_fisher_df.rds")

if (file.exists(df.topgo.week2.expression.rds)) {
    df.topgo.week2.expression.elim.fisher <- readRDS(df.topgo.week2.expression.rds)
    rm(df.topgo.week2.expression.rds)
} else {
    df.results.expression_in_week2 %>%
        filter(adj.P.Val < P.ADJ.SIG.THRESH) %>%
        pull(SYMBOL) -> genes.week2.expression

    geneList <- factor(as.integer(all.genes %in% genes.week2.expression))
    names(geneList) <- all.genes

    myGOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "my data",
                    allGenes = geneList,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "symbol")

    resultsFisher <- runTest(myGOdata, algorithm = "elim", statistic = "fisher")

    df.topgo.week2.expression.elim.fisher <- GenTable(myGOdata, elim_fisher = resultsFisher, topNodes = 80, orderBy = "elim_fisher")

    saveRDS(df.topgo.week2.expression.elim.fisher, df.topgo.week2.expression.rds)

    rm(genes.week2.expression, geneList, myGOdata, resultsFisher, df.topgo.week2.expression.rds)
}

df.topgo.week2.expression.elim.fisher %>%
    mutate(elim_fisher = as.numeric(elim_fisher),
           label = paste(GO.ID, Term, sep = "_")) %>%
    top_n(15, wt = -elim_fisher) %>%
    ggplot(aes(x = reorder(label, elim_fisher), y = -log10(elim_fisher))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    theme(axis.text.x = element_text(angle = 90)) +
    plotTheme("figure")

ggsave("~/Desktop/go_week2_expression.pdf", height = 3, width = 4)



# Targets ##########

df.targetScan.predictions <- read_tsv(targetScan8.0.mir4707.targets.txt)

df.mirdb.predictions <- read_tsv(miRDB6.0.target.predictions.txt, col_names = c("miRNA_Name", "GenBank", "Target_Score", "ENTREZID"))
df.mirdb.predictions %<>%
    filter(miRNA_Name == "hsa-miR-4707-3p") %>%
    filter(!duplicated(ENTREZID)) %>%
    mutate(ENTREZID = as.character(ENTREZID))

anno <- AnnotationDbi::select(org.Hs.eg.db,
                              keys = df.mirdb.predictions$ENTREZID,
                              keytype = "ENTREZID",
                              columns = c("SYMBOL"))

df.mirdb.predictions %<>%
    left_join(anno, by = "ENTREZID")

df.mirtarbase.interactions <- read_xlsx(miRTarBase7.target.interactions.xlsx)
df.mirtarbase.interactions %<>%
    filter(miRNA == "hsa-miR-4707-3p") %>%
    filter(!duplicated(`Target Gene`))



### Differentially Expressed Targets


# week 1, control vs 4707, genes that are up in control (down in 4707)
df.results.expression_in_week1 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
    pull(SYMBOL) -> genes.exp.up.in.control.week1



# week 2, control vs 4707, genes that are up in control (down in 4707)
df.results.expression_in_week2 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
    pull(SYMBOL) -> genes.exp.up.in.control.week2


# genes that down-regulated by 4707 expression in BOTH week1 and week2
genes.downregulated.by.4707 <- genes.exp.up.in.control.week2[genes.exp.up.in.control.week2 %in% genes.exp.up.in.control.week1]



# topGo data frame file for this analysis
df.topgo.downregulated.by4707.rds <- paste0(dir.topgo, "topgo_results_downreg_by_4707_classic_fisher_df.rds")

if (file.exists(df.topgo.downregulated.by4707.rds)) {
    topgo.genes.downregulated.by4707 <- readRDS(df.topgo.downregulated.by4707.rds)
    rm(df.topgo.downregulated.by4707.rds)
} else {
    geneList <- factor(as.integer(all.genes %in% genes.downregulated.by.4707))
    names(geneList) <- all.genes

    myGOdata <- new("topGOdata",
                    ontology = "BP",
                    description = "my data",
                    allGenes = geneList,
                    nodeSize = 10,
                    annotationFun = annFUN.org,
                    mapping = "org.Hs.eg.db",
                    ID = "symbol")

    resultsFisher <- runTest(myGOdata, algorithm = "elim", statistic = "fisher")

    topgo.genes.downregulated.by4707 <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")

    saveRDS(topgo.genes.downregulated.by4707, df.topgo.downregulated.by4707.rds)

    rm(geneList, myGOdata, resultsFisher, df.topgo.downregulated.by4707.rds)
}

topgo.genes.downregulated.by4707 %>%
    mutate(classic_fisher = as.numeric(classic_fisher),
           label = paste(GO.ID, Term, sep = "_")) %>%
    top_n(15, wt = -classic_fisher) %>%
    ggplot(aes(x = reorder(label, classic_fisher), y = -log10(classic_fisher))) +
    geom_bar(stat = "identity") +
    geom_hline(yintercept = -log10(0.05)) +
    theme(axis.text.x = element_text(angle = 90))









#### Genes Downregulated AND Targeted by 4707


genes.downregulated.and.targeted.by.4707 <- genes.downregulated.by.4707[genes.downregulated.by.4707 %in% df.targetScan.predictions$`Target gene`]

genes.downregulated.by.4707[genes.downregulated.by.4707 %in% df.mirtarbase.interactions$`Target Gene`]

genes.downregulated.and.targeted.by.4707

# df.results.expression_in_week1 %>%
#     mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
#     arrange(SIG) %>%
#     ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
#     geom_point() +
#     geom_hline(yintercept = 0, size=1, color="black") +
#     geom_label_repel(data = function(x) subset(x, SYMBOL %in% genes.downregulated.and.targeted.by.4707), mapping = aes(label = SYMBOL), size = 3, max.overlaps = 30, color = "red") +
#     geom_point(data = function(x) subset(x, SYMBOL %in% genes.downregulated.and.targeted.by.4707), size = 3, shape = 1, color = "red") +
#     labs(x="Mean of Normalized Counts",
#          y="Log2 Fold Change",
#          title="Differential Expression: Control vs mir4707 in Week 1",
#          color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
#     annotate(geom = "text", label = "Up in Control", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
#     annotate(geom = "text", label = "Up in mir4707", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6) +
#     theme(legend.position = "bottom") +
#     plotTheme() +
#     scale_color_manual(values=c("grey70", cbPalette[4]))

df.results.expression_in_week2 %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point(size = 1) +
    geom_hline(yintercept = 0, size=1, color="black") +
    geom_text(data = function(x) subset(x, SYMBOL %in% genes.downregulated.and.targeted.by.4707), mapping = aes(label = SYMBOL), size = 1, color = "black") +
    geom_point(data = function(x) subset(x, SYMBOL %in% genes.downregulated.and.targeted.by.4707), size = 2, shape = 21, color = "black") +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    annotate(geom = "text", label = "Up in Control", x = Inf, y = Inf, hjust = 1.1, vjust = 1.5) +
    annotate(geom = "text", label = "Up in mir4707", x = Inf, y = -Inf, hjust = 1.1, vjust = -0.6) +
    theme(legend.position = "none") +
    plotTheme("figure") +
    scale_color_manual(values=c("grey70", paperCatEight[1])) +
    scale_y_continuous(limits = c(-2.8,2.8))

ggsave("~/Desktop/ma_plot_week2_expression.pdf", height = 1.8, width = 4)

