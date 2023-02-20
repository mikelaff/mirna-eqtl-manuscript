
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

pData(eSet)

filt.eSet <- eSet[,!eSet$outlier]

pData(filt.eSet)

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



# MA Plots ####################

P.ADJ.SIG.THRESH <- 0.05

df.results.diff_in_control %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    #filter() %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         title="Differential Expression: Week1 v Week2 in Control",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue"))

df.results.diff_in_4707 %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    #filter() %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         title="Differential Expression: Week1 v Week2 in 4707",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue"))

df.results.expression_in_week1 %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    #filter() %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         title="Differential Expression: Control v 4707 in Week1",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue"))

df.results.expression_in_week2 %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    #filter() %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         title="Differential Expression: Control v 4707 in Week2",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue"))

df.results.interaction %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    #filter() %>%
    arrange(SIG) %>%
    ggplot(aes(x = AveExpr, y = logFC, color = SIG)) +
    geom_point() +
    labs(x="Mean of Normalized Counts",
         y="Log2 Fold Change",
         title="Differential Expression: Interaction",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    geom_hline(yintercept = 0, size=1, color="blue") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue"))

df.results.interaction %>%
    mutate(SIG = adj.P.Val < P.ADJ.SIG.THRESH) %>%
    arrange(SIG) %>%
    ggplot(aes(x = logFC, y = -log10(adj.P.Val), color = SIG)) +
    geom_point() +
    labs(x="Log2 Fold Change",
         y="-Log10(adj. P)",
         title="Differential Expression: Interaction",
         color=paste("p.adj <", P.ADJ.SIG.THRESH, sep="")) +
    #scale_x_log10() +
    theme(legend.position = "bottom") +
    #geom_vline(xintercept = 0, size=1, color="black", linetype = "dashed") +
    plotTheme() +
    scale_color_manual(values=c("grey70", "blue")) +
    geom_label_repel(data = function(x) subset(x, SIG), mapping = aes(label = SYMBOL), size = 3)

# topGO Analysis ######################3

# all expressed genes
all.genes <- fData(filt.eSet)$SYMBOL

# diff up (upregulated in progenitors)
df.results.diff_in_control %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
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

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")

# diff down (upregulated in neurons)
df.results.diff_in_control %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
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

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")

# # background genes by expression
# df.results.diff_in_control %>%
#     filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
#     pull(PROBEID) -> probes.diff.up.control
#
# back_genes_idx <- genefilter::genefinder(filt.eSet, as.character(probes.diff.up.control), method = "manhattan", scale = "none")
#
# back_genes_idx <- sapply(back_genes_idx, function(x)x$indices)
#
# back_genes <- featureNames(filt.eSet)[back_genes_idx]
# back_genes <- setdiff(back_genes, probes.diff.up.control)
#
# intersect(back_genes, probes.diff.up.control)
# length(back_genes)
#
# gene_ids <- df.results.diff_in_control$PROBEID
# in_universe <- gene_ids %in% c(probes.diff.up.control, back_genes)
# in_selection <- gene_ids %in% probes.diff.up.control
#
# all_genes <- in_selection[in_universe]
# all_genes <- factor(as.integer(in_selection[in_universe]))
# names(all_genes) <- gene_ids[in_universe]


# miRNA Targets ############


# multimir_dbInfoVersions()
#
# df.predicted <- get_multimir(org = "hsa",
#                              mirna = "hsa-miR-4707-3p",
#                              table = "predicted")@data
#
# df.validated <- get_multimir(org = "hsa",
#                              mirna = "hsa-miR-4707-3p",
#                              table = "validated")@data
#
# df.pred <- getPredictedTargets(mirna = "miR-4707-3p", species = "hsa")


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


table(df.mirdb.predictions$SYMBOL %in% df.targetScan.predictions$`Target gene`)

table(df.mirdb.predictions$SYMBOL %in% df.mirtarbase.interactions$`Target Gene`)

table(df.mirtarbase.interactions$`Target Gene` %in% df.targetScan.predictions$`Target gene`)


# Something #####

# week 1, control vs 4707, genes that are up in control (down in 4707)
df.results.expression_in_week1 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
    pull(SYMBOL) -> genes.exp.up.in.control.week1

table(genes.exp.up.in.control.week1 %in% df.mirdb.predictions$SYMBOL)

table(genes.exp.up.in.control.week1 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.exp.up.in.control.week1 %in% df.targetScan.predictions$`Target gene`)

# week 2, control vs 4707, genes that are up in control (down in 4707)
df.results.expression_in_week2 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC > 0) %>%
    pull(SYMBOL) -> genes.exp.up.in.control.week2

table(genes.exp.up.in.control.week2 %in% df.mirdb.predictions$SYMBOL)

table(genes.exp.up.in.control.week2 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.exp.up.in.control.week2 %in% df.targetScan.predictions$`Target gene`)

# genes that down-regulated by 4707 expression in BOTH week1 and week2
genes.downregulated.by.4707 <- genes.exp.up.in.control.week2[genes.exp.up.in.control.week2 %in% genes.exp.up.in.control.week1]

table(genes.downregulated.by.4707 %in% df.mirdb.predictions$SYMBOL)

table(genes.downregulated.by.4707 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.downregulated.by.4707 %in% df.targetScan.predictions$`Target gene`)


# combine expression data into one dataframe

df.results.diff_in_control %<>%
    mutate(Sig_at_0.05 = adj.P.Val < 0.05,
           Sig_at_0.1 = adj.P.Val < 0.1)

df.results.diff_in_4707 %<>%
    mutate(Sig_at_0.05 = adj.P.Val < 0.05,
           Sig_at_0.1 = adj.P.Val < 0.1)

df.results.expression_in_week1 %<>%
    mutate(Sig_at_0.05 = adj.P.Val < 0.05,
           Sig_at_0.1 = adj.P.Val < 0.1)

df.results.expression_in_week2 %<>%
    mutate(Sig_at_0.05 = adj.P.Val < 0.05,
           Sig_at_0.1 = adj.P.Val < 0.1)

df.results.diff <- left_join(df.results.diff_in_control, df.results.diff_in_4707,
                             by = c("PROBEID", "SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "duplicate.annotation"),
                             suffix = c(".diff_in_control", ".diff_in_4707"))


df.results.expression <- left_join(df.results.expression_in_week1, df.results.expression_in_week2,
                                   by = c("PROBEID", "SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "duplicate.annotation"),
                                   suffix = c(".expression_in_week1", ".expression_in_week2"))

df.results <- left_join(df.results.diff, df.results.expression,
                                   by = c("PROBEID", "SYMBOL", "GENENAME", "ENSEMBL", "ENTREZID", "duplicate.annotation"))


df.results %<>%
    mutate(miRDB_predicted_target = SYMBOL %in% df.mirdb.predictions$SYMBOL,
           TargetScan_predicted_target = SYMBOL %in% df.targetScan.predictions$`Target gene`,
           miRTarBase_validated_target = SYMBOL %in% df.mirtarbase.interactions$`Target Gene`)

df.results %>%
    filter(Sig_at_0.05.expression_in_week1 & Sig_at_0.05.expression_in_week2) %>%
    filter(logFC.expression_in_week1 > 0 & logFC.expression_in_week2 > 0) -> df.sig

mat.results <- as.matrix(dplyr::select(df.sig, ends_with("_target")))
rownames(mat.results) <- df.sig$SYMBOL

m <- make_comb_mat(mat.results)

UpSet(m)


# topgo on downregulated genes


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

resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

topgo.genes.downregulated.by4707 <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")


geneList <- factor(as.integer(all.genes %in% genes.downregulated.by.4707 & all.genes %in% df.targetScan.predictions$`Target gene`))
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

topgo.genes.downregulated.and.targeted.by4707 <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")


df.results %<>%
    mutate(down_and_target = SYMBOL %in% genes.downregulated.by.4707[genes.downregulated.by.4707 %in% df.targetScan.predictions$`Target gene`])


df.2sources <- read_csv(miRNAtap.predictions.2sources.csv)
df.2sources %<>%
    filter(miRNA == "hsa-miR-4707-5p")

# week 1, control vs 4707, genes that are up in 4707 (down in control)
df.results.expression_in_week1 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
    pull(SYMBOL) -> genes.exp.down.in.control.week1

table(genes.exp.down.in.control.week1 %in% df.mirdb.predictions$SYMBOL)

table(genes.exp.down.in.control.week1 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.exp.down.in.control.week1 %in% df.targetScan.predictions$`Target gene`)

# week 2, control vs 4707, genes that are up in 4707 (down in control)
df.results.expression_in_week2 %>%
    filter(adj.P.Val < P.ADJ.SIG.THRESH & logFC < 0) %>%
    pull(SYMBOL) -> genes.exp.down.in.control.week2

table(genes.exp.down.in.control.week2 %in% df.mirdb.predictions$SYMBOL)

table(genes.exp.down.in.control.week2 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.exp.down.in.control.week2 %in% df.targetScan.predictions$`Target gene`)

# genes that up-regulated by 4707 expression in BOTH week1 and week2
genes.upregulated.by.4707 <- genes.exp.down.in.control.week2[genes.exp.down.in.control.week2 %in% genes.exp.down.in.control.week1]

table(genes.upregulated.by.4707 %in% df.mirdb.predictions$SYMBOL)

table(genes.upregulated.by.4707 %in% df.mirtarbase.interactions$`Target Gene`)

table(genes.upregulated.by.4707 %in% df.targetScan.predictions$`Target gene`)

geneList <- factor(as.integer(all.genes %in% genes.upregulated.by.4707))
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

topgo.genes.upregulated.by.4707 <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 80, orderBy = "classic_fisher")


genes.upregulated.by.4707[genes.upregulated.by.4707 %in% df.targetScan.predictions$`Target gene`]


df.results %>%
    dplyr::filter(SYMBOL %in% genes.upregulated.by.4707 | SYMBOL %in% genes.downregulated.by.4707) -> df.tmp


sum(duplicated(df.tmp$SYMBOL))


# Interaction plots?

exprs(filt.eSet)

tmp <- as.data.frame(t(exprs(filt.eSet[match("IGFBP3", fData(filt.eSet)$SYMBOL),])))
tmp$CEL <- rownames(tmp)
colnames(tmp) <- c("IGFBP3_expression", "CEL")


df.samples <- as_tibble(pData(filt.eSet))

df.samples %<>%
    left_join(tmp, by = "CEL")

df.samples %>%
    ggplot(aes(x = Timepoint, y = IGFBP3_expression, color = Condition)) +
    geom_point()


tmp <- as.data.frame(t(exprs(filt.eSet[match("PDGFRA", fData(filt.eSet)$SYMBOL),])))
tmp$CEL <- rownames(tmp)
colnames(tmp) <- c("PDGFRA_expression", "CEL")


df.samples <- as_tibble(pData(filt.eSet))

df.samples %<>%
    left_join(tmp, by = "CEL")

df.samples %>%
    ggplot(aes(x = Timepoint, y = PDGFRA_expression, color = Expression)) +
    geom_point() +
    geom_boxplot() +
    geom_smooth(method = lm)





