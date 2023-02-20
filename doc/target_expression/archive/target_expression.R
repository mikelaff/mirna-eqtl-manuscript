# diff expression of mirna targets
# mirnas that are diff expressed

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(topGO)
library(miRNAtap)
library(org.Hs.eg.db)
library(AnnotationHub)
library(gridExtra)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_expression/pdfs/")
# write plots
write_plots <- FALSE

pdf("~/Desktop/top_mirna_targets_by_prediction.pdf", height = 8, width = 10)

# INPUT FILES #####################################################################################
# combined deseq2 shrunken results of gene expression in gzcp or gw as a dataframe with labels for 
# significance and category (neuro vs maturation)
gene_expression_combined_df_rds <- here("doc/diff_expression/rdata/combined_gene_expression_deseq2_shrunken_results_df.rds")
# combined deseq2 shrunken results of mirna expression in gzcp or gw as a dataframe with labels for 
# significance and category (neuro vs maturation)
mirna_expression_combined_df_rds <- here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_df.rds")

# Import Data #####################################################################################
gene_df <- readRDS(gene_expression_combined_df_rds)
mirna_df <- readRDS(mirna_expression_combined_df_rds)

# Annotation Hub ##################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]

# Quick Plots #####################################################################################
labelColor <- "black"
circColor <- "black"
# mirna diff expression plot by categories
mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

# gene diff expression plot by categories
gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

# miRNA Targets in Graph #######################################################################
# Neurogenesis: up in GZ
mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", baseMean.gzcp > 1000) %>%
  dplyr::select(mirna, baseMean.gzcp, log2FoldChange.gzcp, padj.gzcp, sig.both, category) %>%
  #dplyr::top_n(n=10, wt=log2FoldChange.gzcp) %>%
  dplyr::arrange(-log2FoldChange.gzcp) -> top_mirnas

p <- list()

for (i in 1:8) {
  # truncate name to not include '/' suffixes or 'hsa-' prefix
  mir <- strsplit(strsplit(top_mirnas$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(mir)
  # get predicted targets
  predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
  # check for no targets
  if (length(predictions) == 0) {
    next
  }
  # select top entrez ids
  predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]
  # get ensg and names for entrez ids
  predicted_rowdata <- AnnotationDbi::select(orgdb, 
                                             keys = predicted_ids_entrez, 
                                             columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                             keytype = "ENTREZID")
  gene_df %>%
    filter(category != "none") %>%
    ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
    geom_point(size=.5) +
    labs(title = paste(mir, "counts:", round(top_mirnas$baseMean.gzcp[i],2), "l2fc:", round(top_mirnas$log2FoldChange.gzcp[i],2))) +
    plotTheme +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=6)) +
    scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
    geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL, category != "none"),
               aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp),
               color="black") -> p[[i]]
}

do.call(grid.arrange,p)

# Neurogenesis: up in CP
mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", baseMean.gzcp > 1000) %>%
  dplyr::select(mirna, baseMean.gzcp, log2FoldChange.gzcp, padj.gzcp, sig.both, category) %>%
  #dplyr::top_n(n=10, wt=log2FoldChange.gzcp) %>%
  dplyr::arrange(log2FoldChange.gzcp) -> top_mirnas

p <- list()

for (i in 1:9) {
  # truncate name to not include '/' suffixes or 'hsa-' prefix
  mir <- strsplit(strsplit(top_mirnas$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(mir)
  # get predicted targets
  predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
  # check for no targets
  if (length(predictions) == 0) {
    next
  }
  # select top entrez ids
  predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]
  # get ensg and names for entrez ids
  predicted_rowdata <- AnnotationDbi::select(orgdb, 
                                             keys = predicted_ids_entrez, 
                                             columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                             keytype = "ENTREZID")
  gene_df %>%
    filter(category != "none") %>%
    ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
    geom_point(size=.5) +
    labs(title = paste(mir, "counts:", round(top_mirnas$baseMean.gzcp[i],2), "l2fc:", round(top_mirnas$log2FoldChange.gzcp[i],2))) +
    plotTheme +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=6)) +
    scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
    geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL, category != "none"),
               aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp),
               color="black") -> p[[i]]
}

do.call(grid.arrange,p)

# Maturation: up early
mirna_df %>%
  dplyr::filter(category == "maturation_early", baseMean.gw > 1000) %>%
  dplyr::select(mirna, baseMean.gw, log2FoldChange.gw, padj.gw, sig.both, category) %>%
  #dplyr::top_n(n=10, wt=log2FoldChange.gw) %>%
  dplyr::arrange(log2FoldChange.gw) -> top_mirnas

p <- list()

for (i in 1:5) {
  # truncate name to not include '/' suffixes or 'hsa-' prefix
  mir <- strsplit(strsplit(top_mirnas$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(mir)
  # get predicted targets
  predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
  # check for no targets
  if (length(predictions) == 0) {
    next
  }
  # select top entrez ids
  predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]
  # get ensg and names for entrez ids
  predicted_rowdata <- AnnotationDbi::select(orgdb, 
                                             keys = predicted_ids_entrez, 
                                             columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                             keytype = "ENTREZID")
  gene_df %>%
    filter(category != "none") %>%
    ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
    geom_point(size=.5) +
    labs(title = paste(mir, "counts:", round(top_mirnas$baseMean.gw[i],2), "l2fc:", round(top_mirnas$log2FoldChange.gw[i],2))) +
    plotTheme +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=6)) +
    scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
    geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL, category != "none"),
               aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp),
               color="black") -> p[[i]]
}

do.call(grid.arrange,p)

# Maturation: up late
mirna_df %>%
  dplyr::filter(category == "maturation_late", baseMean.gw > 1000) %>%
  dplyr::select(mirna, baseMean.gw, log2FoldChange.gw, padj.gw, sig.both, category) %>%
  #dplyr::top_n(n=10, wt=log2FoldChange.gw) %>%
  dplyr::arrange(-log2FoldChange.gw) -> top_mirnas

p <- list()

for (i in 1:9) {
  # truncate name to not include '/' suffixes or 'hsa-' prefix
  mir <- strsplit(strsplit(top_mirnas$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(mir)
  # get predicted targets
  predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
  # check for no targets
  if (length(predictions) == 0) {
    next
  }
  # select top entrez ids
  predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]
  # get ensg and names for entrez ids
  predicted_rowdata <- AnnotationDbi::select(orgdb, 
                                             keys = predicted_ids_entrez, 
                                             columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                             keytype = "ENTREZID")
  gene_df %>%
    filter(category != "none") %>%
    ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
    geom_point(size=.5) +
    labs(title = paste(mir, "counts:", round(top_mirnas$baseMean.gw[i],2), "l2fc:", round(top_mirnas$log2FoldChange.gw[i],2))) +
    plotTheme +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(size=6)) +
    scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
    geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL, category != "none"),
               aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp),
               color="black") -> p[[i]]
}

do.call(grid.arrange,p)

#

dev.off()





























#####################################
#
mir <- "miR-92b-3p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")

#
mir <- "miR-99a-5p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")

#
mir <- "miR-181a-2-3p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")

#
mir <- "miR-25-3p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")

#
mir <- "miR-16-5p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")

#
mir <- "let-7a-5p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)
predicted_ids_entrez <- dimnames(predictions)[[1]][1:10]

predicted_rowdata <- AnnotationDbi::select(orgdb, keys = predicted_ids_entrez, columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), keytype = "ENTREZID")

gene_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=filter(gene_df, ensg %in% predicted_rowdata$ENSEMBL), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), color="black")










