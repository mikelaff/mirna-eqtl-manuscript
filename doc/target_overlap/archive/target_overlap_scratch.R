# look at miRNA targets and their overlap with differentially expressed genes
# scratch for fisher exact and chi squared tests on contingency tables

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(AnnotationHub)
library(miRNAtap)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_overlap/pdfs/")
# write plots
write_plots <- FALSE

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
rm(ah, orgs)


# scratch

mir <- "miR-92b-3p"

# get predicted targets
predictions <- data.frame(getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2))
# number of sources per predicted gene
predictions$sources <- apply(predictions[,1:5], 1, function(x) sum(!is.na(x)))
# get gene ensg and symbol
predictions$entrezid <- rownames(predictions)
tmp <- AnnotationDbi::select(orgdb, 
                             keys = predictions$entrezid, 
                             columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                             keytype = "ENTREZID")
predictions <- left_join(predictions, tmp, by = c("entrezid" = "ENTREZID"))

# Fisher's Exact Test ######

# all genes in analysis are in either gz>cp or cp>gz
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp")]

# the targets that overlap gz>cp or cp>gz
targets_gz_greater_than_cp <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_gz_greater_than_cp)]
targets_cp_greater_than_gz <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_cp_greater_than_gz)]

# cont. table: genes are either targets or not targets, and genes are either gz>cp or cp>gz
tab <- matrix(c(length(genes_gz_greater_than_cp)-length(targets_gz_greater_than_cp), length(genes_cp_greater_than_gz)-length(targets_cp_greater_than_gz),
                length(targets_gz_greater_than_cp), length(targets_cp_greater_than_gz)), byrow = TRUE, nrow = 2)

# labels rows and columns for this table
colnames(tab) <- c("gz>cp", "cp>gz")
rownames(tab) <- c("not_targets", "targets")

tab

ft <- fisher.test(tab)
ft
ft.2 <- fisher.test(tab, alternative = "two.sided")
ft.2

# same data in different order
tab2 <- matrix(c(length(targets_gz_greater_than_cp), length(targets_cp_greater_than_gz),
                length(genes_gz_greater_than_cp)-length(targets_gz_greater_than_cp), length(genes_cp_greater_than_gz)-length(targets_cp_greater_than_gz)), byrow = TRUE, nrow = 2)

colnames(tab2) <- c("gz>cp", "cp>gz")
rownames(tab2) <- c("targets", "not_targets")

tab2

ft2 <- fisher.test(tab2)
ft2

# same data in different order again
tab3 <- t(tab)
tab3

ft3 <- fisher.test(tab3)
ft3

# Chi Squared Test ###################

# all genes are in one of 4 groups: gz>cp, cp>gz, mat_early, or mat_late
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp")]
genes_maturation_early <- gene_df$ensg[which(gene_df$category == "maturation_early")]
genes_maturation_late <- gene_df$ensg[which(gene_df$category == "maturation_late")]

genes_none <- gene_df$ensg[which(gene_df$category == "none")]

# targets that overlap into one of 4 groups
targets_gz_greater_than_cp <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_gz_greater_than_cp)]
targets_cp_greater_than_gz <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_cp_greater_than_gz)]
targets_maturation_early <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_maturation_early)]
targets_maturation_late <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_maturation_late)]

targets_none <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_none)]

ctab <- matrix(nrow = 2, ncol = 4)
colnames(ctab) <- c("gz>cp", "cp>gz", "mat_early", "mat_late")
rownames(ctab) <- c("not_targets", "targets")
ctab

ctab["not_targets", "gz>cp"] <- length(genes_gz_greater_than_cp)-length(targets_gz_greater_than_cp)
ctab["not_targets", "cp>gz"] <- length(genes_cp_greater_than_gz)-length(targets_cp_greater_than_gz)
ctab["not_targets", "mat_early"] <- length(genes_maturation_early)-length(targets_maturation_early)
ctab["not_targets", "mat_late"] <- length(genes_maturation_late)-length(targets_maturation_late)

ctab["targets", "gz>cp"] <- length(targets_gz_greater_than_cp)
ctab["targets", "cp>gz"] <- length(targets_cp_greater_than_gz)
ctab["targets", "mat_early"] <- length(targets_maturation_early)
ctab["targets", "mat_late"] <- length(targets_maturation_late)

ctab

ct <- chisq.test(ctab)
ct

ft.ctab <- fisher.test(ctab, simulate.p.value = TRUE)
ft.ctab

ctab2 <- matrix(c(ctab[1,1], ctab[1,2]+ctab[1,3]+ctab[1,4], ctab[2,1], ctab[2,2]+ctab[2,3]+ctab[2,4]), nrow = 2, byrow = TRUE)
rownames(ctab2) <- c("not_targets", "targets")
colnames(ctab2) <- c("gz>cp", "other")
ctab2
fisher.test(ctab2)

ctab3 <- matrix(c(ctab[1,2], ctab[1,1]+ctab[1,3]+ctab[1,4], ctab[2,2], ctab[2,1]+ctab[2,3]+ctab[2,4]), nrow = 2, byrow = TRUE)
rownames(ctab3) <- c("not_targets", "targets")
colnames(ctab3) <- c("cp>gz", "other")
ctab3
fisher.test(ctab3)

ctab4 <- matrix(c(ctab[1,3], ctab[1,1]+ctab[1,2]+ctab[1,4], ctab[2,3], ctab[2,1]+ctab[2,2]+ctab[2,4]), nrow = 2, byrow = TRUE)
rownames(ctab4) <- c("not_targets", "targets")
colnames(ctab4) <- c("mat_early", "other")
ctab4
fisher.test(ctab4)

ctab5 <- matrix(c(ctab[1,4], ctab[1,1]+ctab[1,2]+ctab[1,3], ctab[2,4], ctab[2,1]+ctab[2,2]+ctab[2,3]), nrow = 2, byrow = TRUE)
rownames(ctab5) <- c("not_targets", "targets")
colnames(ctab5) <- c("mat_late", "other")
ctab5
fisher.test(ctab5)

# add in the "none" category of genes, those not diff. expressed in dataset
ctab6 <- matrix(c(ctab[1,1], ctab[1,2]+ctab[1,3]+ctab[1,4]+length(genes_none)-length(targets_none), 
                  ctab[2,1], ctab[2,2]+ctab[2,3]+ctab[2,4]+length(targets_none)), nrow = 2, byrow = TRUE)
rownames(ctab6) <- c("not_targets", "targets")
colnames(ctab6) <- c("gz>cp", "other")
ctab6
fisher.test(ctab6)

ctab3 <- matrix(c(ctab[1,2], ctab[1,1]+ctab[1,3]+ctab[1,4], ctab[2,2], ctab[2,1]+ctab[2,3]+ctab[2,4]), nrow = 2, byrow = TRUE)
rownames(ctab3) <- c("not_targets", "targets")
colnames(ctab3) <- c("cp>gz", "other")
ctab3
fisher.test(ctab3)

ctab4 <- matrix(c(ctab[1,3], ctab[1,1]+ctab[1,2]+ctab[1,4], ctab[2,3], ctab[2,1]+ctab[2,2]+ctab[2,4]), nrow = 2, byrow = TRUE)
rownames(ctab4) <- c("not_targets", "targets")
colnames(ctab4) <- c("mat_early", "other")
ctab4
fisher.test(ctab4)

ctab5 <- matrix(c(ctab[1,4], ctab[1,1]+ctab[1,2]+ctab[1,3], ctab[2,4], ctab[2,1]+ctab[2,2]+ctab[2,3]), nrow = 2, byrow = TRUE)
rownames(ctab5) <- c("not_targets", "targets")
colnames(ctab5) <- c("mat_late", "other")
ctab5
fisher.test(ctab5)

# version 2:
# Fisher's Exact Test ######

# all genes in analysis are in gz>cp or other
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
genes_other <- gene_df$ensg[which(gene_df$category != "neurogenesis_gz")]

# the targets that overlap gz>cp or other
targets_gz_greater_than_cp <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_gz_greater_than_cp)]
targets_other <- predictions$ENSEMBL[which(predictions$ENSEMBL %in% genes_other)]

# cont. table: genes are either targets or not targets, and genes are either gz>cp or other
tab <- matrix(c(length(genes_other)-length(targets_other), length(genes_gz_greater_than_cp)-length(targets_gz_greater_than_cp),
                length(targets_other), length(targets_gz_greater_than_cp)), byrow = TRUE, nrow = 2)

# labels rows and columns for this table
colnames(tab) <- c("other", "gz>cp")
rownames(tab) <- c("not_targets", "targets")

tab

ft <- fisher.test(tab)
ft



