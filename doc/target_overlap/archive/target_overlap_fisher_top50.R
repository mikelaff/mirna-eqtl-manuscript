# look at miRNA targets and their overlap with differentially expressed genes

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
# data frame with fisher test results
fisher_test_overlap_rds <- here("doc/target_overlap/rdata/mirna_fisher_test_overlap_top50_results.rds")

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

# Find Target Overlaps ######################################
# minimum number of sources required for a target to be considered
min_sources <- 2

# short name used for DB lookup
mirna_df$mir <- NA
# number of targets
mirna_df$num_targets <- NA
# number of targets overlapping gz>cp genes
mirna_df$num_targets_gz_over_cp <- NA
# number of targets overlapping cp>gz genes
mirna_df$num_targets_cp_over_gz <- NA
# number of targets overlapping maturation early genes
mirna_df$num_targets_mat_early <- NA
# number of targets overlapping maturation late genes
mirna_df$num_targets_mat_late <- NA
# number of targets overlapping none category (not diff. expressed, but present in dataset)
mirna_df$num_targets_notDE <- NA

# results of fisher's exact for gz>cp vs cp>gz
mirna_df$FT_gz_cp_odds <- NA
mirna_df$FT_gz_cp_pval <- NA
# results of fisher's exact for mat_early vs mat_late
mirna_df$FT_early_late_odds <- NA
mirna_df$FT_early_late_pval <- NA

# all genes are in one of 5 groups: gz>cp, cp>gz, mat_early, mat_late, or none (notDE)
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp")]
genes_maturation_early <- gene_df$ensg[which(gene_df$category == "maturation_early")]
genes_maturation_late <- gene_df$ensg[which(gene_df$category == "maturation_late")]
genes_notDE <- gene_df$ensg[which(gene_df$category == "none")]

# totals
num_genes_gz_greater_than_cp <- length(genes_gz_greater_than_cp)
num_genes_cp_greater_than_gz <- length(genes_cp_greater_than_gz)
num_genes_maturation_early <- length(genes_maturation_early)
num_genes_maturation_late <- length(genes_maturation_late)
num_genes_notDE <- length(genes_notDE)

# looking at only sig. expressed mirnas (all categories)
# loop and find targets and overlaps
# indexes <- which(mirna_df$category != "none")
# for (i in indexes) {
#   ft1 <- NULL
#   ft2 <- NULL
#   predictions <- NULL
#   rowdata <- NULL
#   # short name used for database lookup
#   mirna_df$mir[i] <- strsplit(strsplit(mirna_df$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
#   print(paste(i, mirna_df$mir[i], sep=": "))
#   # get predicted targets
#   predictions <- data.frame(getPredictedTargets(mirna_df$mir[i], species = "hsa", method = "geom", min_src = min_sources))
#   if (length(predictions) == 0) {
#     print(paste("No Targets:", mirna_df$mirna[i]))
#     next
#   }
#   predictions$ENTREZID <- rownames(predictions)
#   if (length(predictions$ENTREZID) < 50) {
#     print("Not 50 targets")
#     next
#   }
#   # top 50 predictions
#   predictions <- predictions[1:50,]
#   # get ensembl ids
#   rowdata <- AnnotationDbi::select(orgdb,
#                                    keys = predictions$ENTREZID,
#                                    keytype = "ENTREZID",
#                                    columns = c("ENSEMBL"))
#   # join with predictions
#   predictions <- left_join(predictions, rowdata, by = "ENTREZID")
# 
#   # populate table
#   # number of targets
#   mirna_df$num_targets[i] <- length(predictions$ENTREZID)
#   # number of targets overlapping gz>cp genes
#   mirna_df$num_targets_gz_over_cp[i] <- sum(predictions$ENSEMBL %in% genes_gz_greater_than_cp)
#   # number of targets overlapping cp>gz genes
#   mirna_df$num_targets_cp_over_gz[i] <- sum(predictions$ENSEMBL %in% genes_cp_greater_than_gz)
#   # number of targets overlapping maturation early genes
#   mirna_df$num_targets_mat_early[i] <- sum(predictions$ENSEMBL %in% genes_maturation_early)
#   # number of targets overlapping maturation late genes
#   mirna_df$num_targets_mat_late[i] <- sum(predictions$ENSEMBL %in% genes_maturation_late)
#   # number of targets overlapping none category (not diff. expressed, but present in dataset)
#   mirna_df$num_targets_notDE[i] <- sum(predictions$ENSEMBL %in% genes_notDE)
# 
#   # fisher test gz>cp vs cp>gz
#   ft1 <- fisher.test(matrix(c(num_genes_gz_greater_than_cp-mirna_df$num_targets_gz_over_cp[i],
#                              num_genes_cp_greater_than_gz-mirna_df$num_targets_cp_over_gz[i],
#                              mirna_df$num_targets_gz_over_cp[i],
#                              mirna_df$num_targets_cp_over_gz[i]), nrow = 2, byrow = TRUE))
#   mirna_df$FT_gz_cp_odds[i] <- ft1$estimate
#   mirna_df$FT_gz_cp_pval[i] <- ft1$p.value
# 
#   # fisher test mat_early vs mat_late
#   ft2 <- fisher.test(matrix(c(num_genes_maturation_early-mirna_df$num_targets_mat_early[i],
#                              num_genes_maturation_late-mirna_df$num_targets_mat_late[i],
#                              mirna_df$num_targets_mat_early[i],
#                              mirna_df$num_targets_mat_late[i]), nrow = 2, byrow = TRUE))
#   mirna_df$FT_early_late_odds[i] <- ft2$estimate
#   mirna_df$FT_early_late_pval[i] <- ft2$p.value
# 
# 
# }

#saveRDS(mirna_df, fisher_test_overlap_rds)
mirna_df <- readRDS(fisher_test_overlap_rds)

# filter for expressed mirnas with sufficient targets
mirna_df %>%
  filter(baseMean.gw > 100 | baseMean.gzcp > 1000) %>%
  filter(num_targets > 50) -> tmp

# basic l2fc plot
mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

tmp %>%
  filter(category == "neurogenesis_gz") %>%
  ggplot(aes(x=FT_gz_cp_odds, y=-log10(FT_gz_cp_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp %>%
  filter(category == "neurogenesis_cp") %>%
  ggplot(aes(x=FT_gz_cp_odds, y=-log10(FT_gz_cp_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp %>%
  filter(category == "maturation_early") %>%
  ggplot(aes(x=FT_early_late_odds, y=-log10(FT_early_late_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp %>%
  filter(category == "maturation_late") %>%
  ggplot(aes(x=FT_early_late_odds, y=-log10(FT_early_late_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

# filter for significant mirna
tmp %>%
  filter(FT_gz_cp_pval < 0.05 | FT_early_late_pval < 0.05) -> tmp2


mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  geom_point(data=tmp2, aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), shape=21, size=3, color="black", stroke=1)

tmp2 %>%
  filter(category == "neurogenesis_gz") %>%
  ggplot(aes(x=FT_gz_cp_odds, y=-log10(FT_gz_cp_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp2 %>%
  filter(category == "neurogenesis_cp") %>%
  ggplot(aes(x=FT_gz_cp_odds, y=-log10(FT_gz_cp_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp2 %>%
  filter(category == "maturation_early") %>%
  ggplot(aes(x=FT_early_late_odds, y=-log10(FT_early_late_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))

tmp2 %>%
  filter(category == "maturation_late") %>%
  ggplot(aes(x=FT_early_late_odds, y=-log10(FT_early_late_pval))) +
  geom_point() +
  plotTheme +
  geom_hline(yintercept = -log10(0.05))
