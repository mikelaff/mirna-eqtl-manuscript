# look at miRNA targets and their overlap with differentially expressed genes

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(AnnotationHub)
library(miRNAtap)
library(EnsDb.Hsapiens.v86)
library(reshape2)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_overlap/pdfs/")
# write plots
write_plots <- FALSE
# data frame with fisher test results
fisher_test_overlap_rds <- here("doc/target_overlap/rdata/mirna_fisher_test_overlap_results_v2.rds")

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

# Filter Genes #####################################################################################
edb <- EnsDb.Hsapiens.v86
prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

gene_df %<>%
  dplyr::filter(ensg %in% prot_genes_ensg)

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
# number of targets not overlapping gz>cp genes
mirna_df$num_targets_not_gz_over_cp <- NA
# number of targets not overlapping cp>gz genes
mirna_df$num_targets_not_cp_over_gz <- NA
# number of targets not overlapping maturation early genes
mirna_df$num_targets_not_mat_early <- NA
# number of targets not overlapping maturation late genes
mirna_df$num_targets_not_mat_late <- NA

# results of fisher's exact for gz>cp enrichment
mirna_df$FT_gz_cp_odds <- NA
mirna_df$FT_gz_cp_pval <- NA
mirna_df$FT_gz_cp_ci_low <- NA
mirna_df$FT_gz_cp_ci_high <- NA
# results of fisher's exact for cp>gz enrichment
mirna_df$FT_cp_gz_odds <- NA
mirna_df$FT_cp_gz_pval <- NA
mirna_df$FT_cp_gz_ci_low <- NA
mirna_df$FT_cp_gz_ci_high <- NA
# results of fisher's exact for mat_early
mirna_df$FT_mat_early_odds <- NA
mirna_df$FT_mat_early_pval <- NA
mirna_df$FT_mat_early_ci_low <- NA
mirna_df$FT_mat_early_ci_high <- NA
# results of fisher's exact for mat_late
mirna_df$FT_mat_late_odds <- NA
mirna_df$FT_mat_late_pval <- NA
mirna_df$FT_mat_late_ci_low <- NA
mirna_df$FT_mat_late_ci_high <- NA

# all genes in each group tested
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp")]
genes_maturation_early <- gene_df$ensg[which(gene_df$category == "maturation_early")]
genes_maturation_late <- gene_df$ensg[which(gene_df$category == "maturation_late")]

genes_not_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category != "neurogenesis_gz")]
genes_not_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category != "neurogenesis_cp")]
genes_not_maturation_early <- gene_df$ensg[which(gene_df$category != "maturation_early")]
genes_not_maturation_late <- gene_df$ensg[which(gene_df$category != "maturation_late")]

# totals
num_genes_gz_greater_than_cp <- length(genes_gz_greater_than_cp)
num_genes_cp_greater_than_gz <- length(genes_cp_greater_than_gz)
num_genes_maturation_early <- length(genes_maturation_early)
num_genes_maturation_late <- length(genes_maturation_late)

num_genes_not_gz_greater_than_cp <- length(genes_not_gz_greater_than_cp)
num_genes_not_cp_greater_than_gz <- length(genes_not_cp_greater_than_gz)
num_genes_not_maturation_early <- length(genes_not_maturation_early)
num_genes_not_maturation_late <- length(genes_not_maturation_late)

# looking at only sig. expressed mirnas (all categories)
# loop and find targets and overlaps
indexes <- which(mirna_df$category != "none")
for (i in indexes) {
  ft_gzcp <- NULL
  ft_cpgz <- NULL
  ft_early <- NULL
  ft_late <- NULL
  predictions <- NULL
  rowdata <- NULL
  # short name used for database lookup
  mirna_df$mir[i] <- strsplit(strsplit(mirna_df$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(paste(i, mirna_df$mir[i], sep=": "))
  # get predicted targets
  predictions <- data.frame(getPredictedTargets(mirna_df$mir[i], species = "hsa", method = "geom", min_src = min_sources))
  if (length(predictions) == 0) {
    print(paste("No Targets:", mirna_df$mirna[i]))
    next
  }
  predictions$ENTREZID <- rownames(predictions)
  # get ensembl ids
  rowdata <- AnnotationDbi::select(orgdb,
                                   keys = predictions$ENTREZID,
                                   keytype = "ENTREZID",
                                   columns = c("ENSEMBL"))
  # join with predictions
  predictions <- left_join(predictions, rowdata, by = "ENTREZID")

  # populate table
  # number of targets
  mirna_df$num_targets[i] <- length(predictions$ENTREZID)
  # number of targets overlapping gz>cp genes
  mirna_df$num_targets_gz_over_cp[i] <- sum(predictions$ENSEMBL %in% genes_gz_greater_than_cp)
  # number of targets overlapping cp>gz genes
  mirna_df$num_targets_cp_over_gz[i] <- sum(predictions$ENSEMBL %in% genes_cp_greater_than_gz)
  # number of targets overlapping maturation early genes
  mirna_df$num_targets_mat_early[i] <- sum(predictions$ENSEMBL %in% genes_maturation_early)
  # number of targets overlapping maturation late genes
  mirna_df$num_targets_mat_late[i] <- sum(predictions$ENSEMBL %in% genes_maturation_late)
  # number of targets not overlapping gz>cp genes
  mirna_df$num_targets_not_gz_over_cp[i] <- sum(predictions$ENSEMBL %in% genes_not_gz_greater_than_cp)
  # number of targets not overlapping cp>gz genes
  mirna_df$num_targets_not_cp_over_gz[i] <- sum(predictions$ENSEMBL %in% genes_not_cp_greater_than_gz)
  # number of targets not overlapping maturation early genes
  mirna_df$num_targets_not_mat_early[i] <- sum(predictions$ENSEMBL %in% genes_not_maturation_early)
  # number of targets not overlapping maturation late genes
  mirna_df$num_targets_not_mat_late[i] <- sum(predictions$ENSEMBL %in% genes_not_maturation_late)


  # # cont. table: genes are either targets or not targets, and genes are either gz>cp or other
  # tab <- matrix(c(num_genes_not_gz_greater_than_cp-mirna_df$num_targets_not_gz_over_cp[i],
  #                 num_genes_gz_greater_than_cp-mirna_df$num_targets_gz_over_cp[i],
  #                 mirna_df$num_targets_not_gz_over_cp[i],
  #                 mirna_df$num_targets_gz_over_cp[i]), byrow = TRUE, nrow = 2)
  # 
  # # labels rows and columns for this table
  # colnames(tab) <- c("other", "gz>cp")
  # rownames(tab) <- c("not_targets", "targets")
  # tab

  # gz>cp fisher's exact test
  ft_gzcp <- fisher.test(matrix(c(num_genes_not_gz_greater_than_cp-mirna_df$num_targets_not_gz_over_cp[i],
                                  num_genes_gz_greater_than_cp-mirna_df$num_targets_gz_over_cp[i],
                                  mirna_df$num_targets_not_gz_over_cp[i],
                                  mirna_df$num_targets_gz_over_cp[i]), byrow = TRUE, nrow = 2))

  mirna_df$FT_gz_cp_odds[i] <- ft_gzcp$estimate
  mirna_df$FT_gz_cp_pval[i] <- ft_gzcp$p.value
  mirna_df$FT_gz_cp_ci_low[i] <- ft_gzcp$conf.int[1]
  mirna_df$FT_gz_cp_ci_high[i] <- ft_gzcp$conf.int[2]

  # cp>gz fisher's exact test
  ft_cpgz <- fisher.test(matrix(c(num_genes_not_cp_greater_than_gz-mirna_df$num_targets_not_cp_over_gz[i],
                                  num_genes_cp_greater_than_gz-mirna_df$num_targets_cp_over_gz[i],
                                  mirna_df$num_targets_not_cp_over_gz[i],
                                  mirna_df$num_targets_cp_over_gz[i]), byrow = TRUE, nrow = 2))

  mirna_df$FT_cp_gz_odds[i] <- ft_cpgz$estimate
  mirna_df$FT_cp_gz_pval[i] <- ft_cpgz$p.value
  mirna_df$FT_cp_gz_ci_low[i] <- ft_cpgz$conf.int[1]
  mirna_df$FT_cp_gz_ci_high[i] <- ft_cpgz$conf.int[2]

  # maturation early fisher's exact test
  ft_early <- fisher.test(matrix(c(num_genes_not_maturation_early-mirna_df$num_targets_not_mat_early[i],
                                   num_genes_maturation_early-mirna_df$num_targets_mat_early[i],
                                   mirna_df$num_targets_not_mat_early[i],
                                   mirna_df$num_targets_mat_early[i]), byrow = TRUE, nrow = 2))

  mirna_df$FT_mat_early_odds[i] <- ft_early$estimate
  mirna_df$FT_mat_early_pval[i] <- ft_early$p.value
  mirna_df$FT_mat_early_ci_low[i] <- ft_early$conf.int[1]
  mirna_df$FT_mat_early_ci_high[i] <- ft_early$conf.int[2]

  # maturation late fisher's exact test
  ft_late <- fisher.test(matrix(c(num_genes_not_maturation_late-mirna_df$num_targets_not_mat_late[i],
                                   num_genes_maturation_late-mirna_df$num_targets_mat_late[i],
                                   mirna_df$num_targets_not_mat_late[i],
                                   mirna_df$num_targets_mat_late[i]), byrow = TRUE, nrow = 2))

  mirna_df$FT_mat_late_odds[i] <- ft_late$estimate
  mirna_df$FT_mat_late_pval[i] <- ft_late$p.value
  mirna_df$FT_mat_late_ci_low[i] <- ft_late$conf.int[1]
  mirna_df$FT_mat_late_ci_high[i] <- ft_late$conf.int[2]
}

# save data
#saveRDS(mirna_df, fisher_test_overlap_rds)
# load data
mirna_df <- readRDS(fisher_test_overlap_rds)

# basic l2fc plot
mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

pdf("~/Desktop/mirna_target_enrichment.pdf")

# num_gzcp_test <- length(dplyr::filter(mirna_df, category == "neurogenesis_gz", !is.na(FT_gz_cp_odds))$mirna)
# num_cpgz_test <- length(dplyr::filter(mirna_df, category == "neurogenesis_cp", !is.na(FT_cp_gz_odds))$mirna)
# num_early_test <- length(dplyr::filter(mirna_df, category == "maturation_early", !is.na(FT_mat_early_odds))$mirna)
# num_late_test <- length(dplyr::filter(mirna_df, category == "maturation_late", !is.na(FT_mat_late_odds))$mirna)

# mirna upregulated gz: GZ>CP
mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated cp: CP>GZ
mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated early: maturation_early
mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated late: maturation_late
mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "BH")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

dev.off()

###########################
mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_gz_cp_odds,
                pval = FT_gz_cp_pval,
                ci_low = FT_gz_cp_ci_low,
                ci_high = FT_gz_cp_ci_high) %>%
  dplyr::mutate(target_cat = "gz>cp") -> tmp1

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_cp_gz_odds,
                pval = FT_cp_gz_pval,
                ci_low = FT_cp_gz_ci_low,
                ci_high = FT_cp_gz_ci_high) %>%
  dplyr::mutate(target_cat = "cp>gz") -> tmp2

tmp3 <- bind_rows(tmp1, tmp2)

tmp3 %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), color = target_cat)) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

#################
mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_gz_cp_odds,
                pval = FT_gz_cp_pval,
                ci_low = FT_gz_cp_ci_low,
                ci_high = FT_gz_cp_ci_high) %>%
  dplyr::mutate(target_cat = "gz>cp") -> tmp1

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_cp_gz_odds,
                pval = FT_cp_gz_pval,
                ci_low = FT_cp_gz_ci_low,
                ci_high = FT_cp_gz_ci_high) %>%
  dplyr::mutate(target_cat = "cp>gz") -> tmp2

tmp3 <- bind_rows(tmp1, tmp2)

tmp3 %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), color = target_cat)) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# SfN Plotting ##########################

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="blue") +
  geom_text(aes(label=mir), hjust=-0.5)+
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","1.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="red") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","2.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="green") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","3.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="orange") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","4.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="purple") +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","5.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="navy") +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","6.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="pink") +
  plotTheme +
  geom_text(aes(label=mir), hjust=-1, size=2) +
  labs(title="miRNA Up Late: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","7.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "BH")))) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="yellow") +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","8.pdf"), width = 5, height = 5, useDingbats = FALSE)

# basic l2fc plot
circColor <- "black"
labelColor <- "black"

mirna_df %>%
  dplyr::filter(category != "none") %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=3) +
  presTheme +
  labs(x="Log2 Fold Change Gest. Week",
       y="Log2 Fold Change GZ/CP") +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  theme(legend.position = "none") +
  geom_label(data=subset(mirna_df, mir == "miR-122-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-122-5p"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-122-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mir == "miR-106b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-106b-5p"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-106b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mir == "miR-19a-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-19a-3p"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-19a-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

ggsave(file.path("~/Desktop/","mirna_category_mir122.pdf"), width = 11, height = 9, useDingbats = FALSE)


ggsave(file.path("~/Desktop/","mirna_category.pdf"), width = 11, height = 9, useDingbats = FALSE)

gene_df %>%
  dplyr::filter(category != "none") %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=1) +
  presTheme +
  labs(x="Log2 Fold Change Gest. Week",
       y="Log2 Fold Change GZ/CP") +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  theme(legend.position = "none") +
  geom_label(data=subset(gene_df, ensg == "ENSG00000128573"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="FOXP2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000128573"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000111249"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="CUX2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000111249"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000184486"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="BRN2"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000184486"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000127152"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="CTIP2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000127152"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)
  
ggsave(file.path("~/Desktop/","gene_category.pdf"), width = 9, height = 9, useDingbats = FALSE)

mir <- "miR-106b-5p"

# get predicted targets
predictions <- data.frame(getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2))
predictions$ENTREZID <- rownames(predictions)
# get ensembl ids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = predictions$ENTREZID,
                                 keytype = "ENTREZID",
                                 columns = c("ENSEMBL", "SYMBOL"))
# join with predictions
predictions <- left_join(predictions, rowdata, by = "ENTREZID")

# join with genes
predictions <- left_join(predictions, gene_df, by = c("ENSEMBL" = "ensg"))

predictions %<>% arrange(-log2FoldChange.gw)
write_lines(predictions$ENSEMBL[1:200], "~/Desktop/mir_106b_targets_up_late_200.txt")


mir <- "miR-19a-3p"

# get predicted targets
predictions <- data.frame(getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2))
predictions$ENTREZID <- rownames(predictions)
# get ensembl ids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = predictions$ENTREZID,
                                 keytype = "ENTREZID",
                                 columns = c("ENSEMBL", "SYMBOL"))
# join with predictions
predictions <- left_join(predictions, rowdata, by = "ENTREZID")

# join with genes
predictions <- left_join(predictions, gene_df, by = c("ENSEMBL" = "ensg"))

predictions %<>% arrange(-log2FoldChange.gw)
write_lines(predictions$ENSEMBL[1:200], "~/Desktop/mir_19a_targets_up_late_200.txt")


