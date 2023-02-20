# look at miRNA targets and their overlap with differentially expressed genes

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(AnnotationHub)
library(miRNAtap)
library(EnsDb.Hsapiens.v86)
library(readxl)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_overlap/pdfs/")
# write plots
write_plots <- FALSE
# data frame with fisher test results
fisher_test_overlap_rds <- here("doc/target_overlap/rdata/mirna_fisher_test_overlap_results_valTargets_v2.rds")

# INPUT FILES #####################################################################################
# combined deseq2 shrunken results of gene expression in gzcp or gw as a dataframe with labels for 
# significance and category (neuro vs maturation)
gene_expression_combined_df_rds <- here("doc/diff_expression/rdata/combined_gene_expression_deseq2_shrunken_results_df.rds")
# combined deseq2 shrunken results of mirna expression in gzcp or gw as a dataframe with labels for 
# significance and category (neuro vs maturation)
mirna_expression_combined_df_rds <- here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_df.rds")
# mirtarbase v7 hsa targets
mirna_targets_xlsx <- here("data/mirtarbase_v7/hsa_MTI.xlsx")

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

# Load Targets ####################################################################################
targets <- read_xlsx(mirna_targets_xlsx)
# trim dataframe
targets %<>% dplyr::select(mirna = miRNA, target_gene = `Target Gene`, ENTREZID = `Target Gene (Entrez Gene ID)`)
# only keep unique mirna-target interactions (remove rows from multiple references)
targets %<>% dplyr::filter(!duplicated(paste(targets$mirna, targets$ENTREZID, sep="_")))
targets$ENTREZID <- as.character(targets$ENTREZID)

entrezIDs <- as.character(targets$ENTREZID)
entrezIDs <- entrezIDs[!duplicated(entrezIDs)]

# get ensembl ids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = entrezIDs,
                                 keytype = "ENTREZID",
                                 columns = c("ENSEMBL"))

targets <- left_join(targets, rowdata, by = "ENTREZID")

# Find Target Overlaps ######################################

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
  # short name used for database lookup
  mirna_df$mir[i] <- strsplit(mirna_df$mirna[i], '/')[[1]][1]
  print(paste(i, mirna_df$mir[i], sep=": "))
  # get valid targets
  predictions <- dplyr::filter(targets, mirna == mirna_df$mir[i])$ENSEMBL
  predictions <- predictions[!duplicated(predictions)]
  if (length(predictions) == 0) {
    print(paste("No Targets:", mirna_df$mirna[i]))
    next
  }

  # populate table
  # number of targets
  mirna_df$num_targets[i] <- length(predictions)
  # number of targets overlapping gz>cp genes
  mirna_df$num_targets_gz_over_cp[i] <- sum(predictions %in% genes_gz_greater_than_cp)
  # number of targets overlapping cp>gz genes
  mirna_df$num_targets_cp_over_gz[i] <- sum(predictions %in% genes_cp_greater_than_gz)
  # number of targets overlapping maturation early genes
  mirna_df$num_targets_mat_early[i] <- sum(predictions %in% genes_maturation_early)
  # number of targets overlapping maturation late genes
  mirna_df$num_targets_mat_late[i] <- sum(predictions %in% genes_maturation_late)
  # number of targets not overlapping gz>cp genes
  mirna_df$num_targets_not_gz_over_cp[i] <- sum(predictions %in% genes_not_gz_greater_than_cp)
  # number of targets not overlapping cp>gz genes
  mirna_df$num_targets_not_cp_over_gz[i] <- sum(predictions %in% genes_not_cp_greater_than_gz)
  # number of targets not overlapping maturation early genes
  mirna_df$num_targets_not_mat_early[i] <- sum(predictions %in% genes_not_maturation_early)
  # number of targets not overlapping maturation late genes
  mirna_df$num_targets_not_mat_late[i] <- sum(predictions %in% genes_not_maturation_late)


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
saveRDS(mirna_df, fisher_test_overlap_rds)
# load data
mirna_df <- readRDS(fisher_test_overlap_rds)

# basic l2fc plot
mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

pdf("~/Desktop/mirna_target_enrichment_valTargets.pdf")

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
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

dev.off()

############
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
