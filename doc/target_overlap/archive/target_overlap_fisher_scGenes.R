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
library(gridExtra)

source(here("src/utils/lafferty_utils.R"))

localTheme <- theme(title = element_text(size = 9),
                    axis.title = element_text(size = 9),
                    axis.text = element_text(size = 8),
                    plot.caption = element_text(size = 9))
# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_overlap/pdfs/")
# write plots
write_plots <- FALSE

# INPUT FILES #####################################################################################
# combined deseq2 shrunken results of mirna expression in gzcp or gw as a dataframe with labels for 
# significance and category (neuro vs maturation)
mirna_expression_combined_df_rds <- here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_df.rds")
# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes_by_cluster_csv <- here("data/ucla_single_cell/TableS6 Cluster enriched genes.csv")

# Import Data #####################################################################################
mirna_df_perm <- readRDS(mirna_expression_combined_df_rds)

genes_by_cluster_df <- read_csv(genes_by_cluster_csv)

# Annotation Hub ##################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Clusters ###########################
clusters <- unique(genes_by_cluster_df$Cluster)


# Filter Genes #####################################################################################
edb <- EnsDb.Hsapiens.v86
prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

genes_by_cluster_df %<>%
  dplyr::filter(Ensembl %in% prot_genes_ensg)




# loop over each cluster
#cluster <- "ExN"

for (cluster in clusters) {
  
  mirna_df <- mirna_df_perm
  
  # Find Target Overlaps ######################################
  # minimum number of sources required for a target to be considered
  min_sources <- 2
  
  # short name used for DB lookup
  mirna_df$mir <- NA
  # number of targets
  mirna_df$num_targets <- NA
  # number of targets overlapping cell type specific genes
  mirna_df$num_targets_celltype_specific <- NA
  # number of targets not overlapping cell type specific genes
  mirna_df$num_targets_not_celltype_specific <- NA
  
  # results of fisher's exact for cell type specific gene enrichment
  mirna_df$FT_odds <- NA
  mirna_df$FT_pval <- NA
  mirna_df$FT_ci_low <- NA
  mirna_df$FT_ci_high <- NA
  
  # all genes in each group tested
  genes_celltype_specific <- genes_by_cluster_df$Ensembl[which(genes_by_cluster_df$Cluster == cluster)]
  
  genes_not_celltype_specific <- genes_by_cluster_df$Ensembl[which(genes_by_cluster_df$Cluster != cluster)]
  
  # totals
  num_genes_celltype_specific <- length(genes_celltype_specific)
  
  num_genes_not_celltype_specific <- length(genes_not_celltype_specific)
  
  # looking at only sig. expressed mirnas (all categories)
  # loop and find targets and overlaps
  indexes <- which(mirna_df$category != "none")
  for (i in indexes) {
    ft <- NULL
    predictions <- NULL
    rowdata <- NULL
    # short name used for database lookup
    mirna_df$mir[i] <- strsplit(strsplit(mirna_df$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
    print(paste(i, mirna_df$mir[i], cluster, sep=": "))
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
    # number of targets overlapping cell type specific genes
    mirna_df$num_targets_celltype_specific[i] <- sum(predictions$ENSEMBL %in% genes_celltype_specific)
    # number of targets not overlapping cell type specific genes
    mirna_df$num_targets_not_celltype_specific[i] <- sum(predictions$ENSEMBL %in% genes_not_celltype_specific)
  
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
  
    # fisher's exact test
    ft <- fisher.test(matrix(c(num_genes_not_celltype_specific-mirna_df$num_targets_not_celltype_specific[i],
                                    num_genes_celltype_specific-mirna_df$num_targets_celltype_specific[i],
                                    mirna_df$num_targets_not_celltype_specific[i],
                                    mirna_df$num_targets_celltype_specific[i]), byrow = TRUE, nrow = 2))
  
    mirna_df$FT_odds[i] <- ft$estimate
    mirna_df$FT_pval[i] <- ft$p.value
    mirna_df$FT_ci_low[i] <- ft$conf.int[1]
    mirna_df$FT_ci_high[i] <- ft$conf.int[2]
  }
  
  mirna_df %>%
    dplyr::filter(category == "neurogenesis_gz", !is.na(FT_odds)) %>%
    dplyr::select(odds = FT_odds, pval = FT_pval, ci_low = FT_ci_low, ci_high = FT_ci_high, mir=mir) %>%
    ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
    geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
    geom_point(color="blue") +
    plotTheme + localTheme +
    labs(title=paste("miRNA GZ>CP: Target Enrichment in", cluster),
         caption=paste("odds > 0: targets enriched among", cluster, "genes\nodds < 0: targets depleted among", cluster, "genes"),
         y="-log10(adjusted pval)") +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_x_continuous(limits = c(-1,1)) -> p1
  
  mirna_df %>%
    dplyr::filter(category == "neurogenesis_cp", !is.na(FT_odds)) %>%
    dplyr::select(odds = FT_odds, pval = FT_pval, ci_low = FT_ci_low, ci_high = FT_ci_high, mir=mir) %>%
    ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
    geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
    geom_point(color="red") +
    plotTheme + localTheme +
    labs(title=paste("miRNA CP>GZ: Target Enrichment in", cluster),
         caption=paste("odds > 0: targets enriched among", cluster, "genes\nodds < 0: targets depleted among", cluster, "genes"),
         y="-log10(adjusted pval)") +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_x_continuous(limits = c(-1,1)) -> p2
  
  mirna_df %>%
    dplyr::filter(category == "maturation_early", !is.na(FT_odds)) %>%
    dplyr::select(odds = FT_odds, pval = FT_pval, ci_low = FT_ci_low, ci_high = FT_ci_high, mir=mir) %>%
    ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
    geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
    geom_point(color="green") +
    plotTheme + localTheme +
    labs(title=paste("miRNA Up Early: Target Enrichment in", cluster),
         caption=paste("odds > 0: targets enriched among", cluster, "genes\nodds < 0: targets depleted among", cluster, "genes"),
         y="-log10(adjusted pval)") +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_x_continuous(limits = c(-1,1)) -> p3
  
  mirna_df %>%
    dplyr::filter(category == "maturation_late", !is.na(FT_odds)) %>%
    dplyr::select(odds = FT_odds, pval = FT_pval, ci_low = FT_ci_low, ci_high = FT_ci_high, mir=mir) %>%
    ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
    geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
    geom_point(color="orange") +
    plotTheme + localTheme +
    labs(title=paste("miRNA Up Late: Target Enrichment in", cluster),
         caption=paste("odds > 0: targets enriched among", cluster, "genes\nodds < 0: targets depleted among", cluster, "genes"),
         y="-log10(adjusted pval)") +
    geom_hline(yintercept = -log10(0.05), linetype = 2) +
    geom_vline(xintercept = 0, linetype = 2) +
    scale_x_continuous(limits = c(-1,1)) -> p4
  
  g <- do.call(grid.arrange,list(p1,p2,p3,p4))
  ggsave(filename = paste(pdf_dir, cluster, ".pdf", sep=""), plot = g)

}
