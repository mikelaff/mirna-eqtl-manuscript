# create files for GSEA

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(AnnotationHub)
library(miRNAtap)

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/target_overlap/pdfs/")
# write plots
write_plots <- FALSE
# GSEA Gene Set Database File (.gmt)
#gmt_file <- here("results/GSEA/gene_sets.gmt")
# rnk files directory
rnk_dir <- here("results/GSEA/rnk_files/")

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

#mirna_df$mirna_clean <- sapply(strsplit(mirna_df$mirna, '/'), `[`, 1)

# Annotation Hub ##################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)


# Gene Set Database
# http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#Gene_Set_Database_Formats
# GMT (.gmt)
# tab delimited file format, each row represents a gene set
# first column are gene set names
# second column is description (can be na)
# following columns are genes

# # all genes are in one of 4 groups: gz>cp, cp>gz, mat_early, mat_late4
# genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz")]
# genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp")]
# genes_maturation_early <- gene_df$ensg[which(gene_df$category == "maturation_early")]
# genes_maturation_late <- gene_df$ensg[which(gene_df$category == "maturation_late")]
# 
# # write .grp files
# write_lines(genes_gz_greater_than_cp, here("results/GSEA/gz_greater_than_cp.grp"))
# write_lines(genes_cp_greater_than_gz, here("results/GSEA/cp_greater_than_gz.grp"))
# write_lines(genes_maturation_early, here("results/GSEA/maturation_early.grp"))
# write_lines(genes_maturation_late, here("results/GSEA/maturation_late.grp"))

# all genes are in one of 4 groups: gz>cp, cp>gz, mat_early, mat_late4
genes_gz_greater_than_cp <- gene_df$ensg[which(gene_df$category == "neurogenesis_gz" & gene_df$baseMean.gzcp > 1000)]
genes_cp_greater_than_gz <- gene_df$ensg[which(gene_df$category == "neurogenesis_cp" & gene_df$baseMean.gzcp > 1000)]
genes_maturation_early <- gene_df$ensg[which(gene_df$category == "maturation_early" & gene_df$baseMean.gw > 1000)]
genes_maturation_late <- gene_df$ensg[which(gene_df$category == "maturation_late" & gene_df$baseMean.gw > 1000)]

# write .grp files
write_lines(genes_gz_greater_than_cp, here("results/GSEA/gz_greater_than_cp_high_expressers.grp"))
write_lines(genes_cp_greater_than_gz, here("results/GSEA/cp_greater_than_gz_high_expressers.grp"))
write_lines(genes_maturation_early, here("results/GSEA/maturation_early_high_expressers.grp"))
write_lines(genes_maturation_late, here("results/GSEA/maturation_late_high_expressers.grp"))

min_sources <- 2

# write rnk files for each mirna
indexes <- which(mirna_df$category != "none")
for (i in indexes) {
  predictions <- NULL
  rowdata <- NULL
  mir <- NULL
  # short name used for database lookup
  mir <- strsplit(strsplit(mirna_df$mirna[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(paste(i, mir, sep=": "))
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

  # write rnk file
  write_tsv(data.frame(predictions$ENSEMBL, length(predictions$rank_final)-predictions$rank_final),
            file.path(rnk_dir, mirna_df$category[i], paste(mir, ".rnk", sep="")), col_names = FALSE)
}



# # for GSEA pre ranked
# # list of targets ranked
# write_tsv(data.frame(predictions$ENSEMBL, length(predictions$rank_final)-predictions$rank_final),
#           here("results/GSEA/miR-92b-3b_ranked_targets.rnk"),col_names = FALSE)






















# scratch
gene_df %>%
  #filter(category == "neurogenesis_gz") %>%
  arrange(-log2FoldChange.gzcp) -> tmp

write_lines(tmp$ensg, "~/Desktop/gz_over_cp.txt")


write_lines(mirna_df$mirna[which(mirna_df$category == "neurogenesis_gz")], "~/Desktop/mirna_gz_over_cp.txt")
write_lines(mirna_df$mirna[which(mirna_df$category == "neurogenesis_cp")], "~/Desktop/mirna_cp_over_gz.txt")

write_lines(mirna_df$mirna[which(mirna_df$category == "maturation_early")], "~/Desktop/mirna_mat_early.txt")
write_lines(mirna_df$mirna[which(mirna_df$category == "maturation_late")], "~/Desktop/mirna_mat_late.txt")


mirna_df %>%
  arrange(log2FoldChange.gzcp) -> tmp

write_lines(tmp$mirna_clean, "~/Desktop/mirna_ranked_cp_over_gz.txt")

mirna_df %>%
  arrange(-log2FoldChange.gzcp) -> tmp

write_lines(tmp$mirna_clean, "~/Desktop/mirna_ranked_gz_over_cp.txt")

