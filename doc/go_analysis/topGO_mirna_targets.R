# use topGO to look at mirna targets

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(readxl)
library(topGO)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/go_analysis/pdfs/")
# write plots
write_plots <- FALSE

# INPUT FILES #####################################################################################
# mirtarbase v7 hsa targets
mirna_targets_xlsx <- here("data/mirtarbase_v7/hsa_MTI.xlsx")
# combined deseq2 shrunken results as a dataframe with labels for significance and category (neuro vs maturation)
combined_df_rds <- here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_df.rds")

# Load DF #########################################################################################
combined_results_df <- readRDS(combined_df_rds)

# Load Targets ####################################################################################
targets <- read_xlsx(mirna_targets_xlsx)
# trim dataframe
targets %<>% dplyr::select(mirna = miRNA, target_gene = `Target Gene`, target_gene_entrez_id = `Target Gene (Entrez Gene ID)`)
# only keep unique mirna-target interactions (remove rows from multiple references)
targets %<>% dplyr::filter(!duplicated(paste(targets$mirna, targets$target_gene_entrez_id, sep="_")))

# all gene target entrez ids
all_entrez_ids <- as.character(targets$target_gene_entrez_id[!duplicated(targets$target_gene_entrez_id)])

# Neuro Up in CP #####################
combined_results_df %>%
  filter(category == "neurogenesis_cp", baseMean.gzcp > 1000) %>%
  dplyr::select(miRNA = mirna, baseMean.gzcp, log2FoldChange.gzcp, padj.gzcp, sig.both, category) %>%
  arrange(log2FoldChange.gzcp) %>%
  top_n(n=5, wt=-log2FoldChange.gzcp) -> tmp

all(tmp$miRNA %in% targets$mirna)

tmp_targets <- as.character(dplyr::filter(targets, mirna %in% tmp$miRNA)$target_gene_entrez_id)
tmp_targets <- tmp_targets[!duplicated(tmp_targets)]

#write_lines(tmp_targets, "~/Desktop/cp_mirna_targets.txt")

geneList <- factor(as.integer(all_entrez_ids %in% tmp_targets))
names(geneList) <- all_entrez_ids

myGOdata <- new("topGOdata",
                ontology = "BP",
                description = "my data",
                allGenes = geneList,
                nodeSize = 10,
                annotationFun = annFUN.org,
                mapping = "org.Hs.eg.db")

resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")

# Neuro Up in GZ #####################
combined_results_df %>%
  filter(category == "neurogenesis_gz", baseMean.gzcp > 1000) %>%
  dplyr::select(miRNA = mirna, baseMean.gzcp, log2FoldChange.gzcp, padj.gzcp, sig.both, category) %>%
  arrange(log2FoldChange.gzcp) %>%
  top_n(n=5, wt=log2FoldChange.gzcp) -> tmp

all(tmp$miRNA %in% targets$mirna)

tmp_targets <- as.character(dplyr::filter(targets, mirna %in% tmp$miRNA)$target_gene_entrez_id)
tmp_targets <- tmp_targets[!duplicated(tmp_targets)]

write_lines(tmp_targets, "~/Desktop/gz_mirna_targets.txt")

geneList <- factor(as.integer(all_entrez_ids %in% tmp_targets))
names(geneList) <- all_entrez_ids

myGOdata <- new("topGOdata",
                ontology = "BP",
                description = "my data",
                allGenes = geneList,
                nodeSize = 10,
                annotationFun = annFUN.org,
                mapping = "org.Hs.eg.db")

resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 20, orderBy = "classic_fisher")






# SCRATCH #########################################################################################
# one mirna:
mirna_one <- "hsa-miR-128-2-5p"

# target entrez ids associated with this mirna
my_target_entrez_ids <- as.character(dplyr::filter(targets, mirna == mirna_one)$target_gene_entrez_id)

geneList <- factor(as.integer(all_entrez_ids %in% my_target_entrez_ids))
names(geneList) <- all_entrez_ids

myGOdata <- new("topGOdata",
                ontology = "BP",
                description = "my data",
                allGenes = geneList,
                nodeSize = 10,
                annotationFun = annFUN.org,
                mapping = "org.Hs.eg.db")

resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 50, orderBy = "classic_fisher")

# one mirna:
mirna_one <- "hsa-miR-129-5p"

# target entrez ids associated with this mirna
my_target_entrez_ids <- as.character(dplyr::filter(targets, mirna == mirna_one)$target_gene_entrez_id)

geneList <- factor(as.integer(all_entrez_ids %in% my_target_entrez_ids))
names(geneList) <- all_entrez_ids

myGOdata <- new("topGOdata",
                ontology = "BP",
                description = "my data",
                allGenes = geneList,
                nodeSize = 10,
                annotationFun = annFUN.org,
                mapping = "org.Hs.eg.db")

resultsFisher <- runTest(myGOdata, algorithm = "classic", statistic = "fisher")

allRes <- GenTable(myGOdata, classic_fisher = resultsFisher, topNodes = 50, orderBy = "classic_fisher")

write_lines(my_target_entrez_ids, "~/Desktop/mir_target_entrez_ids.txt")
