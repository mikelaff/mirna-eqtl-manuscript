# investigate cell-type specific gene lengths and 3' UTR lengths
# question of if neuron cell-type specific genes have longer lengths and longer 3' UTRs
# is this the reason for increased miRNA target enrichment among neuronal cell types?
# 28 May 2019

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(reshape2)

library(EnsDb.Hsapiens.v86)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ########################################################################################################
# output directory for png files
dir.png <- here("doc/eda/pngs/")

# INPUT FILES #########################################################################################################
# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.xlsx <- here("data/ucla_single_cell/TableS4 Cluster analysis.xlsx")

# GLOBALS #############################################################################################################
# write plots
WRITE.PLOTS <- TRUE

# Import Single Cell Gene Data ########################################################################################

df.genes.by.cluster <- read_xlsx(genes.by.cluster.xlsx, sheet = "Cluster enriched genes")
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

# Annotate Genes ######################################################################################################
# get biotypes by ensembl id
#df.genes <- tidyr::unnest(genes(EnsDb.Hsapiens.v86, columns = c("gene_biotype", "entrezid"), return.type = "data.frame"))
df.genes <- genes(EnsDb.Hsapiens.v86, columns = c("gene_biotype"), return.type = "data.frame")
df.genes %<>% dplyr::rename(Ensembl = gene_id)
# join with genes by cluster df
df.genes.by.cluster %<>%
  left_join(df.genes, by = "Ensembl")
rm(df.genes)

# Get Gene Lengths ####################################################################################################
# get transcript lengths, cds, 5'utr, and 3'utr lenghts
df.transcripts <- transcriptLengths(EnsDb.Hsapiens.v86,
                                    with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE)
# summarize lenghts by gene
df.transcripts %<>%
  dplyr::filter(gene_id %in% df.genes.by.cluster$Ensembl) %>%
  dplyr::group_by(gene_id) %>%
  dplyr::summarise_at(.vars = c("nexon", "tx_len", "cds_len", "utr5_len", "utr3_len"),
                      .funs = c(mean=mean, median=median, min=min, max=max, var=var)) %>%
  dplyr::rename(Ensembl = gene_id)

# join with cluster genes
df.genes.by.cluster %<>%
  left_join(df.transcripts, by = "Ensembl")

rm(df.transcripts)

# Plots ###############################################################################################################

# df for cluster labels
CellType <- factor(c("Other", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "ExNeuron", "InNeuron", "InNeuron",
              "Progenitor", "Other", "Other", "Progenitor", "Other", "Progenitor", "Progenitor", "Progenitor"),
              levels = c("Progenitor", "ExNeuron", "InNeuron", "Other"), ordered = TRUE)
df.cluster <- data.frame(CellType = CellType)
df.cluster$Cluster <- as.factor(clusters)

df.genes.by.cluster %<>%
  left_join(df.cluster, by = "Cluster")
# reorder clusters for plotting
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster,
                                      levels = df.cluster$Cluster[order(df.cluster$CellType)],
                                      ordered = TRUE)

df.genes.by.cluster %>%
  dplyr::select(Ensembl, Cluster, gene_biotype, CellType,
                starts_with("nexon_"),
                starts_with("tx_"),
                starts_with("cds_"),
                starts_with("utr5_"),
                starts_with("utr3_")) %>%
  melt(id.vars = c("Ensembl", "Cluster", "gene_biotype", "CellType")) %>%
  dplyr::filter(!grepl("_var", variable)) %>%
  dplyr::filter(!grepl("nexon", variable)) %>%
  dplyr::filter(grepl("mean", variable)) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(value != 0) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  ggplot(aes(x=Cluster, y=value, fill=CellType)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9)) +
  labs(y="Log10(length in bp)",
       title="Mean Lengths by Cell Type (protein-coding only)") +
  scale_fill_manual(values = c(Other = "#1B9E77",
                               ExNeuron = "#D95F02",
                               InNeuron = "#7570B3",
                               Progenitor = "#E7298A"))

if(WRITE.PLOTS){ggsave(file.path(dir.png,"mean_lengths_by_cell_type_protein-coding.png"), width = 8, height = 7)}

df.genes.by.cluster %>%
  dplyr::select(Ensembl, Cluster, gene_biotype, CellType,
                starts_with("nexon_"),
                starts_with("tx_"),
                starts_with("cds_"),
                starts_with("utr5_"),
                starts_with("utr3_")) %>%
  melt(id.vars = c("Ensembl", "Cluster", "gene_biotype", "CellType")) %>%
  dplyr::filter(!grepl("_var", variable)) %>%
  dplyr::filter(!grepl("nexon", variable)) %>%
  dplyr::filter(grepl("median", variable)) %>%
  dplyr::filter(!is.na(value)) %>%
  dplyr::filter(value != 0) %>%
  dplyr::filter(gene_biotype == "protein_coding") %>%
  ggplot(aes(x=Cluster, y=value, fill=CellType)) +
  geom_boxplot() +
  facet_wrap(~variable) +
  scale_y_log10() +
  theme(axis.text.x = element_text(angle = 45, vjust = 0.9)) +
  labs(y="Log10(length in bp)",
       title="Median Lengths by Cell Type (protein-coding only)") +
  scale_fill_manual(values = c(Other = "#1B9E77",
                               ExNeuron = "#D95F02",
                               InNeuron = "#7570B3",
                               Progenitor = "#E7298A"))

if(WRITE.PLOTS){ggsave(file.path(dir.png,"median_lengths_by_cell_type_protein-coding.png"), width = 8, height = 7)}
