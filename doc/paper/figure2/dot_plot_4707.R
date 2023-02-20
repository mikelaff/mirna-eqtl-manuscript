
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(mikelaffr)

# FIGURE ###############################################################################################################
# figure 2B: dot plot 4707 locus
output.pdf <- paste0(here("doc/paper/figure2/pdfs/"), "figure2B_4707_locus.pdf")

# OUTPUT FILES #########################################################################################################
# output directory for pdf files
#dir.pdf <- here("doc/paper/figure2/pdfs/")

# INPUT FILES ##########################################################################################################
# primary association results
primary.results.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# primary eqtls
primary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/primary/20200120_mirQTLor_primary_eQTLs_dataFrame.rds")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# nominal p-value from eigenMT-BH procedure
nom.p.value.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################


# Import Summarized Results ############################################################################################
# mirQTL eqtls
df.eqtls <- readRDS(primary.eqtls.rds)

# all variants mirQTL
df.results <- as_tibble(readRDS(primary.results.rds))
# filter for only 4707
df.results %<>%
    filter(UniName == "hsa-mir-4707_hsa-miR-4707-3p")

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.value.txt))

# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]


# Dot Plot #############################################################################################################






#ggsave(output.pdf, height = 2, width = 6.5, units = "in", useDingbats=FALSE)
