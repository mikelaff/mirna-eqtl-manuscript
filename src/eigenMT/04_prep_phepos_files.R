# prep "PHEPOS" files for eigenMT
# 'gene_id', 'chrom_probe', 's1', 's2'
# one file per chrom

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES #########################################################################################################
# directory for phepos files
output.dir <- paste0(here("results/eigenmt/"), results.name, "/phepos_files/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# list of tested genes
genes.txt <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_genes_hg38.txt")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


# Import Genes #########################################################################################################
df.genes <- read_tsv(genes.txt, col_names = c("gene_id", "chrom_probe", "s1", "s2"))

# Export Files #########################################################################################################

for (chrom in CHROMS) {

    df.output <- NULL

    df.genes %>%
        dplyr::filter(chrom_probe == chrom) -> df.output

    # write tab sep file
    write_tsv(df.output,
              path = paste0(output.dir, chrom, ".phepos.txt"))

}

