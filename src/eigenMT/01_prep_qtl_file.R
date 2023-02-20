# prep "QTL" files for eigenMT
# need 3 columns: 'snps', 'gene', and 'pvalue'
# using EMMAX nominal p-values
# one file per chrom

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES #########################################################################################################
# directory for qtl files
output.dir <- paste0(here("results/eigenmt/"), results.name, "/qtl_files/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# emmax results on all variants
results.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Data ##########################################################################################################

df.results <- as_tibble(readRDS(results.rds))


# Export Files #########################################################################################################

for (chrom in CHROMS) {

    df.output <- NULL

    # filter for this chr, select and relable columns
    df.results %>%
        dplyr::filter(CHR == chrom) %>%
        dplyr::select(snps = SNP,
                      gene = UniName,
                      pvalue = P) -> df.output

    # write tab sep file
    write_tsv(df.output,
              path = paste0(output.dir, chrom, ".qtl.txt"))

}

