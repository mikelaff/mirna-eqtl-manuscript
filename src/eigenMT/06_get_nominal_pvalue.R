# compile eigenMT results
# apply BH and get nominal p-value

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES #########################################################################################################
# output directory
output.dir <- paste0(here("results/eigenmt/"), results.name, "/compiled/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# nominal p-value output
output.txt <- paste0(here("results/eigenmt/"), results.name, "/compiled/", results.name, "_eigenMT-BH_nomPvalue_5percent.txt")

# nominal p-value output (10% FDR)
output.high.txt <- paste0(here("results/eigenmt/"), results.name, "/compiled/", results.name, "_eigenMT-BH_nomPvalue_10percent.txt")

# INPUT FILES ##########################################################################################################
# output directory from eigenMT
output.files.dir <- here("results/eigenmt/20200120_mirQTLor/output/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

FDR.THRESHOLD <- 0.05

FDR.THRESHOLD.HIGH <- 0.1

# Import Results #######################################################################################################

df.results <- tibble()

for (chrom in CHROMS) {

    df.input <- read_tsv(paste0(output.files.dir, chrom))

    df.results <- bind_rows(df.results, df.input)

}

# calculate FDR 0.05
df.results %>%
    mutate(BF.adj = p.adjust(BF, method = "fdr")) %>%
    filter(BF.adj < FDR.THRESHOLD) %>%
    top_n(1, pvalue) %>%
    pull(pvalue) -> p.value.nominal.threshold

p.value.nominal.threshold <- p.value.nominal.threshold[1]

# calculate FDR 0.1
df.results %>%
    mutate(BF.adj = p.adjust(BF, method = "fdr")) %>%
    filter(BF.adj < FDR.THRESHOLD.HIGH) %>%
    top_n(1, pvalue) %>%
    pull(pvalue) -> p.value.nominal.threshold.high

p.value.nominal.threshold.high <- p.value.nominal.threshold.high[1]

# Export Files #########################################################################################################

write_lines(p.value.nominal.threshold, output.txt)
write_lines(p.value.nominal.threshold.high, output.high.txt)

