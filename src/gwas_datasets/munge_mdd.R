# format MDD GWAS data

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
output.granges.rds <- here("data/gwas_datasets/mdd/mdd.hg38.GRanges.rds")

# INPUT FILES ##########################################################################################################

raw.data <- here("data/gwas_datasets/mdd/download/PGC_UKB_depression_genome-wide.txt.gz")

# chain file for hg19 to hg38 liftOver
path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")

# use parkinson's disorder munged data to get genomic positions for each rsid
gr.pd.rds <- here("data/gwas_datasets/pd/pd.hg38.GRanges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

HG38.SEQINFO <- Seqinfo(genome = "hg38")

# Import Data ##########################################################################################################
print("Importing GWAS data...")

df.raw <- read_delim(raw.data, delim = " ")

# columns from educational attainment used as model for the rest
# RSID CHR POS A1.effect A2 EUR.freq Beta SE Pval

# rename columns for consistency
df.raw %<>%
    dplyr::select(RSID = MarkerName,
                  A1.effect = A1,
                  A2,
                  LogOR,
                  SE.LogOR = StdErrLogOR,
                  Pval = P)


# import pd data
gr.pd <- readRDS(gr.pd.rds)

df.pd <- as_tibble(gr.pd)
df.pd %<>%
    dplyr::select(RSID,
                  CHR = seqnames,
                  POS = start)

df.raw %<>%
    left_join(df.pd, by = "RSID")

df.raw %<>%
    filter(!is.na(CHR))

df.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

print("Converting to GRanges...")

# make GRanges on hg19
gr.raw <- makeGRangesFromDataFrame(df.raw,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE,
                                   seqnames.field = "CHR",
                                   start.field = "POS",
                                   end.field = "POS")

# modify seqlevels and seqinfo
seqlevels(gr.raw) <- CHROMS
seqinfo(gr.raw) <- HG38.SEQINFO
gr.raw <- keepSeqlevels(gr.raw,
                        CHROMS,
                        pruning.mode = "coarse")

print("Finished converting to GRanges.")

rm(df.raw, df.pd, gr.pd)

# Save RDS #############################################################################################################
print("Saving RDS...")

saveRDS(gr.raw, output.granges.rds)

print("Finished saving RDS.")

rm(gr.raw)
