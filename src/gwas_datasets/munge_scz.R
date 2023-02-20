# format SCZ GWAS data

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
output.granges.rds <- here("data/gwas_datasets/scz/scz.hg38.GRanges.rds")

# INPUT FILES ##########################################################################################################

raw.data <- here("data/gwas_datasets/scz/download/clozuk_pgc2.meta.sumstats.txt.gz")

# chain file for hg19 to hg38 liftOver
path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")

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
    dplyr::select(SNPID = SNP,
                  CHR,
                  POS = BP,
                  A1.effect = A1,
                  A2,
                  OR,
                  SE,
                  Pval = P)

df.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

levels(factor(df.raw$CHR))

# modify to UCSC naming, this dataset does have X
df.raw %<>%
    mutate(CHR = paste0("chr", CHR))

df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

# remove NA positions
# df.raw %<>%
#     filter(CHR != "chrNA")

all(levels(factor(df.raw$CHR)) %in% CHROMS)

print("Converting to GRanges...")

# make GRanges on hg19
gr.raw <- makeGRangesFromDataFrame(df.raw,
                                   keep.extra.columns = TRUE,
                                   ignore.strand = TRUE,
                                   seqnames.field = "CHR",
                                   start.field = "POS",
                                   end.field = "POS")

print("Finished converting to GRanges.")

rm(df.raw)

# LiftOver to hg38 #####################################################################################################
print("LiftOver to hg38...")

# chain for hg19 to hg38 conversion
ch <- import.chain(path.hg19.to.hg38)

# GRangesList of GRanges conversion
lo <- liftOver(gr.raw, ch)

# unlist to get GRanges
gr.raw.hg38 <- unlist(lo)

# modify seqlevels and seqinfo
seqlevels(gr.raw.hg38) <- CHROMS
seqinfo(gr.raw.hg38) <- HG38.SEQINFO
gr.raw.hg38 <- keepSeqlevels(gr.raw.hg38,
                             CHROMS,
                             pruning.mode = "coarse")

print("Finished liftOver.")

rm(ch, lo, gr.raw)

# Save RDS #############################################################################################################
print("Saving RDS...")

saveRDS(gr.raw.hg38, output.granges.rds)

print("Finished saving RDS.")

rm(gr.raw.hg38)
