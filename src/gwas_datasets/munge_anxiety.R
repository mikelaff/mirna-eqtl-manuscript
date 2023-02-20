# format Anxiety GWAS data

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
cc.output.granges.rds <- here("data/gwas_datasets/anxiety/anxiety.case_control.hg38.GRanges.rds")
fs.output.granges.rds <- here("data/gwas_datasets/anxiety/anxiety.factor_score.hg38.GRanges.rds")

# INPUT FILES ##########################################################################################################

cc.raw.data <- here("data/gwas_datasets/anxiety/download/anxiety.meta.full.cc.tbl.gz")
fs.raw.data <- here("data/gwas_datasets/anxiety/download/anxiety.meta.full.fs.tbl.gz")

# chain file for hg19 to hg38 liftOver
path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

HG38.SEQINFO <- Seqinfo(genome = "hg38")

# Import Data ##########################################################################################################
print("Importing GWAS data...")

# Case Control
df.cc.raw <- read_tsv(cc.raw.data)

# columns from educational attainment used as model for the rest
# RSID CHR POS A1.effect A2 EUR.freq Beta SE Pval

# rename columns for consistency
df.cc.raw %<>%
    dplyr::select(RSID = SNPID,
                  CHR,
                  POS = BP,
                  A1.effect = Allele1,
                  A2 = Allele2,
                  Freq1,
                  Beta = Effect,
                  SE = StdErr,
                  Pval = P.value,
                  TotalN)

df.cc.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

levels(factor(df.cc.raw$CHR))
# modify to UCSC naming, this dataset does not have X
df.cc.raw %<>%
    mutate(CHR = paste0("chr", CHR))

#df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

all(levels(factor(df.cc.raw$CHR)) %in% CHROMS)

# Factor Score
df.fs.raw <- read_tsv(fs.raw.data)

# columns from educational attainment used as model for the rest
# RSID CHR POS A1.effect A2 EUR.freq Beta SE Pval

# rename columns for consistency
df.fs.raw %<>%
    dplyr::select(RSID = SNPID,
                  CHR,
                  POS = BP,
                  A1.effect = Allele1,
                  A2 = Allele2,
                  Freq1,
                  Beta = Effect,
                  SE = StdErr,
                  Pval = P.value,
                  TotalN)

df.fs.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

levels(factor(df.fs.raw$CHR))
# modify to UCSC naming, this dataset does not have X
df.fs.raw %<>%
    mutate(CHR = paste0("chr", CHR))

#df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

all(levels(factor(df.fs.raw$CHR)) %in% CHROMS)

print("Converting to GRanges...")

# make GRanges on hg19
gr.cc.raw <- makeGRangesFromDataFrame(df.cc.raw,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = "CHR",
                                      start.field = "POS",
                                      end.field = "POS")

gr.fs.raw <- makeGRangesFromDataFrame(df.fs.raw,
                                      keep.extra.columns = TRUE,
                                      ignore.strand = TRUE,
                                      seqnames.field = "CHR",
                                      start.field = "POS",
                                      end.field = "POS")

print("Finished converting to GRanges.")

rm(df.cc.raw, df.fs.raw)

# LiftOver to hg38 #####################################################################################################
print("LiftOver to hg38...")

# chain for hg19 to hg38 conversion
ch <- import.chain(path.hg19.to.hg38)

# GRangesList of GRanges conversion
lo.cc <- liftOver(gr.cc.raw, ch)
lo.fs <- liftOver(gr.fs.raw, ch)

# unlist to get GRanges
gr.cc.raw.hg38 <- unlist(lo.cc)
gr.fs.raw.hg38 <- unlist(lo.fs)

# modify seqlevels and seqinfo
seqlevels(gr.cc.raw.hg38) <- CHROMS
seqlevels(gr.fs.raw.hg38) <- CHROMS

seqinfo(gr.cc.raw.hg38) <- HG38.SEQINFO
seqinfo(gr.fs.raw.hg38) <- HG38.SEQINFO

gr.cc.raw.hg38 <- keepSeqlevels(gr.cc.raw.hg38,
                                CHROMS,
                                pruning.mode = "coarse")
gr.fs.raw.hg38 <- keepSeqlevels(gr.fs.raw.hg38,
                                CHROMS,
                                pruning.mode = "coarse")

print("Finished liftOver.")

rm(ch, lo.cc, lo.fs, gr.cc.raw, gr.fs.raw)

# Save RDS #############################################################################################################
print("Saving RDS...")

saveRDS(gr.cc.raw.hg38, cc.output.granges.rds)
saveRDS(gr.fs.raw.hg38, fs.output.granges.rds)

print("Finished saving RDS.")
