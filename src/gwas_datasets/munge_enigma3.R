# format ENIGMA3 GWAS data

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
dir.global.output.granges.rds <- here("data/gwas_datasets/enigma3/ENIGMA3_Global/")
dir.regional.output.granges.rds <- here("data/gwas_datasets/enigma3/ENIGMA3_withGlobal/")

# INPUT FILES ##########################################################################################################

dir.global.raw.data <- here("data/gwas_datasets/enigma3/ENIGMA3_Global/download/")
dir.regional.raw.data <- here("data/gwas_datasets/enigma3/ENIGMA3_withGlobal/download/")

# chain file for hg19 to hg38 liftOver
path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

HG38.SEQINFO <- Seqinfo(genome = "hg38")


# Mean Surface Area ####################################################################################################
# Import Data ##########################################################################################################
sa.global.raw.data <- paste0(dir.global.raw.data, "ENIGMA3_mixed_se_wo_Mean_Full_SurfArea_20190429.txt.gz")

df.sa.global.raw <- read_delim(sa.global.raw.data, delim = " ")

# columns from educational attainment used as model for the rest
# RSID CHR POS A1.effect A2 EUR.freq Beta SE Pval

# rename columns for consistency
df.sa.global.raw %<>%
    dplyr::select(RSID = SNP,
                  CHR,
                  POS = BP,
                  A1.effect = A1,
                  A2,
                  BETA = BETA1,
                  SE,
                  Pval = P,
                  FREQ1,
                  N)

df.sa.global.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

levels(factor(df.sa.global.raw$CHR))
# modify to UCSC naming, this dataset does not have X
df.sa.global.raw %<>%
    mutate(CHR = paste0("chr", CHR))

#df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

all(levels(factor(df.sa.global.raw$CHR)) %in% CHROMS)

print("Converting to GRanges...")

# make GRanges on hg19
gr.sa.global.raw <- makeGRangesFromDataFrame(df.sa.global.raw,
                                             keep.extra.columns = TRUE,
                                             ignore.strand = TRUE,
                                             seqnames.field = "CHR",
                                             start.field = "POS",
                                             end.field = "POS")

print("Finished converting to GRanges.")

rm(df.sa.global.raw)

# LiftOver to hg38 #####################################################################################################
print("LiftOver to hg38...")

# chain for hg19 to hg38 conversion
ch <- import.chain(path.hg19.to.hg38)

# GRangesList of GRanges conversion
lo <- liftOver(gr.sa.global.raw, ch)

# unlist to get GRanges
gr.sa.global.raw.hg38 <- unlist(lo)

# modify seqlevels and seqinfo
seqlevels(gr.sa.global.raw.hg38) <- CHROMS
seqinfo(gr.sa.global.raw.hg38) <- HG38.SEQINFO
gr.sa.global.raw.hg38 <- keepSeqlevels(gr.sa.global.raw.hg38,
                                       CHROMS,
                                       pruning.mode = "coarse")

print("Finished liftOver.")

rm(ch, lo, gr.sa.global.raw)

# Save RDS #############################################################################################################
print("Saving RDS...")

output.granges.sa.global.rds <- paste0(dir.global.output.granges.rds, "enigma3.surface_area.global.hg38.GRanges.rds")

saveRDS(gr.sa.global.raw.hg38, output.granges.sa.global.rds)

print("Finished saving RDS.")

rm(gr.sa.global.raw.hg38)
# Mean Thickness #######################################################################################################
# Import Data ##########################################################################################################
th.global.raw.data <- paste0(dir.global.raw.data, "ENIGMA3_mixed_se_wo_Mean_Full_Thickness_20190429.txt.gz")

df.th.global.raw <- read_delim(th.global.raw.data, delim = " ")

# columns from educational attainment used as model for the rest
# RSID CHR POS A1.effect A2 EUR.freq Beta SE Pval

# rename columns for consistency
df.th.global.raw %<>%
    dplyr::select(RSID = SNP,
                  CHR,
                  POS = BP,
                  A1.effect = A1,
                  A2,
                  BETA = BETA1,
                  SE,
                  Pval = P,
                  FREQ1,
                  N)

df.th.global.raw %<>%
    mutate(A1.effect = toupper(A1.effect),
           A2 = toupper(A2))

levels(factor(df.th.global.raw$CHR))
# modify to UCSC naming, this dataset does not have X
df.th.global.raw %<>%
    mutate(CHR = paste0("chr", CHR))

#df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

all(levels(factor(df.th.global.raw$CHR)) %in% CHROMS)

print("Converting to GRanges...")

# make GRanges on hg19
gr.th.global.raw <- makeGRangesFromDataFrame(df.th.global.raw,
                                             keep.extra.columns = TRUE,
                                             ignore.strand = TRUE,
                                             seqnames.field = "CHR",
                                             start.field = "POS",
                                             end.field = "POS")

print("Finished converting to GRanges.")

rm(df.th.global.raw)

# LiftOver to hg38 #####################################################################################################
print("LiftOver to hg38...")

# chain for hg19 to hg38 conversion
ch <- import.chain(path.hg19.to.hg38)

# GRangesList of GRanges conversion
lo <- liftOver(gr.th.global.raw, ch)

# unlist to get GRanges
gr.th.global.raw.hg38 <- unlist(lo)

# modify seqlevels and seqinfo
seqlevels(gr.th.global.raw.hg38) <- CHROMS
seqinfo(gr.th.global.raw.hg38) <- HG38.SEQINFO
gr.th.global.raw.hg38 <- keepSeqlevels(gr.th.global.raw.hg38,
                                       CHROMS,
                                       pruning.mode = "coarse")

print("Finished liftOver.")

rm(ch, lo, gr.th.global.raw)

# Save RDS #############################################################################################################
print("Saving RDS...")

output.granges.th.global.rds <- paste0(dir.global.output.granges.rds, "enigma3.thickness.global.hg38.GRanges.rds")

saveRDS(gr.th.global.raw.hg38, output.granges.th.global.rds)

print("Finished saving RDS.")

rm(gr.th.global.raw.hg38)

# Regional #############################################################################################################

regional.files <- list.files(dir.regional.raw.data)

# loop over each file
for (file in regional.files) {
    print(file)

    region <- strsplit(file, "_")[[1]][6]

    # output file name
    output.file <- NULL
    if (grepl("surfavg", file)) {

        output.file <- paste0(dir.regional.output.granges.rds, "enigma3.surface_area.", region, ".hg38.GRanges.rds")

    } else if (grepl("thickavg", file)) {

        output.file <- paste0(dir.regional.output.granges.rds, "enigma3.thickness.", region, ".hg38.GRanges.rds")

    } else {
        stop("something wrong")
    }

    # Import Data ####################

    df.regional.raw <- read_delim(paste0(dir.regional.raw.data, file), delim = " ")

    # rename columns for consistency
    df.regional.raw %<>%
        dplyr::select(RSID = SNP,
                      CHR,
                      POS = BP,
                      A1.effect = A1,
                      A2,
                      BETA = BETA1,
                      SE,
                      Pval = P,
                      FREQ1,
                      N)

    df.regional.raw %<>%
        mutate(A1.effect = toupper(A1.effect),
               A2 = toupper(A2))

    # modify to UCSC naming, this dataset does not have X
    df.regional.raw %<>%
        mutate(CHR = paste0("chr", CHR))

    #df.raw$CHR[df.raw$CHR == "chr23"] <- "chrX"

    stopifnot(all(levels(factor(df.regional.raw$CHR)) %in% CHROMS))

    print("Converting to GRanges...")

    # make GRanges on hg19
    gr.regional.raw <- makeGRangesFromDataFrame(df.regional.raw,
                                                 keep.extra.columns = TRUE,
                                                 ignore.strand = TRUE,
                                                 seqnames.field = "CHR",
                                                 start.field = "POS",
                                                 end.field = "POS")

    print("Finished converting to GRanges.")

    rm(df.regional.raw)

    # LiftOver to hg38 ######################
    print("LiftOver to hg38...")

    # chain for hg19 to hg38 conversion
    ch <- import.chain(path.hg19.to.hg38)

    # GRangesList of GRanges conversion
    lo <- liftOver(gr.regional.raw, ch)

    # unlist to get GRanges
    gr.regional.raw.hg38 <- unlist(lo)

    # modify seqlevels and seqinfo
    seqlevels(gr.regional.raw.hg38) <- CHROMS
    seqinfo(gr.regional.raw.hg38) <- HG38.SEQINFO
    gr.regional.raw.hg38 <- keepSeqlevels(gr.regional.raw.hg38,
                                           CHROMS,
                                           pruning.mode = "coarse")

    print("Finished liftOver.")

    rm(ch, lo, gr.regional.raw)

    # Save RDS #############################################################################################################
    print("Saving RDS...")

    saveRDS(gr.regional.raw.hg38, output.file)

    print("Finished saving RDS.")

    rm(gr.regional.raw.hg38)

}





