# format eQTL data for garfield enrichment analysis
# write pval files
# one per chrom
# hg38 positions, pvalues

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(readxl)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################

dir.output <- here("data/garfield-data/pval/huan2015_onlySig/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################

# directory with maf and tss distance variants.
maftssd.dir <- here("data/garfield-data/maftssd/hg38-1000G-EUR/")

# Huan 2015 miRNA-eQTLs (all snps passing their FDR threshold, not clumped)
blood.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# chain file for hg19 to hg38 liftOver
#path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
path.hg19.to.hg38 <- here("data/hg19ToHg38.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

p.val.threshold <- 6.62e-5

# Import Blood miRNA-eQTLs #############################################################################################

df.eqtls <- read_xlsx(blood.eqtls.xlsx, skip = 1)

df.eqtls %<>%
    dplyr::select(SNP.pos,
                  chr.SNP,
                  Pval)

# GRanges from blood eqtls
gr.eqtls <- makeGRangesFromDataFrame(df = df.eqtls,
                                     start.field = "SNP.pos",
                                     end.field = "SNP.pos",
                                     seqnames.field = "chr.SNP",
                                     keep.extra.columns = TRUE)

# liftover to hg38
# chain for hg19 to hg38 conversion
ch <- import.chain(path.hg19.to.hg38)

# GRangesList of GRanges conversion
lo <- liftOver(gr.eqtls, ch)

# unlist to get GRanges
gr.eqtls.hg38 <- unlist(lo)

# modify seqlevels and seqinfo
seqlevels(gr.eqtls.hg38) <- CHROMS
gr.eqtls.hg38 <- keepSeqlevels(gr.eqtls.hg38,
                               CHROMS,
                               pruning.mode = "coarse")

rm(ch, lo, gr.eqtls)


# Format and Export ####################################################################################################

df.eqtls <- as_tibble(gr.eqtls.hg38)

df.eqtls %<>%
    dplyr::select(BP.hg38 = start,
                  P = Pval,
                  CHR = seqnames)

df.eqtls %<>%
    dplyr::mutate(pos = paste(CHR, BP.hg38, sep = "_"))

df.eqtls %>%
    group_by(pos) %>%
    filter(P == min(P)) -> df.eqtls.minP


df.eqtls.minP %>%
    ungroup() %>%
    dplyr::select(pos = BP.hg38,
                  pval = P,
                  chr = CHR) %>%
    distinct() -> df.hg38

for (chrom in CHROMS) {

    df.output <- dplyr::filter(df.hg38, chr == chrom)

    write_lines(paste(df.output$pos, signif(df.output$pval, digits = 6)), paste0(dir.output, chrom))

}

# # check against garfield prep file to see what happens with repeated locations
# df.results$pos <- paste(df.results$CHR, df.results$BP.hg38, sep = "_")
# sum(duplicated(df.results$pos))
#
# df.results %>%
#     group_by(pos) %>%
#     filter(P == min(P)) -> df.results.minP
#
# df.prep <- read_table2("~/Downloads/garfield.prep.hg38-mirQTL-MIXED_mirQTL_chromHMM.out",
#                        col_names = c("chr", "bp", "pval", "tags", "maf", "dist", "anno"),
#                        col_types = "cididic")
#
# df.prep.minP <- read_table2("~/Downloads/garfield.prep.hg38-mirQTL-MIXED_mirQTL-minP_chromHMM.out",
#                             col_names = c("chr", "bp", "pval", "tags", "maf", "dist", "anno"),
#                             col_types = "cididic")
#
# df.prep$pos <- paste(df.prep$chr, df.prep$bp, sep = "_")
# sum(duplicated(df.prep$pos))



