# format mirQTL data for garfield enrichment analysis
# write pval files
# one per chrom
# hg38 positions, pvalues

library(here)
library(dplyr)
library(magrittr)
library(readr)
#library(GenomicRanges)
#library(rtracklayer)
#library(liftOver)

# OUTPUT FILES #########################################################################################################

dir.output <- here("data/garfield-data/pval/mirQTL_minP/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################

# mirQTL all variants
mirqtl.results.variants.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# chain file for hg38 to hg19 liftOver
#path.hg38.to.hg19 <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Data ##########################################################################################################

df.results <- as_tibble(readRDS(mirqtl.results.variants.dataframe.rds))

df.results$pos <- paste(df.results$CHR, df.results$BP.hg38, sep = "_")
#sum(duplicated(df.results$pos))

df.results %>%
    group_by(pos) %>%
    filter(P == min(P)) -> df.results.minP

#sum(duplicated(df.results.minP$pos))

# gr.results <-  makeGRangesFromDataFrame(df.results,
#                                         keep.extra.columns = TRUE,
#                                         ignore.strand = TRUE,
#                                         seqnames.field = "CHR",
#                                         start.field = "BP.hg38",
#                                         end.field = "BP.hg38")

# LiftOver #############################################################################################################

# # chain for hg38 to hg19 conversion
# ch <- import.chain(path.hg38.to.hg19)
#
# # GRangesList of GRanges conversion
# lo <- liftOver(gr.results, ch)
#
# # unlist to get GRanges
# gr.hg19 <- unlist(lo)
#
# # modify seqlevels and seqinfo
# gr.hg19 <- keepSeqlevels(gr.hg19,
#                          CHROMS,
#                          pruning.mode = "coarse")
#
# df.hg19 <- as_tibble(gr.hg19)
#
# rm(ch, df.results, gr.results, lo, gr.hg19)

# Write Output #########################################################################################################

# df.hg19 %<>%
#     dplyr::select(chr = seqnames,
#                   pos = start,
#                   pval = P)


# for (chrom in CHROMS) {
#
#     df.output <- dplyr::filter(df.hg19, chr == chrom)
#
#     write_lines(paste(df.output$pos, signif(df.output$pval, digits = 6)), paste0(dir.output, chrom))
#
# }


df.results.minP %>%
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



