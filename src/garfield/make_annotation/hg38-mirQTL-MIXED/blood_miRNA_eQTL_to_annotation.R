# process Huan 2015 blood miRNA-eQTLs for use as garfield annotations

library(here)
library(readr)
library(dplyr)
library(readxl)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
output.dir <- here("data/garfield-data/annotation/hg38-mirQTL-MIXED/blood_miRNA-eQTL/")

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# directory with maf and tss distance variants. only these variants will be used to create annotation files
maftssd.dir <- here("data/garfield-data/maftssd/hg38-mirQTL-MIXED/")

# Huan 2015 miRNA-eQTLs (all snps passing their FDR threshold, not clumped)
blood.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# chain file for hg19 to hg38 liftOver
#path.hg19.to.hg38 <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
path.hg19.to.hg38 <- here("data/hg19ToHg38.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Blood miRNA-eQTLs #############################################################################################

df.eqtls <- read_xlsx(blood.eqtls.xlsx, skip = 1)

# GRanges from blood eqtls
gr.eqtls <- makeGRangesFromDataFrame(df = df.eqtls,
                                     start.field = "SNP.pos",
                                     end.field = "SNP.pos",
                                     seqnames.field = "chr.SNP")

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

# Format Annotations ###################################################################################################

# annotation names
df.anno <- tibble(Annotation = "Blood_miRNA-eQTL")
gr.eqtls.hg38$name <- "Blood_miRNA-eQTL"

# add additional columns
df.anno %<>%
    mutate(Index = seq(0, nrow(df.anno)-1),
           Celltype = "PrimaryTissue",
           Tissue = "Blood",
           Type = "eQTL",
           Category = "miRNA-eQTL")

# select columns for link_file.txt
df.anno %<>%
    dplyr::select(Index,
                  Annotation,
                  Celltype,
                  Tissue,
                  Type,
                  Category)

# Process Variants #####################################################################################################

# loop over chroms
for (chrom in CHROMS) {
    print(paste("Processing variants from chromosome:", chrom))

    # load mafftssd file
    print("Loading variants...")
    df.vars <- read_table2(paste0(maftssd.dir, chrom), col_names = c("BP", "MAF", "TSSD"), col_types = "idi")
    df.vars$seqnames <- chrom

    df.vars %<>%
        dplyr::select(seqnames,
                      BP)

    # convert to GRanges
    gr.vars <- makeGRangesFromDataFrame(df.vars,
                                        ignore.strand = TRUE,
                                        seqnames.field = "seqnames",
                                        start.field = "BP",
                                        end.field = "BP")

    print("Finding overlaps...")
    # overlap with annotations
    hits <- GenomicRanges::findOverlaps(gr.vars, gr.eqtls.hg38, select = "first", ignore.strand = TRUE)
    df.vars$annotation <- gr.eqtls.hg38$name[hits]

    # get annotation index
    df.vars$anno_index <- df.anno$Index[match(df.vars$annotation, df.anno$Annotation)]

    # make -1 index for NA values
    df.vars$anno_index[is.na(df.vars$anno_index)] <- -1

    print("Formatting annotation strings...")
    # create annotation string
    df.vars$anno_string <- paste(rep("0", max(df.anno$Index) + 1), collapse = "")
    # replace 0 with 1 at anno index, -1's will replace at position 0, meaning no replacement
    substr(df.vars$anno_string, df.vars$anno_index + 1, df.vars$anno_index + 1) <- "1"

    print("Writing annotation file...")
    # export position and annotation string
    df.vars %<>%
        dplyr::select(BP, anno_string)
    write_delim(df.vars, path = paste0(output.dir, chrom), col_names = FALSE, delim = " ")

}

# Export Link File #####################################################################################################

print("Writing link_file.txt")
write_delim(df.anno, path = paste0(output.dir, "link_file.txt"), delim = " ")
