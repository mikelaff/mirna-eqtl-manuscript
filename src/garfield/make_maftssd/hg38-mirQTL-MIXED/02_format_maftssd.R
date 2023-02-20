# convert a plink .frq to files suitable to garfield maftssd

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(biomaRt)
library(GenomicRanges)

# OUTPUT FILES #########################################################################################################
output.dir <- here("data/garfield-data/maftssd/hg38-mirQTL-MIXED/")

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
input.dir <- "/pine/scr/m/i/mikelaff/garfield-mirQTL-MIXED/"

#input.dir <- "~/Downloads/"

#gencode.gtf <- here("data/genome/hg38_UCSC/gencode.v25.primary_assembly.annotation.gtf")

# GLOBALS ##############################################################################################################
CHROMS <- c(seq(1,22), "X")

# Get TSS ##############################################################################################################
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

df.tss <- as_tibble(getBM(attributes = c("transcription_start_site",
                                         "chromosome_name",
                                         "strand",
                                         "ensembl_transcript_id",
                                         "transcript_gencode_basic"),
                          filters = "transcript_gencode_basic", values = c(TRUE),
                          mart = genemart))

df.tss %<>%
    dplyr::filter(chromosome_name %in% CHROMS)

df.tss$strand[df.tss$strand == 1] <- "+"
df.tss$strand[df.tss$strand == -1] <- "-"

gr.tss <- makeGRangesFromDataFrame(df.tss,
                                   seqnames.field = "chromosome_name",
                                   start.field = "transcription_start_site",
                                   end.field = "transcription_start_site",
                                   strand.field = "strand")

rm(genemart, df.tss)

# Import/Export ########################################################################################################
# Loop over each chr
for (chrom in CHROMS) {

    df.maf <- NULL
    gr.nearest <- NULL
    gr.snp <- NULL

    print(paste("Working on chrom:", chrom))

    # import .frq
    print("Importing...")
    df.maf <- read_table2(paste0(input.dir, "chr", chrom, "_updated_var_ids_maf.frq"), col_types = "ccccdi")

    # filter for duplicates remove excess columns
    print("Filtering...")
    df.maf %>%
        dplyr::select(SNP, MAF) %>%
        dplyr::filter(!duplicated(SNP)) -> df.maf

    # remove chr and :
    print("Formatting...")
    df.maf$SNP <- as.integer(gsub(paste0(chrom, ":"), "", df.maf$SNP))

    print("Getting TSS distance...")
    # GRanges of all the SNPs
    gr.snp <- GRanges(seqnames = chrom,
                      ranges = IRanges(start = df.maf$SNP, end = df.maf$SNP),
                      strand = "*")

    # GRanges of nearest tss
    gr.nearest <- gr.tss[nearest(gr.snp, gr.tss)]

    # Distance to nearest tss (upstream of tss: negative, downstream of tss: positive)
    df.maf$TSSD <- start(gr.snp) - start(gr.nearest)

    # export formated tags file
    print("Exporting...")
    write_delim(df.maf, paste0(output.dir, "chr", chrom), delim = " ", col_names = FALSE)

}
