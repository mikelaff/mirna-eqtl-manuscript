# build ranged summarized experiment for fetal tissue samples of RNA splicing using leafcutter
# Data analysis done by Nil

# Per sample splice percent
# Quantile normalized
# /proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/bulk_unique_perind.counts.gz.qqnorm_chr$i (fetal bulk)
# Raw ratios
# /proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/bulk_unique_perind.counts.gz.phen_chr$i (fetal bulk)

library(here)
library(readr)
library(dplyr)
library(magrittr)
#library(ensembldb)
#library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)

# OUTPUT FILES ####################################################################################
output.rse.rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_rse_leafcutter_splicing.rds", sep="")

# INPUT FILES #####################################################################################
# gene rse for phenotype information
gene.rse.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# Raw and normalized ratios dir
ratios.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/"

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import RSE for Gene Expression #######################################################################################
rse.genes <- read_rds(gene.rse.rds)

# Splicing Data ########################################################################################################

df.raw.ratios <- tibble()
df.norm.ratios <- tibble()


for (chr in CHROMS) {

    print(paste("Working on", chr))



    print("Importing...")

    if (chr == "chrX") {
        df.this.chrom.raw <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.phen_chr23"))

        df.this.chrom.norm <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.qqnorm_chr23"))

    } else {
        df.this.chrom.raw <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.phen_", chr))

        df.this.chrom.norm <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.qqnorm_", chr))

    }

    colnames(df.this.chrom.raw)[1] <- "chr"
    colnames(df.this.chrom.norm)[1] <- "chr"

    df.this.chrom.raw$chr <- chr
    df.this.chrom.norm$chr <- chr

    df.raw.ratios %<>%
        bind_rows(df.this.chrom.raw)

    df.norm.ratios %<>%
        bind_rows(df.this.chrom.norm)

}

rm(df.this.chrom.norm, df.this.chrom.raw)

colnames(df.raw.ratios) <- sapply(strsplit(colnames(df.raw.ratios), "\\."), `[`, 1)
colnames(df.norm.ratios) <- sapply(strsplit(colnames(df.norm.ratios), "\\."), `[`, 1)

all(df.raw.ratios$ID %in% df.norm.ratios$ID)
all(df.norm.ratios$ID %in% df.raw.ratios$ID)

df.raw.ratios %<>%
    arrange(ID)

df.norm.ratios %<>%
    arrange(ID)

all(df.raw.ratios$ID == df.norm.ratios$ID)

# Create RSE ######################################################################################
# create rowranges
gr.rows <- GRanges(seqnames = df.raw.ratios$chr,
                   ranges = IRanges(start = df.raw.ratios$start,
                                    end = df.raw.ratios$end,
                                    names = df.raw.ratios$ID),
                   strand = sapply(strsplit(df.raw.ratios$ID, "_"), `[`, 3))

# coldata
coldata <- colData(rse.genes)

# matrix
df.raw.ratios %>%
    dplyr::select(starts_with("RNAID")) %>%
    as.matrix() -> mat.raw
rownames(mat.raw) <- names(gr.rows)

df.norm.ratios %>%
    dplyr::select(starts_with("RNAID")) %>%
    as.matrix() -> mat.norm
rownames(mat.norm) <- names(gr.rows)

# order matrix

mode(mat.raw) <- "double"
mode(mat.norm) <- "double"

all(rownames(mat.raw) == rownames(mat.norm))
all(colnames(mat.raw) == colnames(mat.norm))

coldata <- coldata[colnames(mat.raw),]

stopifnot(all(rownames(coldata) == colnames(mat.raw)))
stopifnot(all(rownames(coldata) == colnames(mat.norm)))

# build RangedSummarizedExperiment
rse <- SummarizedExperiment(assays = list(splicing.raw = mat.raw, splicing.norm = mat.norm),
                            rowRanges = gr.rows,
                            colData = coldata)

# save rse to rds file
saveRDS(rse, output.rse.rds)
