
library(here)
library(dplyr)
library(magrittr)
#library(plyranges)
library(DESeq2)
library(ggplot2)
library(rtracklayer)
library(EnsDb.Hsapiens.v86)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(mikelaffr)

## GLOBALS
# pdf directory
dir.pdf <- here("doc/eda/pdfs/")

# RangedSummarizedExperiment of miRge quantified miRNA expression data
mirna.rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GRanges of compiled known and novel miRNAs (miRBase_v22, Friedlander2014, Nowakowski2018, miRge, and miRDeep2)
mirna.granges.rds <- here("data/gtf_and_granges/20190917_all_known_and_novel_mirna_non_overlapping_granges.rds")

CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrY")

## Load Data

# miRge quantified miRNA expression data
rse <- readRDS(mirna.rse.rds)

# all known and novel miRNA annotations
gr.mirnas <- readRDS(mirna.granges.rds)
gr.mirnas <- keepSeqlevels(gr.mirnas, CHROMS, pruning.mode = "coarse")

# ENSEMBL gene database
gr.ensembl <- genes(EnsDb.Hsapiens.v86)
seqlevelsStyle(gr.ensembl) <- seqlevelsStyle(gr.mirnas)
gr.ensembl <- keepSeqlevels(gr.ensembl, CHROMS, pruning.mode = "coarse")
seqinfo(gr.ensembl) <- seqinfo(gr.mirnas)


## Summarize Data

#gr.mirnas %>%
#    plyranges::filter(type == "miRNA")

df.mirnas <- as.data.frame(gr.mirnas)

df.mirnas$mature <- ifelse(df.mirnas$type == "miRNA" | df.mirnas$type == "miRNA_putative_mature" | df.mirnas$type == "miRNA_putative_star", TRUE, FALSE)

df.mirnas %>%
    dplyr::filter(mature) %>%
    group_by(source) %>%
    add_tally() %>%
    ggplot(aes(x="", fill = reorder(source, n))) +
    geom_bar() +
    coord_polar("y", start = 0)

## Genomic Overlaps

### miRBase v22 miRNAs

# only miRBase miRNAs
gr.mirnas.mirbase <- gr.mirnas[gr.mirnas$source == "miRBase_v22"]

# remove miRNAs from ENSEMBL annotations
gr.ensembl.sub <- gr.ensembl[gr.ensembl$gene_biotype != "miRNA"]

# overlaps of known mirnas with any annotation except mirna
sum(!duplicated(queryHits(findOverlaps(gr.mirnas.mirbase, gr.ensembl.sub))))

# protein coding genes
gr.ensembl.prot.coding <- gr.ensembl[gr.ensembl$gene_biotype == "protein_coding"]

# overlaps of known mirnas with protein coding genes
sum(!duplicated(queryHits(findOverlaps(gr.mirnas.mirbase, gr.ensembl.prot.coding))))



df.mirnas %>%
    group_by(Name) %>%
    add_tally() -> df.mod

df.mirnas %>%
    dplyr::filter(Name == "hsa-miR-1179")

# expressed mirnas
# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

mirna.names.expressed <- rowRanges(rse)[rowRanges(rse)$Name %in% df.mirnas$Name]$Name

#pdf(file.path(dir.pdf, "mirna_from_multiple_loci_all_miRBase.pdf"))

df.mirnas %>%
    dplyr::filter(source == "miRBase_v22", type == "miRNA") %>%
    group_by(Name) %>%
    summarise(count = n()) %>%
    group_by(count) %>%
    summarise(n = n()) %>%
    ggplot(aes(x=count, y=n)) +
    geom_col() +
    geom_text(aes(y=n+10,label=n), vjust=0, size=6, fontface="bold") +
    scale_y_continuous(expand=expand_scale(mult = c(0.01,0.08))) +
    scale_x_continuous(breaks = c(1,2,3,4,5,11)) +
    plotTheme("presentation") +
    labs(y="Number of miRNAs",
         x="Number of Primary Sequences",
         title="miRBase miRNAs")

ggsave(paste(dir.pdf, "mirna_mult_loci_mirbase.pdf", sep = "/"), height = 6, width = 6)

#dev.off()

#pdf(file.path(dir.pdf, "mirna_from_multiple_loci_expressed.pdf"))

df.mirnas %>%
    dplyr::filter(Name %in% mirna.names.expressed) %>%
    group_by(Name) %>%
    summarise(count = n()) %>%
    group_by(count) %>%
    summarise(n = n()) %>%
    ggplot(aes(x=count, y=n)) +
    geom_col() +
    geom_text(aes(y=n+10,label=n), vjust=0, size=6, fontface="bold") +
    scale_y_continuous(expand=expand_scale(mult = c(0.01,0.08))) +
    scale_x_continuous(breaks = c(1,2,3,4,5,11)) +
    plotTheme("presentation") +
    labs(y="Number of miRNAs",
         x="Number of Primary Sequences",
         title="Expressed miRNAs")

ggsave(paste(dir.pdf, "mirna_mult_loci_expressed.pdf", sep = "/"), height = 6, width = 6)

#dev.off()

# check each mature sequence for number of genomic loci

df.mirnas %>%
    dplyr::filter(type == "miRNA_putative_mature" | type == "miRNA_putative_star", Name %in% mirna.names.expressed) %>%
    group_by(Name) %>%
    add_tally() -> df.mirnas.novel.mature.expressed

df.mirnas.novel.mature.expressed$num_loci <- 0

mirna.dna.string.set <- DNAStringSet(RNAStringSet(df.mirnas.novel.mature.expressed$sequence))

genome <- BSgenome.Hsapiens.UCSC.hg38

#matchPattern(DNAString(RNAString(df.mirnas.mature$sequence[1])), genome[["chr1"]])

# loop over chroms
for (chr in CHROMS) {
    print(chr)

    chrom <- genome[[chr]]

    # loop over all mirna
    for (i in 1:length(mirna.dna.string.set)) {
        print(paste0(i, " of ", length(mirna.dna.string.set), ": ", df.mirnas.novel.mature.expressed$Name[i]))
        df.mirnas.novel.mature.expressed$num_loci[i] <- df.mirnas.novel.mature.expressed$num_loci[i] + countPattern(mirna.dna.string.set[[i]], chrom)
    }

}


# for (i in 1:length(df.mirnas.mature$ID)) {
#     print(paste0(i, " of ", length(df.mirnas.mature$ID), ": ", df.mirnas.mature$Name[i]))
#     numMatches <- 0
#     # loop over all chroms
#     for (k in 1:length(CHROMS)) {
#         print(CHROMS[k])
#         numMatches <- numMatches + countPattern(DNAString(RNAString(df.mirnas.mature$sequence[i])), genome[[CHROMS[k]]])
#     }
#     print(paste0("Total matches for ", df.mirnas.mature$Name[i], ": ", numMatches))
#     df.mirnas.mature$num_loci[i] <- numMatches
# }





