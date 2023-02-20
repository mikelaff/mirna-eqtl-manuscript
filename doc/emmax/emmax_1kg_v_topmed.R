# compare 1kg v topmed emmax association test

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(rtracklayer)
library(liftOver)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/1kg_v_topmed/")
#dir.pngs <- here("doc/emmax/pngs/")

# INPUT FILES ##########################################################################################################
# GRangesList containing SNP-miR associations within 1MB of each miR
grl.1kg.rds <- here("results/emmax/association_results/20191101_mirQTL_1kg/compiled/20191101_mirQTL_1kg_GRangesList.rds")
grl.topmed.rds <- here("results/emmax/association_results/20191101_mirQTL_topmed/compiled/20191101_mirQTL_topmed_GRangesList.rds")
# Dataframe containing the association results
df.1kg.rds <- here("results/emmax/association_results/20191101_mirQTL_1kg/compiled/20191101_mirQTL_1kg_dataFrame.rds")
df.topmed.rds <- here("results/emmax/association_results/20191101_mirQTL_topmed/compiled/20191101_mirQTL_topmed_dataFrame.rds")

# genes (mirs) used in this analysis
genes.1kg.hg19.txt <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/20191101_mirQTL_genes_hg19.txt")
genes.topmed.hg38.txt <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/20191101_mirQTL_genes_hg38.txt")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

FDR.THRESHOLD <- 0.05
# Import Data ##########################################################################################################

# df.1kg <- readRDS(df.1kg.rds)
# df.topmed <- readRDS(df.topmed.rds)
#
# df.1kg %<>% as_tibble()
# df.topmed %<>% as_tibble()

# 1kg GRangesList
grl.1kg <- readRDS(grl.1kg.rds)

#seqlevelsStyle(grl.1kg)
seqlevels(grl.1kg) <- CHROMS
seqinfo(grl.1kg) <- Seqinfo(genome = "hg19")

# 1kg GRanges
gr.1kg <- unlist(grl.1kg)

# chain for hg19 to hg38 conversion
path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch <- import.chain(path)

# GRangesList of GRanges conversion
lo <- liftOver(gr.1kg, ch)

# keep only GRanges that successfully lifted over, convert to GRanges
gr.1kg.hg38 <- unlist(lo[lengths(lo) == 1])
names(gr.1kg.hg38) <- paste(gr.1kg.hg38$UniName, gr.1kg.hg38$SNP, sep = ".")

rm(lo, gr.1kg, grl.1kg, ch)

df.1kg.hg38 <- as_tibble(gr.1kg.hg38)
rm(gr.1kg.hg38)

# import topmed results
grl.topmed.hg38 <- readRDS(grl.topmed.rds)
gr.topmed.hg38 <- unlist(grl.topmed.hg38)

df.topmed.hg38 <- as_tibble(gr.topmed.hg38)
rm(grl.topmed.hg38, gr.topmed.hg38)

# Format ? #########

mirs.common <- unique(df.1kg.hg38$UniName)[unique(df.1kg.hg38$UniName) %in% unique(df.topmed.hg38$UniName)]


df.1kg.hg38 %>%
    filter(A1.HOM.count != 1, HET.count >= 2) -> df.1kg.hg38.filt

df.topmed.hg38 %>%
    filter(A1.HOM.count != 1, HET.count >= 2) -> df.topmed.hg38.filt

df.1kg.hg38.filt$PVAL.adj <- p.adjust(df.1kg.hg38.filt$PVAL, method = "fdr")
df.topmed.hg38.filt$PVAL.adj <- p.adjust(df.topmed.hg38.filt$PVAL, method = "fdr")

# summarize mirs by minimum adjusted p-value, threshold to find emiRs
df.1kg.hg38.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.1kg

df.topmed.hg38.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.topmed

df.1kg.hg38.filt %>%
    filter(PVAL.adj < FDR.THRESHOLD) %>%
    top_n(1, PVAL) %>%
    pull(PVAL) -> p.value.nominal.threshold.1kg

df.topmed.hg38.filt %>%
    filter(PVAL.adj < FDR.THRESHOLD) %>%
    top_n(1, PVAL) %>%
    pull(PVAL) -> p.value.nominal.threshold.topmed

df.1kg.hg38.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < 0.1) %>%
    pull(UniName) -> emiRs.1kg.high.fdr

df.topmed.hg38.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < 0.1) %>%
    pull(UniName) -> emiRs.topmed.high.fdr

# Plot ? ##########

emiRs.1kg.only <- emiRs.1kg[!emiRs.1kg %in% emiRs.topmed]
emiRs.topmed.only <- emiRs.topmed[!emiRs.topmed %in% emiRs.1kg]

emiRs.1kg.only.high.fdr <- emiRs.1kg.high.fdr[!emiRs.1kg.high.fdr %in% emiRs.topmed.high.fdr]
emiRs.topmed.only.high.fdr <- emiRs.topmed.high.fdr[!emiRs.topmed.high.fdr %in% emiRs.1kg.high.fdr]

emiRs.combined <- unique(c(emiRs.1kg, emiRs.topmed))
emiRs.combined.high.fdr <- unique(c(emiRs.1kg.high.fdr, emiRs.topmed.high.fdr))

emiRs.combined.only.high <- emiRs.combined.high.fdr[!emiRs.combined.high.fdr %in% emiRs.combined]

df.genes <- read_tsv(genes.topmed.hg38.txt, col_names = c("uniqename", "chr", "start", "end"))
df.emirs <- dplyr::filter(df.genes, uniqename %in% emiRs.combined)

#df.emirs <- dplyr::filter(df.genes, uniqename %in% emiRs.combined.only.high)

# loop over each emir

genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

for (i in 1:length(df.emirs$uniqename)) {
#for (i in 294:306) {
    cat(df.emirs$uniqename[i], "\n")

    mir <- df.emirs$uniqename[i]

    if (mir %in% emiRs.1kg.only) {
        filename <- file.path(dir.pdfs, paste0("/1kg_only/", mir, ".pdf"))
    } else if (mir %in% emiRs.topmed.only) {
        filename <- file.path(dir.pdfs, paste0("/topmed_only/", mir, ".pdf"))
    } else {
        filename <- file.path(dir.pdfs, paste0("/both/", mir, ".pdf"))
    }


    MB.WINDOW <- 1e6

    startRange <- NULL
    endRange <- NULL

    chr <- df.emirs$chr[i]

    # genomic range for SNPs associated with this miRNA
    startRange <- df.emirs$start[i] - MB.WINDOW
    endRange <- df.emirs$end[i] + MB.WINDOW

    # 1kg Data GRanges ##############################

    df.1kg.hg38.filt %>%
        filter(UniName == !!mir) -> df.1kg.sub

    gr.1kg <- GRanges(seqnames = df.1kg.sub$seqnames,
                      ranges = IRanges(start = df.1kg.sub$start,
                                       end = df.1kg.sub$end,
                                       names = df.1kg.sub$SNP),
                      neglog10p = -log10(df.1kg.sub$PVAL))

    # topmed Data GRanges #####################################

    df.topmed.hg38.filt %>%
        filter(UniName == !!mir) -> df.topmed.sub

    gr.topmed <- GRanges(seqnames = df.topmed.sub$seqnames,
                       ranges = IRanges(start = df.topmed.sub$start,
                                        end = df.topmed.sub$end,
                                        names = df.topmed.sub$SNP),
                       neglog10p = -log10(df.topmed.sub$PVAL))



    # Make Tracks

    col.1kg <- "black"
    if (mir %in% emiRs.1kg) {
        col.1kg <- "darkgreen"
    } else {
        col.1kg <- "darkred"
    }

    if (length(gr.1kg) > 0) {
    dTrack.1kg <- DataTrack(gr.1kg,
                             name = "1KG",
                             type = "p",
                             baseline=-log10(p.value.nominal.threshold.1kg),
                             ylim=c(0,max(gr.1kg$neglog10p)),
                             col.baseline="black",
                             lty.baseline=2,
                             lwd.baseline=1,
                             legend=FALSE,
                             col=col.1kg)
    } else {
    dTrack.1kg <- DataTrack()
}
    col.topmed <- "black"
    if (mir %in% emiRs.topmed) {
        col.topmed <- "darkgreen"
    } else {
        col.topmed <- "darkred"
    }

    dTrack.topmed <- DataTrack(gr.topmed,
                            name = "TOPMed",
                            type = "p",
                            baseline=-log10(p.value.nominal.threshold.topmed),
                            ylim=c(0,max(gr.topmed$neglog10p)),
                            col.baseline="black",
                            lty.baseline=2,
                            lwd.baseline=1,
                            legend=FALSE,
                            col=col.topmed)

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    # # Gene region track
    # cat("Adding gene region track\n")
    #
    # grtrack <- BiomartGeneRegionTrack(genome = "hg19",
    #                                   biomart = genemart,
    #                                   start = startRange,
    #                                   end = endRange,
    #                                   chromosome = chr,
    #                                   showId = TRUE,
    #                                   geneSymbols = TRUE,
    #                                   rotation.title = 90,
    #                                   transcriptAnnotation = "symbol",
    #                                   name = "ENSEMBL",
    #                                   collapseTranscripts = TRUE,
    #                                   just.group="below")
    #
    # mirTrack <- AnnotationTrack(gr.mirs,
    #                             genome = "hg38",
    #                             chromosome = chr,
    #                             name = "miRNA",
    #                             rotation.title = 90,
    #                             id = gr.mirs$Name,
    #                             showFeatureId = TRUE,
    #                             cex = .6,
    #                             fontcolor.item = "black",
    #                             shape = "box")

    pdf(filename)

    plotTracks(list(itrack, gtrack, dTrack.1kg, dTrack.topmed),
               chromosome = chr,
               from = startRange,
               to = endRange,
               transcriptAnnotation="symbol",
               add53=TRUE,
               showBandID=TRUE,
               cex.bands=0.7,
               stackHeight=0.8,
               background.title = "white",
               col.axis="black",
               col.title="black",
               cex.title=1,
               cex.axis=.8,
               just.group="below",
               main = mir)

    dev.off()
}

