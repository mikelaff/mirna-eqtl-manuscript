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
dir.pdfs <- here("doc/emmax/pdfs/topmed_threshold/")
#dir.pngs <- here("doc/emmax/pngs/")

# INPUT FILES ##########################################################################################################
# GRangesList containing SNP-miR associations within 1MB of each miR
#grl.1kg.rds <- here("results/emmax/association_results/20191101_mirQTL_1kg/compiled/20191101_mirQTL_1kg_GRangesList.rds")
grl.topmed.rds <- here("results/emmax/association_results/20191101_mirQTL_topmed/compiled/20191101_mirQTL_topmed_GRangesList.rds")
# Dataframe containing the association results
#df.1kg.rds <- here("results/emmax/association_results/20191101_mirQTL_1kg/compiled/20191101_mirQTL_1kg_dataFrame.rds")
df.topmed.rds <- here("results/emmax/association_results/20191101_mirQTL_topmed/compiled/20191101_mirQTL_topmed_dataFrame.rds")

# genes (mirs) used in this analysis
#genes.1kg.hg19.txt <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/20191101_mirQTL_genes_hg19.txt")
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

# # 1kg GRangesList
# grl.1kg <- readRDS(grl.1kg.rds)
#
# #seqlevelsStyle(grl.1kg)
# seqlevels(grl.1kg) <- CHROMS
# seqinfo(grl.1kg) <- Seqinfo(genome = "hg19")
#
# # 1kg GRanges
# gr.1kg <- unlist(grl.1kg)
#
# # chain for hg19 to hg38 conversion
# path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
# ch <- import.chain(path)
#
# # GRangesList of GRanges conversion
# lo <- liftOver(gr.1kg, ch)
#
# # keep only GRanges that successfully lifted over, convert to GRanges
# gr.1kg.hg38 <- unlist(lo[lengths(lo) == 1])
# names(gr.1kg.hg38) <- paste(gr.1kg.hg38$UniName, gr.1kg.hg38$SNP, sep = ".")
#
# rm(lo, gr.1kg, grl.1kg, ch)
#
# df.1kg.hg38 <- as_tibble(gr.1kg.hg38)
# rm(gr.1kg.hg38)

# import topmed results
grl.topmed.hg38 <- readRDS(grl.topmed.rds)
gr.topmed.hg38 <- unlist(grl.topmed.hg38)

df.topmed.hg38 <- as_tibble(gr.topmed.hg38)
rm(grl.topmed.hg38, gr.topmed.hg38)

# Format ? #########

df.topmed.hg38 %>%
    filter(A1.HOM.count != 1, HET.count > 1) -> df.topmed.hg38.filt.strict

df.topmed.hg38 %>%
    filter(A1.HOM.count > 1 | HET.count > 1) -> df.topmed.hg38.filt.loose

df.topmed.hg38.filt.strict$PVAL.adj <- p.adjust(df.topmed.hg38.filt.strict$PVAL, method = "fdr")
df.topmed.hg38.filt.loose$PVAL.adj <- p.adjust(df.topmed.hg38.filt.loose$PVAL, method = "fdr")

# summarize mirs by minimum adjusted p-value, threshold to find emiRs
df.topmed.hg38.filt.strict %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.strict

df.topmed.hg38.filt.loose %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(PVAL.adj), min.pval = min(PVAL)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.loose

df.topmed.hg38.filt.strict %>%
    filter(PVAL.adj < FDR.THRESHOLD) %>%
    top_n(1, PVAL) %>%
    pull(PVAL) -> p.value.nominal.threshold.strict

df.topmed.hg38.filt.loose %>%
    filter(PVAL.adj < FDR.THRESHOLD) %>%
    top_n(1, PVAL) %>%
    pull(PVAL) -> p.value.nominal.threshold.loose



# Plot ? ##########

emiRs.strict.only <- emiRs.strict[!emiRs.strict %in% emiRs.loose]
emiRs.loose.only <- emiRs.loose[!emiRs.loose %in% emiRs.strict]

emiRs.combined <- unique(c(emiRs.strict, emiRs.loose))

df.genes <- read_tsv(genes.topmed.hg38.txt, col_names = c("uniqename", "chr", "start", "end"))
df.emirs <- dplyr::filter(df.genes, uniqename %in% emiRs.combined)

# loop over each emir

genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

for (i in 1:length(df.emirs$uniqename)) {
    #for (i in 294:306) {
    cat(df.emirs$uniqename[i], "\n")

    mir <- df.emirs$uniqename[i]

    if (mir %in% emiRs.strict.only) {
        filename <- file.path(dir.pdfs, paste0("/strict_only/", mir, ".pdf"))
    } else if (mir %in% emiRs.loose.only) {
        filename <- file.path(dir.pdfs, paste0("/loose_only/", mir, ".pdf"))
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

    # Strict Data GRanges ##############################

    df.topmed.hg38.filt.strict %>%
        filter(UniName == !!mir) -> df.strict.sub

    gr.strict <- GRanges(seqnames = df.strict.sub$seqnames,
                         ranges = IRanges(start = df.strict.sub$start,
                                          end = df.strict.sub$end,
                                          names = df.strict.sub$SNP),
                         neglog10p = -log10(df.strict.sub$PVAL))

    # Loose Data GRanges #####################################

    df.topmed.hg38.filt.loose %>%
        filter(UniName == !!mir) -> df.loose.sub

    gr.loose <- GRanges(seqnames = df.loose.sub$seqnames,
                        ranges = IRanges(start = df.loose.sub$start,
                                         end = df.loose.sub$end,
                                         names = df.loose.sub$SNP),
                        neglog10p = -log10(df.loose.sub$PVAL))



    # Make Tracks

    col.strict <- "black"
    if (mir %in% emiRs.strict) {
        col.strict <- "darkgreen"
    } else {
        col.strict <- "darkred"
    }

    if (length(gr.strict) > 0) {
        dTrack.strict <- DataTrack(gr.strict,
                                   name = "Strict",
                                   type = "p",
                                   baseline=-log10(p.value.nominal.threshold.strict),
                                   ylim=c(0,max(gr.strict$neglog10p)),
                                   col.baseline="black",
                                   lty.baseline=2,
                                   lwd.baseline=1,
                                   legend=FALSE,
                                   col=col.strict)
    } else {
        dTrack.strict <- DataTrack()
    }

    col.loose <- "black"
    if (mir %in% emiRs.loose) {
        col.loose <- "darkgreen"
    } else {
        col.loose <- "darkred"
    }

    dTrack.loose <- DataTrack(gr.loose,
                              name = "Loose",
                              type = "p",
                              baseline=-log10(p.value.nominal.threshold.loose),
                              ylim=c(0,max(gr.loose$neglog10p)),
                              col.baseline="black",
                              lty.baseline=2,
                              lwd.baseline=1,
                              legend=FALSE,
                              col=col.loose)

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

    plotTracks(list(itrack, gtrack, dTrack.strict, dTrack.loose),
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

