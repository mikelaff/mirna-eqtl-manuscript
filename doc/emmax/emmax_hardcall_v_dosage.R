# compare 1kg v topmed emmax association test

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/hardcall_v_dosage")
#dir.pngs <- here("doc/emmax/pngs/")

# INPUT FILES ##########################################################################################################
# GRangesList containing SNP-miR associations within 1MB of each miR
grl.hardcall.rds <- here("results/emmax/association_results/20191101_mirQTL_hardcall/compiled/20191101_mirQTL_topmed_hardcall_GRangesList.rds")
grl.dosage.rds <- here("results/emmax/association_results/20191101_mirQTL_dosage/compiled/20191101_mirQTL_topmed_dosage_GRangesList.rds")
# Dataframe containing the association results
df.hardcall.rds <- here("results/emmax/association_results/20191101_mirQTL_hardcall/compiled/20191101_mirQTL_topmed_hardcall_dataFrame.rds")
df.dosage.rds <- here("results/emmax/association_results/20191101_mirQTL_dosage/compiled/20191101_mirQTL_topmed_dosage_dataFrame.rds")

# genes (mirs) used in this analysis
genes.topmed.hg38.txt <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/20191101_mirQTL_genes_hg38.txt")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

FDR.THRESHOLD <- 0.05

# Import SNP Associations ##############################################################################################
df.hardcall <- readRDS(df.hardcall.rds)
df.dosage <- readRDS(df.dosage.rds)

# Import Genes #########################################################################################################
df.genes <- read_tsv(genes.topmed.hg38.txt, col_names = c("uniqename", "chr", "start", "end"))

# Format Tables ########################################################################################################

df.hardcall %<>% as_tibble()
df.hardcall$A1.count <- (df.hardcall$A1.HOM.count * 2) + df.hardcall$HET.count
df.hardcall$A2.count <- (df.hardcall$A2.HOM.count * 2) + df.hardcall$HET.count
df.hardcall$A1.freq <- df.hardcall$A1.count / (df.hardcall$A1.count + df.hardcall$A2.count)

df.dosage %<>% as_tibble()
df.dosage$ALT.freq <- df.dosage$ALT_CTS / df.dosage$OBS_CT
df.dosage$REF.freq <- 1 - df.dosage$ALT.freq
df.dosage$MAF <- ifelse(df.dosage$ALT.freq < df.dosage$REF.freq, df.dosage$ALT.freq, df.dosage$REF.freq)


# Find emiRs ###########################################################################################################

# MAF > 0.01 and (HOM.MINOR != 1 & HET > 1)
# use hardcall data to filter variants in the dosage data
df.hardcall %>%
    filter(A1.freq > 0.01) %>%
    filter(A1.HOM.count != 1 & HET.count > 1) -> df.hardcall.filt

df.dosage %>%
    filter(SNP %in% df.hardcall.filt$SNP) -> df.dosage.filt


# adjusted p-value
df.hardcall.filt$P.adj <- p.adjust(df.hardcall.filt$P, method = "fdr")
df.dosage.filt$P.adj <- p.adjust(df.dosage.filt$P, method = "fdr")

# get emirs
df.hardcall.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(P.adj), min.pval = min(P)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.hardcall

df.dosage.filt %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(P.adj), min.pval = min(P)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.dosage

# get nominal p-value threshold
df.hardcall.filt %>%
    filter(P.adj < FDR.THRESHOLD) %>%
    top_n(1, P) %>%
    pull(P) -> p.value.nominal.hardcall

df.dosage.filt %>%
    filter(P.adj < FDR.THRESHOLD) %>%
    top_n(1, P) %>%
    pull(P) -> p.value.nominal.dosage

# Plot ##########

emiRs.dosage.only <- emiRs.dosage[!emiRs.dosage %in% emiRs.hardcall]
emiRs.hardcall.only <- emiRs.hardcall[!emiRs.hardcall %in% emiRs.dosage]

emiRs.combined <- unique(c(emiRs.dosage, emiRs.hardcall))

df.emirs <- dplyr::filter(df.genes, uniqename %in% emiRs.combined)

# loop over each emir

#genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

for (i in 1:length(df.emirs$uniqename)) {
    #for (i in 294:306) {
    cat(df.emirs$uniqename[i], "\n")

    mir <- df.emirs$uniqename[i]

    if (mir %in% emiRs.hardcall.only) {
        filename <- file.path(dir.pdfs, paste0("/hardcall_only/", mir, ".pdf"))
    } else if (mir %in% emiRs.dosage.only) {
        filename <- file.path(dir.pdfs, paste0("/dosage_only/", mir, ".pdf"))
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

    if (startRange < 0) {
        startRange <- 0
    }

    # Dosage Data ##############################

    df.dosage.filt %>%
        filter(UniName == !!mir) -> df.dosage.sub

    gr.dosage <- GRanges(seqnames = df.dosage.sub$CHR,
                         ranges = IRanges(start = df.dosage.sub$BP.hg38,
                                          end = df.dosage.sub$BP.hg38,
                                          names = df.dosage.sub$SNP),
                         neglog10p = -log10(df.dosage.sub$P))

    # Hardcall Data #####################################

    df.hardcall.filt %>%
        filter(UniName == !!mir) -> df.hardcall.sub

    gr.hardcall <- GRanges(seqnames = df.hardcall.sub$CHR,
                           ranges = IRanges(start = df.hardcall.sub$BP.hg38,
                                            end = df.hardcall.sub$BP.hg38,
                                            names = df.hardcall.sub$SNP),
                           neglog10p = -log10(df.hardcall.sub$P))



    # Make Tracks

    col.dosage <- "black"
    if (mir %in% emiRs.dosage) {
        col.dosage <- "darkgreen"
    } else {
        col.dosage <- "darkred"
    }

    if (length(gr.dosage) > 0) {
        dTrack.dosage <- DataTrack(gr.dosage,
                                 name = "Dosage",
                                 type = "p",
                                 baseline=-log10(p.value.nominal.dosage),
                                 ylim=c(0,max(gr.dosage$neglog10p)),
                                 col.baseline="black",
                                 lty.baseline=2,
                                 lwd.baseline=1,
                                 legend=FALSE,
                                 col=col.dosage)
    } else {
        dTrack.dosage <- DataTrack()
    }

    col.hardcall <- "black"
    if (mir %in% emiRs.hardcall) {
        col.hardcall <- "darkgreen"
    } else {
        col.hardcall <- "darkred"
    }

    dTrack.hardcall <- DataTrack(gr.hardcall,
                            name = "Hardcall",
                            type = "p",
                            baseline=-log10(p.value.nominal.hardcall),
                            ylim=c(0,max(gr.hardcall$neglog10p)),
                            col.baseline="black",
                            lty.baseline=2,
                            lwd.baseline=1,
                            legend=FALSE,
                            col=col.hardcall)

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

    plotTracks(list(itrack, gtrack, dTrack.dosage, dTrack.hardcall),
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

