# plot mirQTL and blood miRNA-eQTLs together

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(readxl)
library(rtracklayer)
library(liftOver)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.output <- here("doc/mirQTL_plots/pdfs/mirQTL_blood_overlap/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
esnps.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold
nom.p.val.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# Huan et. al. 2015, eQTLs (not corrected for cell type?)
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


# Import #########################################################################################################

df.huan.hg19 <- read_xlsx(huan.eqtls.xlsx, skip = 1)

# convert hg19 to hg38 coordinates

# chain for hg19 to hg38 conversion
path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch <- import.chain(path)

# GRanges of huan eqtls
gr <- makeGRangesFromDataFrame(df.huan.hg19,
                               seqnames.field = "chr.SNP",
                               start.field = "SNP.pos",
                               end.field = "SNP.pos",
                               strand.field = "SNP.strand",
                               ignore.strand = TRUE,
                               keep.extra.columns = TRUE)
# GRangesList of GRanges conversion
lo <- liftOver(gr, ch)

gr.huan.hg38 <- unlist(lo)

rm(lo, gr, ch, path)

# mirQTLs
df.emirs <- readRDS(emirs.dataframe.rds)

# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# ideogram
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna annotations
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# loop over each emir
for (i in 1:nrow(df.emirs)) {

    printMessage(paste("Plotting", i, "of", nrow(df.emirs), ":", df.emirs$emiR[i]))

    emir <- df.emirs$emiR[i]
    chr <- as.character(df.emirs$CHR[i])

    # output plot file
    output.file <- paste0(dir.output, emir, ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # subset gr.huan for these coords
    gr.huan.emir <- gr.huan.hg38[seqnames(gr.huan.hg38) == chr]
    gr.huan.emir <- gr.huan.emir[start(gr.huan.emir) >= min.bp & end(gr.huan.emir) <= max.bp]

    mcols(gr.huan.emir) <- -log10(gr.huan.emir$Pval)

    if(length(gr.huan.emir) == 0) {
        print("no overlap")
        next
    }

    # granges mirQTL
    gr.emir <- GRanges(seqnames = df.results.emir$CHR,
                       ranges = IRanges(start = df.results.emir$BP.hg38,
                                        end = df.results.emir$BP.hg38,
                                        names = df.results.emir$SNP),
                       neglog10p = -log10(df.results.emir$P))

    # data track for this emir
    cat("Adding data track\n")
    dTrack.brain <- DataTrack(gr.emir,
                              name = "Brain",
                              type = "p",
                              baseline=-log10(nom.p.val),
                              ylim=c(0,max(gr.emir$neglog10p)),
                              col.baseline="black",
                              lty.baseline=2,
                              lwd.baseline=1,
                              legend=FALSE,
                              col="black",
                              cex=0.5)

    dTrack.blood <- DataTrack(gr.huan.emir,
                              name = "Blood",
                              type = "p",
                              legend = FALSE,
                              col = "black",
                              cex = 0.5)

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr, bands = ideo)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    pdf(output.file)

    plotTracks(list(itrack, gtrack, dTrack.brain, dTrack.blood),
               chromosome = chr,
               from = min.bp,
               to = max.bp,
               transcriptAnnotation="symbol",
               add53=TRUE,
               showBandID=TRUE,
               cex.bands=0.5,
               stackHeight=0.5,
               background.title = "white",
               col.axis="black",
               col.title="black",
               cex.title=0.7,
               cex.axis=0.7,
               add=TRUE,
               just.group = "right",
               main = emir)

    dev.off()
}

# eSNP overlap?

gr.eqtls <- readRDS(eqtls.granges.rds)

hits <- findOverlaps(gr.eqtls, gr.huan.hg38)

gr.eqtls[queryHits(hits)]
gr.huan.hg38[subjectHits(hits)]



