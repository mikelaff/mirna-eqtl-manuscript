# plot haus4 trans


library(here)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/mir4707_trans/pdfs/")

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# HAUS4 trans raw ps file
haus4.ps <- here("results/emmax_transTotalRNA/association_results/20200803_transTotalRNA_mir4707/raw/chr14.dosage.prefiltered.mir4707.ENSG00000092036.ps")

# HAUS4 trans raw ps file expression pcs
haus4.pcs.ps <- here("results/emmax_transTotalRNA/association_results/20200928_transTotalRNA_mir4707_HAUS4/raw/chr14.dosage.prefiltered.mir4707.ENSG00000092036.expressionPCs.ps")

#HAUS4 trans raw ps with pcs and batchVars
haus4.pcs.batchvars.ps <- here("results/emmax_transTotalRNA/association_results/20200928_transTotalRNA_mir4707_HAUS4_batchVars/raw/chr14.dosage.prefiltered.mir4707.ENSG00000092036.pcs.batchVars.ps")
# GLOBALS ##############################################################################################################
FDR.THRESHOLD <- 0.05 / (2 * 14225)

# Import ###############################################################################################################
# mirQTL eqtls
#df.eqtls <- readRDS(eqtls.dataframe.rds)

# all variants mirQTL
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))

# subset for mir4707
df.results %>%
    dplyr::filter(UniName == "hsa-mir-4707_hsa-miR-4707-3p") -> df.results.emir

df.results.emir %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

esnp <- "chr14:22953244:A:G"
pos <- 22953244
chr <- "chr14"

eqtl <- df.results.emir$UniName_SNP[match(esnp, df.results.emir$SNP)]


# max and min base pair which can be plotted
min.bp <- min(df.results.emir$BP.hg38)
max.bp <- max(df.results.emir$BP.hg38)

# plot window
buffer.bp <- 1e5

# always plot the same window, make sure all vars are in the window
startRange.esnp <- pos - buffer.bp
endRange.esnp <- pos + buffer.bp


# ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
ld.file <- paste0(dir.ld, "chr14.hardcall.prefiltered.mirQTLor.chr14:22953244:A:G.ld")
# ld for this esnp
df.ld <- read_table(ld.file, col_types = cols())
df.ld %<>%
    dplyr::select(SNP = SNP_B, R2)

# add ld to results table
df.results.emir %<>%
    left_join(df.ld, by = "SNP")

# dot color based on LD
df.results.emir$color <- "navy"
df.results.emir$color[which(df.results.emir$R2 >= 0.8)] <- "red"
df.results.emir$color[which(df.results.emir$R2 >= 0.6 & df.results.emir$R2 < 0.8)] <- "orange"
df.results.emir$color[which(df.results.emir$R2 >= 0.4 & df.results.emir$R2 < 0.6)] <- "green"
df.results.emir$color[which(df.results.emir$R2 >= 0.2 & df.results.emir$R2 < 0.4)] <- "lightblue"
df.results.emir$color[which(df.results.emir$SNP == esnp)] <- "purple"

# haus4 results
df.haus <- read_tsv(haus4.ps)

# combine tables
df.results.emir %>%
    left_join(df.haus, by = "SNP", suffix = c(".cis", ".trans.haus")) -> df.combo

# haus4 results expression pcs
df.haus.pcs <- read_tsv(haus4.pcs.ps)
colnames(df.haus.pcs) <- c("SNP", "BETA.trans.haus.pcs", "P.trans.haus.pcs")

# combine tables
df.combo %<>%
    left_join(df.haus.pcs, by = "SNP")

# haus4 results expression pcs plus batch
df.haus.pcs.batchvars <- read_tsv(haus4.pcs.batchvars.ps)
colnames(df.haus.pcs.batchvars) <- c("SNP", "BETA.trans.haus.pcs.batch", "P.trans.haus.pcs.batch")

# combine tables
df.combo %<>%
    left_join(df.haus.pcs.batchvars, by = "SNP")

# granges for mirQTL
gr.mirqtl <- GRanges(seqnames = df.combo$CHR,
                     ranges = IRanges(start = df.combo$BP.hg38,
                                      end = df.combo$BP.hg38,
                                      names = df.combo$SNP),
                     neglog10p = -log10(df.combo$P.cis))

# purple plot
purpleTrack <- DataTrack(gr.mirqtl[esnp],
                         name = "miRNA -log10(P)",
                         type = "p",
                         legend=FALSE,
                         ylim=c(0,max(gr.mirqtl$neglog10p)),
                         col="purple",
                         cex=1.25,
                         pch=18)

# red plot
if (sum(df.combo$color == "red") == 0) {
    redTrack <- DataTrack()
} else {
    redTrack <- DataTrack(gr.mirqtl[df.combo$color == "red"],
                          name = "miRNA -log10(P)",
                          type = "p",
                          legend=FALSE,
                          ylim=c(0,max(gr.mirqtl$neglog10p)),
                          col="red")
}

# orange plot
if (sum(df.combo$color == "orange") == 0) {
    orangeTrack <- DataTrack()
} else {
    orangeTrack <- DataTrack(gr.mirqtl[df.combo$color == "orange"],
                             name = "miRNA -log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(gr.mirqtl$neglog10p)),
                             col="orange")
}

# green plot
if (sum(df.combo$color == "green") == 0) {
    greenTrack <- DataTrack()
} else {
    greenTrack <- DataTrack(gr.mirqtl[df.combo$color == "green"],
                            name = "miRNA -log10(P)",
                            type = "p",
                            legend=FALSE,
                            ylim=c(0,max(gr.mirqtl$neglog10p)),
                            col="green")
}

# lightblue plot
if (sum(df.combo$color == "lightblue") == 0) {
    lightblueTrack <- DataTrack()
} else {
    lightblueTrack <- DataTrack(gr.mirqtl[df.combo$color == "lightblue"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.mirqtl$neglog10p)),
                                col="lightblue")
}

# navy plot
navyTrack <- DataTrack(gr.mirqtl[df.combo$color == "navy"],
                       name = "miRNA -log10(P)",
                       type = "p",
                       baseline=-log10(nom.p.val),
                       col.baseline="black",
                       lty.baseline=2,
                       lwd.baseline=1,
                       legend=FALSE,
                       ylim=c(0,max(gr.mirqtl$neglog10p)),
                       col="navy")

overlayTrack.mirqtl <- OverlayTrack(trackList = list(navyTrack,
                                                     lightblueTrack,
                                                     greenTrack,
                                                     orangeTrack,
                                                     redTrack,
                                                     purpleTrack))


# granges for haus
gr.haus <- GRanges(seqnames = df.combo$CHR,
                     ranges = IRanges(start = df.combo$BP.hg38,
                                      end = df.combo$BP.hg38,
                                      names = df.combo$SNP),
                     neglog10p = -log10(df.combo$P.trans.haus))

# black plot
blackTrack <- DataTrack(gr.haus,
                       name = "mRNA -log10(P)",
                       type = "p",
                       baseline=-log10(FDR.THRESHOLD),
                       col.baseline="black",
                       lty.baseline=2,
                       lwd.baseline=1,
                       legend=FALSE,
                       ylim=c(0,max(c(gr.haus$neglog10p, -log10(FDR.THRESHOLD)), na.rm = TRUE)),
                       col="black")

# granges for haus pcs
gr.haus.pcs <- GRanges(seqnames = df.combo$CHR,
                       ranges = IRanges(start = df.combo$BP.hg38,
                                        end = df.combo$BP.hg38,
                                        names = df.combo$SNP),
                       neglog10p = -log10(df.combo$P.trans.haus.pcs))

# black plot
blackTrack.pcs <- DataTrack(gr.haus.pcs,
                            name = "mRNA -log10(P)\nExpressionPC10",
                            type = "p",
                            baseline=-log10(FDR.THRESHOLD),
                            col.baseline="black",
                            lty.baseline=2,
                            lwd.baseline=1,
                            legend=FALSE,
                            ylim=c(0,max(c(gr.haus.pcs$neglog10p, -log10(FDR.THRESHOLD)), na.rm = TRUE)),
                            col="black")

# granges for haus pcs plus batch
gr.haus.pcs.batch <- GRanges(seqnames = df.combo$CHR,
                       ranges = IRanges(start = df.combo$BP.hg38,
                                        end = df.combo$BP.hg38,
                                        names = df.combo$SNP),
                       neglog10p = -log10(df.combo$P.trans.haus.pcs.batch))

# black plot
blackTrack.pcs.batch <- DataTrack(gr.haus.pcs.batch,
                            name = "mRNA -log10(P)\nExpressionPC10+BatchVars",
                            type = "p",
                            baseline=-log10(FDR.THRESHOLD),
                            col.baseline="black",
                            lty.baseline=2,
                            lwd.baseline=1,
                            legend=FALSE,
                            ylim=c(0,max(c(gr.haus.pcs.batch$neglog10p, -log10(FDR.THRESHOLD)), na.rm = TRUE)),
                            col="black")

# Ideogram track
cat("Adding ideogram track\n")

itrack <- IdeogramTrack(genome = "hg38", chromosome = chr, bands = ideo)

# Genome axis track
cat("Adding genome axis track\n")

gtrack <- GenomeAxisTrack(exponent = 0)

# Gene region track
cat("Adding gene region track\n")

grtrack.esnp <- BiomartGeneRegionTrack(genome = "hg38",
                                       biomart = genemart,
                                       start = startRange.esnp,
                                       end = endRange.esnp,
                                       chromosome = chr,
                                       showId = TRUE,
                                       geneSymbols = TRUE,
                                       rotation.title = 0,
                                       transcriptAnnotation = "symbol",
                                       name = "REFSEQ",
                                       collapseTranscripts = "meta",
                                       filter=list(with_refseq_mrna = TRUE),
                                       cex=0.1)
# miRNA track
cat("Adding miRNA track\n")

miRTrack.esnp <- AnnotationTrack(gr.mirna.primary,
                                 genome = "hg38",
                                 chromosome = chr,
                                 name = "miRNA",
                                 rotation.title = 0,
                                 #id = gr.mirna.primary$Name,
                                 group = gr.mirna.primary$Name,
                                 groupAnnotation = "group",
                                 stacking = "squish",
                                 #featureAnnotation="id",
                                 #showFeatureId = TRUE,
                                 fontcolor.item = "black",
                                 shape = "box")

pdf(paste0(dir.pdf, "mir4707_HAUS4_trans_expressionPCs_batchVars.pdf"), width = 12, height = 12)

plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.mirqtl, blackTrack, blackTrack.pcs, blackTrack.pcs.batch),
           chromosome = chr,
           from = startRange.esnp,
           to = endRange.esnp,
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
           add=FALSE,
           just.group = "right",
           main = paste(eqtl, "HAUS4 : ENSG00000092036", sep = "\n"),
           cex.main=1)

dev.off()

