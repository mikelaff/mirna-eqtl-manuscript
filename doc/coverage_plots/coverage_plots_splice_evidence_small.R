# plot small rna seq coverage over mirna and genome annotations

library(here)
library(dplyr)
library(readxl)
library(readr)
library(magrittr)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(rtracklayer)

# OUTPUT ###############################################################################################################


# INPUT ################################################################################################################
# directory with small rna-seq bam and bai index files
dir.small.bam <- here("results/small_rna_seq_bams/")
dir.small.bam <- "~/scr/small_rna_seq_bams/"

# directory of total rna-seq bam and bai index files
dir.total.bam <- "/proj/steinlab/projects/FetalTissueQTL/merged_bams/cw/"

# genotypes directory
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# gene expression rse
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")
# mirna expression rse
mirna.rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# vsd for the mirQTL
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
# snp for genotypes
esnp <- "chr10:103394332:A:G"
# emir
emir <- "hsa-mir-1307_hsa-miR-1307-5p"

intron1 <- "10:103394394:103396409:clu_12029_+"
intron2 <- "10:103392466:103394199:clu_12029_+"

chr <- "chr10"
gen <- "hg38"

padding <- 1000
start <- 103392466 - padding
end <- 103396409 + padding

# plotting sequences for data tracks
startseq <- seq(start, end)
endseq <- seq(start, end)

CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Load RSE #############################################################################################################

# mirna expression
vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd.mirna)]

# # set seqnames for gene expression to UCSC
# seqlevelsStyle(rse.gene) <- "UCSC"
# genome(rse.gene) <- "hg38"
# # remove non-standard chroms
# seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS

# threshold for expression
rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

rm(rse.gene, dds.gene)


ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]


# get genotypes at this snp
df.genotypes <- read_table2(paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".raw"))
df.genotypes %<>%
  dplyr::select(donor_id = 1,
                dnaid = 2,
                genotypes_number_of_G = 7)


samples <- as_tibble(colData(vsd.mirna))
samples %<>%
  mutate(donor_id = paste0("D", donor_id))

samples %<>%
  left_join(df.genotypes, by = "donor_id")

# RNAID of AA genotypes
rnaid.aa <- samples$rnaid[samples$genotypes_number_of_G == 0]
rnaid.ag <- samples$rnaid[samples$genotypes_number_of_G == 1]
rnaid.gg <- samples$rnaid[samples$genotypes_number_of_G == 2]


# build gg track
# Import Bams
chr <- "chr10"
samples %>%
  dplyr::filter(rnaid %in% rnaid.gg) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$rnaid

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bam"),
                         index=paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bai"),
                         param=param)
  # convert to coverage
  bamcoverage <- coverage(bam)
  # index into bamcoverage object
  chr_index <- which(names(bamcoverage) == chr)
  # save coverage values in matrix
  values[k,] <- as.numeric(bamcoverage[chr_index][[1]])[(start):(end)]
}

rm(bam, bamcoverage, chr_index, param, k)

# Scale coverage values #############
# divide the values by the normalization factor
scaledValues <- apply(values, 2, function (x) x / sample.data$sizeFactor)
# mean values (if only one sample, then scaledValues is already a vector)
if(nrow(sample.data) > 1) {
  meanValues <- colMeans(scaledValues[,])
} else {
  meanValues <- scaledValues
}
# total counts
#totalValues <- colSums(values[,])


gg.meanValues <- meanValues

# build ag track
# Import Bams
chr <- "chr10"
samples %>%
  dplyr::filter(rnaid %in% rnaid.ag) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$rnaid

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bam"),
                         index=paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bai"),
                         param=param)
  # convert to coverage
  bamcoverage <- coverage(bam)
  # index into bamcoverage object
  chr_index <- which(names(bamcoverage) == chr)
  # save coverage values in matrix
  values[k,] <- as.numeric(bamcoverage[chr_index][[1]])[(start):(end)]
}

rm(bam, bamcoverage, chr_index, param, k)

# Scale coverage values #############
# divide the values by the normalization factor
scaledValues <- apply(values, 2, function (x) x / sample.data$sizeFactor)
# mean values (if only one sample, then scaledValues is already a vector)
if(nrow(sample.data) > 1) {
  meanValues <- colMeans(scaledValues[,])
} else {
  meanValues <- scaledValues
}
# total counts
#totalValues <- colSums(values[,])


ag.meanValues <- meanValues

# build aa track
# Import Bams
chr <- "chr10"
samples %>%
  dplyr::filter(rnaid %in% rnaid.aa) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$rnaid

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bam"),
                         index=paste0(dir.small.bam, sample.data$rnaid[k], ".smallRNAseq.duplicateMarked.sortedByCoord.bai"),
                         param=param)
  # convert to coverage
  bamcoverage <- coverage(bam)
  # index into bamcoverage object
  chr_index <- which(names(bamcoverage) == chr)
  # save coverage values in matrix
  values[k,] <- as.numeric(bamcoverage[chr_index][[1]])[(start):(end)]
}

rm(bam, bamcoverage, chr_index, param, k)

# Scale coverage values #############
# divide the values by the normalization factor
scaledValues <- apply(values, 2, function (x) x / sample.data$sizeFactor)
# mean values (if only one sample, then scaledValues is already a vector)
if(nrow(sample.data) > 1) {
  meanValues <- colMeans(scaledValues[,])
} else {
  meanValues <- scaledValues
}
# total counts
#totalValues <- colSums(values[,])


aa.meanValues <- meanValues


max.val <- max(aa.meanValues, ag.meanValues, gg.meanValues)

# AA coverage track
cat("Adding total-RNA AA coverage track\n")

color <- "red"

fillcolor <- c(color,color)
linecolor <- color
total.aa <- DataTrack(start = startseq,
                      end = endseq,
                      data = ag.meanValues,
                      chromosome = chr,
                      genome = gen,
                      type = "polygon",
                      ylim = c(0, max.val),
                      name = paste("AA\nMean Counts\n(n=", length(rnaid.aa), ")", sep=""),
                      fill.mountain = fillcolor,
                      col = linecolor,
                      baseline = 0,
                      col.baseline = "black")


# AG coverage track
cat("Adding total-RNA AG coverage track\n")

color <- "purple"

fillcolor <- c(color,color)
linecolor <- color
total.ag <- DataTrack(start = startseq,
                      end = endseq,
                      data = ag.meanValues,
                      chromosome = chr,
                      genome = gen,
                      type = "polygon",
                      ylim = c(0, max.val),
                      name = paste("AG\nMean Counts\n(n=", length(rnaid.ag), ")", sep=""),
                      fill.mountain = fillcolor,
                      col = linecolor,
                      baseline = 0,
                      col.baseline = "black")

# GG coverage track
cat("Adding total-RNA GG coverage track\n")

color <- "blue"

fillcolor <- c(color,color)
linecolor <- color
total.gg <- DataTrack(start = startseq,
                      end = endseq,
                      data = gg.meanValues,
                      chromosome = chr,
                      genome = gen,
                      type = "polygon",
                      ylim = c(0, max.val),
                      name = paste("GG\nMean Counts\n(n=", length(rnaid.gg), ")", sep=""),
                      fill.mountain = fillcolor,
                      col = linecolor,
                      baseline = 0,
                      col.baseline = "black")







# Create tracks ##################


# Ideogram track
cat("Adding ideogram track\n")

itrack <- IdeogramTrack(genome = gen, chromosome = chr, bands = ideo)

# Genome axis track
cat("Adding genome axis track\n")

gtrack <- GenomeAxisTrack(exponent = 0)

# Gene region track
cat("Adding gene region track\n")

genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

grtrack <- BiomartGeneRegionTrack(genome = gen,
                                  biomart = genemart,
                                  start = start,
                                  end = end,
                                  chromosome = chr,
                                  showId = TRUE,
                                  geneSymbols = TRUE,
                                  rotation.title = 1,
                                  transcriptAnnotation = "symbol",
                                  name = "REFSEQ",
                                  collapseTranscripts = FALSE,
                                  filter=list(with_refseq_mrna = TRUE),
                                  just.group="below")

miRTrack <- AnnotationTrack(gr.mirna.primary,
                                 genome = "hg38",
                                 chromosome = chr,
                                 name = "miRNA",
                                 rotation.title = 0,
                                 #id = gr.mirna.primary$Name,
                                 group = gr.mirna.primary$Name,
                                 #groupAnnotation = "group",
                                 stacking = "dense",
                                 featureAnnotation="group",
                                 #showID = TRUE,
                                 showFeatureId = TRUE,
                                 fontcolor.item = "black",
                                 shape = "box",
                                 cex = 0.7,
                                 just.id = "left")


overlayTrack <- OverlayTrack(trackList = list(total.aa, total.ag, total.gg))

# Make Plot #################
cat("Plotting...\n")


pdf("~/small.seq.pdf", useDingbats = FALSE)

plotTracks(trackList = list(itrack, gtrack, grtrack, miRTrack, total.aa, total.ag, total.gg),
           chromosome = chr,
           from = start,
           to = end,
           main = "coverage",
           background.title = "white",
           col.axis="black",
           col.title="black",
           cex.title=0.5,
           cex.axis=0.5)

dev.off()


