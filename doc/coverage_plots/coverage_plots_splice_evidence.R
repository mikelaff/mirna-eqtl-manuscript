# plot small rna seq coverage over mirna and genome annotations

library(here)
library(dplyr)
library(readxl)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(rtracklayer)
library(mikelaffr)

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

# leafcutter splicing data
rse.splice.rds <- here("results/rdata_files/20211217_rse_leafcutter_splicing.rds")

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

clusterNumber <- "clu_12029_+"

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
  dplyr::select(DonorID = 1,
                DNAID = 2,
                genotypes_number_of_G = 7)


samples <- as_tibble(colData(vsd.gene))
# samples %<>%
#   mutate(donor_id = paste0("D", donor_id))

samples %<>%
  left_join(df.genotypes, by = "DonorID")

# RNAID of AA genotypes
rnaid.aa <- samples$RNAID[samples$genotypes_number_of_G == 0]
rnaid.ag <- samples$RNAID[samples$genotypes_number_of_G == 1]
rnaid.gg <- samples$RNAID[samples$genotypes_number_of_G == 2]


# build gg track
# Import Bams
chr <- "10"
samples %>%
  dplyr::filter(RNAID %in% rnaid.gg) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$RNAID

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bam"),
                         index=paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bai"),
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
chr <- "10"
samples %>%
  dplyr::filter(RNAID %in% rnaid.ag) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$RNAID

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bam"),
                         index=paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bai"),
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
chr <- "10"
samples %>%
  dplyr::filter(RNAID %in% rnaid.aa) -> sample.data

param <- ScanBamParam(what=c("pos","qwidth"),
                      which=GRanges(chr,
                                    IRanges(start, end)),
                      flag=scanBamFlag(isUnmappedQuery=FALSE))

areawidth <- end - start + 1

# empty values matrix
values <- matrix(NA, nrow=nrow(sample.data), ncol=areawidth)
rownames(values) <- sample.data$RNAID

# loop over bams
for (k in 1:nrow(sample.data)) {
  cat("Loading bam file:", k, "of", nrow(sample.data), "\n")

  # read in bam file
  bam <- readGAlignments(paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bam"),
                         index=paste0(dir.total.bam, sample.data$RNAID[k], ".totalRNAseq.duplicateMarked.Aligned.sortedByCoord.bai"),
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


pdf("~/total.seq.pdf", useDingbats = FALSE)

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

# Splicing Numbers #####################################################################################################
# ATP5MK
ensg <- "ENSG00000173915"

# mirna expression data
vsd.mirna <- read_rds(vsd.mirna.rds)

# splicing data
rse.splicing <- read_rds(rse.splice.rds)

# get genotypes at this snp
df.genotypes <- read_table(paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".raw"))
df.genotypes %<>%
  dplyr::select(donor_id = 1,
                DNAID = 2,
                genotypes_number_of_G = 7)


samples <- as_tibble(colData(vsd.mirna))
samples %<>%
  mutate(donor_id = paste0("D", donor_id))

samples %<>%
  left_join(df.genotypes, by = "donor_id")

samples %<>%
    mutate(genotypes_number_of_G = factor(genotypes_number_of_G, levels = c(0,1,2), labels = c("AA", "AG", "GG"), ordered = TRUE))

samples %<>%
    dplyr::rename(genotype = genotypes_number_of_G)

# RNAID of genotypes
rnaid.aa <- samples$rnaid[samples$genotype == "AA"]
rnaid.ag <- samples$rnaid[samples$genotype == "AG"]
rnaid.gg <- samples$rnaid[samples$genotype == "GG"]


# subset for this cluster and samples in the mirna analysis
rse.cluster <- rse.splicing[grepl(clusterNumber, rownames(rse.splicing)), colnames(vsd.mirna)]

#assay(rse.cluster, 2)

sum(rowMeans(assay(rse.cluster, 2)))

splice.labels <- tibble(levels = c("10:103392466:103394199:clu_12029_+",
                        "10:103392466:103395746:clu_12029_+",
                        "10:103392466:103396409:clu_12029_+",
                        "10:103394394:103396409:clu_12029_+",
                        "10:103396032:103396409:clu_12029_+"),
                        labels = c("spliceA", "spliceB", "spliceC", "spliceD", "spliceE"))

# rowMeans of AA genotypes
rowMeans(assay(rse.cluster, 2)[,rnaid.aa])
# rowMeans of AG genotypes
rowMeans(assay(rse.cluster, 2)[,rnaid.ag])
# rowMeans of GG genotypes
rowMeans(assay(rse.cluster, 2)[,rnaid.gg])



df.norm.splicing <- as.data.frame(t(assay(rse.cluster, 2)))
colnames(df.norm.splicing) <- splice.labels$labels[match(colnames(df.norm.splicing), table = splice.labels$levels)]
df.norm.splicing$rnaid <- rownames(df.norm.splicing)

df.norm.splicing <- as_tibble(df.norm.splicing)

df.norm.splicing %<>%
    left_join(samples, by = "rnaid")


# ATP5MK vst expression
gene.expression.vst <- as.data.frame(t(assay(vsd.gene[ensg,])))
colnames(gene.expression.vst) <- "ATP5MK.vst"
gene.expression.vst$rnaid <- rownames(gene.expression.vst)

expression.vst <- as.data.frame(t(assay(vsd.mirna["hsa-mir-1307_hsa-miR-1307-5p",], 1)))
colnames(expression.vst) <- "expression.vst"
expression.vst$rnaid <- rownames(expression.vst)

expression.vst.residual <- as.data.frame(t(assay(vsd.mirna["hsa-mir-1307_hsa-miR-1307-5p",], 2)))
colnames(expression.vst.residual) <- "expression.vst.residual"
expression.vst.residual$rnaid <- rownames(expression.vst.residual)

df.expression <- as_tibble(left_join(left_join(expression.vst, expression.vst.residual, by = "rnaid"), gene.expression.vst, by = "rnaid"))
df.expression %<>%
    left_join(samples, by = "rnaid")


df.expression %>%
    ggplot(aes(x = genotype, y = expression.vst.residual)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/mir-1307-5p_boxplot.pdf", height = 2, width = 2)

df.norm.splicing %>%
    ggplot(aes(x = genotype, y = spliceA)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/spliceA_boxplot.pdf", height = 2, width = 2)

df.norm.splicing %>%
    ggplot(aes(x = genotype, y = spliceD)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/splieD_boxplot.pdf", height = 2, width = 2)

df.expression %>%
    ggplot(aes(x = genotype, y = ATP5MK.vst)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1, size = 0.5) +
    plotTheme("figure")

df.expression %>%
    ggplot(aes(x = expression.vst, y = ATP5MK.vst)) +
    geom_point() +
    geom_smooth(method = "lm") +
    plotTheme("figure")






