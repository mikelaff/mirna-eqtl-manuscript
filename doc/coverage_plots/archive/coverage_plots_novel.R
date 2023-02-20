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

# OUTPUT ################
pdf_dir <- here("doc/coverage_plots/pdfs/novel/")
dir.create(file.path(pdf_dir))

# INPUT ################
# directory with bam and bai index files
bam_dir <- here("results/small_rna_seq_bams/")
# summarized experiment of mirge count data for phenotype data and size factors
mirge_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# mirbase v22 gff file
mirbase_gff <- here("data/mirbase22/hsa.gff3")
# friedlander 2014 paper
friedlander_xls <- here("data/friedlander_novel_2014/13059_2013_3254_MOESM3_ESM.xls")
# mirdeep2 novel
mirdeep2_novel_csv <- here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.csv")
# mirge2.0 novel
mirge_novel_tsv <- here("results/mirge2.0/20190130_mirge2.0_novel_miRNAs.tsv")
# compiled novel mirnas to plot
novel_tsv <- here("results/known_and_novel_mirna/20190214_novel_mirnas.tsv")

# Load SE ###############
mirge_se <- readRDS(mirge_se_rds)

# get size factors
dds <- DESeqDataSet(mirge_se, design = ~1)
dds <- estimateSizeFactors(dds)
sampleData <- as.data.frame(colData(dds))

rm(mirge_se, mirge_se_rds, dds)

# Get sample bam files #################
# get bam and bai files, anything with "bam" at end of string ($)
bams <- list.files(bam_dir, ".bam$")
bais <- list.files(bam_dir, ".bai$")

# join files with samples
files <- data.frame(bams=bams, bais=bais, stringsAsFactors=FALSE)
files$rnaid <- sapply(strsplit(files$bams, "\\."), `[`, 1)
sampleData <- dplyr::left_join(sampleData, files, by = "rnaid")

rm(files)

# filter for samples with bams/bais
sampleData <- sampleData[!is.na(sampleData$bams),]

nsamples <- length(sampleData$rnaid)

# Load miRBase miRNA gene models ############################
df_mirbase <- as.data.frame(readGFF(mirbase_gff))
# create width column
df_mirbase$width <- df_mirbase$end - df_mirbase$start + 1
# feature, exon, transcript columns
df_mirbase$feature <- "miRNA"
df_mirbase$exon <- NA
df_mirbase$transcript <- NA
#df_mirbase <- df_mirbase[!is.na(df_mirbase$Derives_from),]
# format df
df_mirbase %<>%
  dplyr::select(chromosome = seqid, start, end, width, feature, gene = ID, exon, transcript, id = Name, derives_from = Derives_from)

# Load friedlander miRNA gene models ############################
df_friedlander <- read_xls(friedlander_xls)

df_friedlander$width <- df_friedlander$`End position (pos included)` - df_friedlander$`Begin pos (1-indexed)` + 1
df_friedlander %<>%
  dplyr::select(chromosome = Chromosome, start = `Begin pos (1-indexed)`, end = `End position (pos included)`, width, id = Identifier)

# Load mirdeep2 miRNA gene models ############################
df_mirdeep <- read_csv(mirdeep2_novel_csv)

df_mirdeep$seqid <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 1)
df_mirdeep$strand <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 3)
df_mirdeep$pos <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 2)
df_mirdeep$start <- as.integer(sapply(strsplit(df_mirdeep$pos, "[.]"), `[`, 1))
df_mirdeep$end <- as.integer(sapply(strsplit(df_mirdeep$pos, "[.]"), `[`, 3))

df_mirdeep %<>%
  dplyr::select(chromosome = seqid, start, end, strand, id = provisional_id) %>%
  dplyr::mutate(width = end - start + 1) %>%
  dplyr::mutate(id = paste("mirdp_", id, sep=""))

# Load mirge2.0 miRNA gene models ############################
df_mirge <- read_tsv(mirge_novel_tsv)

df_mirge %<>% 
  dplyr::filter(!duplicated(`Mature miRNA sequence`)) %>%
  dplyr::select(chromosome = Chr, start = `Start Pos`, end = `End Pos`, strand = Strand, id = provisional_id) %>%
  dplyr::mutate(width = end - start + 1) %>%
  dplyr::mutate(id = paste("mirge_", id, sep=""))

# Specify mirs to plot ################
df_mirs <- read_tsv(novel_tsv)

nmirs <- length(df_mirs$id)

# loop through mirs
for (mir_index in 1:nmirs) {
  cat("####################", "Working on mir", mir_index, "of", nmirs, "####################" , "\n")

  # Specify plotting region #############
  name <- df_mirs$id[mir_index]
  chr <- as.character(df_mirs$chromosome[mir_index])
  gen <- "hg38"
  start <- df_mirs$start[mir_index]
  end <- df_mirs$end[mir_index]
  source <- df_mirs$source[mir_index]
  
  if(source == "mirge") {
    upstream <- 75
    downstream <- 75
    color <- "red"
  } else {
    upstream <- 25
    downstream <- 25
    color <- "blue"
  }
  
  
  gene.range <- GRanges(seqnames = chr,
                        ranges = IRanges(start = start - upstream, end = end + downstream))
  
  # Import bams ###################
  
  param <- ScanBamParam(what=c("pos","qwidth"),
                        which=GRanges(seqnames(gene.range),
                                      IRanges(start(gene.range),
                                              end(gene.range))),
                        flag=scanBamFlag(isUnmappedQuery=FALSE))
  
  areawidth <- width(gene.range)
  
  # empty values matrix
  values <- matrix(NA, nrow=nsamples, ncol=areawidth)

  # loop over bams
  for (i in 1:nsamples) {
    cat("Loading bam file:", i, "of", nsamples, "\n")
    
    # read in bam file
    bam <- readGAlignments(file.path(bam_dir, sampleData$bams[i]),
                           index=file.path(bam_dir, sampleData$bais[i]),
                           param=param)
    # convert to coverage
    bamcoverage <- coverage(bam)
    # index into bamcoverage object
    chr_index <- which(names(bamcoverage) == chr)
    # save coverage values in matrix
    values[i,] <- as.numeric(bamcoverage[chr_index][[1]])[(start(gene.range)):(end(gene.range))]
  }
  
  rm(bam, bamcoverage, chr_index, param)
  
  # Scale coverage values #############
  # divide the values by the normalization factor
  scaledValues <- apply(values, 2, function (x) x / sampleData$sizeFactor)
  # mean values (if only one sample, then scaledValues is already a vector)
  if(nsamples > 1) {
    meanValues <- colMeans(scaledValues[,])
  } else {
    meanValues <- scaledValues
  }
  # total counts
  totalValues <- colSums(values[,])
  
  # Create tracks ##################
  # plotting sequences for data tracks
  startseq <- seq(start(gene.range), end(gene.range))
  endseq <- seq(start(gene.range), end(gene.range))
  
  # Ideogram track
  cat("Adding ideogram track\n")
  
  itrack <- IdeogramTrack(genome = gen, chromosome = chr)
  
  # Genome axis track
  cat("Adding genome axis track\n")
  
  gtrack <- GenomeAxisTrack(exponent = 0)
  
  # Gene region track
  cat("Adding gene region track\n")
  
  genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")
  
  grtrack <- BiomartGeneRegionTrack(genome = gen,
                                    biomart = genemart,
                                    start = start(gene.range),
                                    end = end(gene.range),
                                    chromosome = chr,
                                    showId = TRUE,
                                    geneSymbols = TRUE,
                                    transcriptAnnotation = "symbol",
                                    name = "ENSEMBL",
                                    collapseTranscripts = TRUE,
                                    just.group="below")
  
  # mean coverage track
  cat("Adding mean coverage track\n")
  
  fillcolor <- c(color,color)
  linecolor <- color
  meantrack <- DataTrack(start = startseq,
                         end = endseq,
                         data = meanValues,
                         chromosome = chr,
                         genome = gen,
                         type = "polygon",
                         name = paste("Mean Counts\n(n=", nsamples, ")", sep=""),
                         fill.mountain = fillcolor,
                         col = linecolor,
                         baseline = 0,
                         col.baseline = "black")
  
  # total coverage track
  cat("Adding total coverage track\n")
  
  fillcolor <- c(color,color)
  linecolor <- color
  totaltrack <- DataTrack(start = startseq,
                          end = endseq,
                          data = totalValues,
                          chromosome = chr,
                          genome = gen,
                          type = "polygon",
                          name = paste("Total Counts\n(n=", nsamples, ")", sep=""),
                          fill.mountain = fillcolor,
                          col = linecolor,
                          baseline = 0,
                          col.baseline = "black")
  
  # mirbase track
  cat("Adding mirbase track\n")
  
  mirbasetrack <- AnnotationTrack(df_mirbase,
                                  genome = gen,
                                  chromosome = chr,
                                  name = "miRBase V22",
                                  featureAnnotation = "id",
                                  showFeatureId = TRUE,
                                  cex = .6,
                                  fontcolor.item = "black")
  
  # friedlander track
  cat("Adding friedlander track\n")
  
  friedtrack <- AnnotationTrack(df_friedlander,
                                genome = gen,
                                chromosome = chr,
                                name = "\n\nFried",
                                featureAnnotation = "id",
                                showFeatureId = TRUE,
                                cex = .6,
                                fontcolor.item = "black")
  
  # mirdeep track
  cat("Adding mirdeep track\n")
  
  mirdeeptrack <- AnnotationTrack(df_mirdeep,
                                  genome = gen,
                                  chromosome = chr,
                                  name = "\n\n\nmirdeep",
                                  featureAnnotation = "id",
                                  showFeatureId = TRUE,
                                  cex = .6,
                                  fontcolor.item = "black")
  
  # mirge track
  cat("Adding mirge track\n")
  
  mirgetrack <- AnnotationTrack(df_mirge,
                                genome = gen,
                                chromosome = chr,
                                name = "\n\n\n\nmirge",
                                featureAnnotation = "id",
                                showFeatureId = TRUE,
                                cex = .6,
                                fontcolor.item = "black")
  
  # Make Plot #################
  cat("Plotting...\n")
  
  plotname <- name
  
  pdf(paste(pdf_dir, plotname, ".pdf", sep=""), useDingbats = FALSE)
  
  plotTracks(trackList = list(itrack, gtrack, meantrack, totaltrack, grtrack, mirbasetrack, friedtrack, mirdeeptrack, mirgetrack),
             chromosome = chr,
             from = start(gene.range),
             to = end(gene.range),
             main = plotname,
             background.title = "white",
             col.axis="black",
             col.title="black",
             cex.title=0.5,
             cex.axis=0.5)
  
  dev.off()
  
}
