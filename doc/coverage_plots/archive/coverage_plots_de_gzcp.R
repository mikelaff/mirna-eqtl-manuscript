# plot small rna seq coverage over mirna and genome annotations

library(here)
library(dplyr)
library(magrittr)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(rtracklayer)

# OUTPUT ################
pdf_dir <- here("doc/coverage_plots/pdfs/de_gzcp/")
dir.create(file.path(pdf_dir))

# INPUT ################
# directory with bam and bai index files
bam_dir <- here("results/small_rna_seq_bams/")
# summarized experiment of mirge count data for phenotype data and size factors
mirge_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# mirbase v22 gff file
mirbase_gff <- here("data/mirbase22/hsa.gff3")

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

#nsamples <- length(sampleData$rnaid)
#normfactors <- sampleData$sizeFactor

gz_samples_index <- which(sampleData$tissue_section == "GZ")
cp_samples_index <- which(sampleData$tissue_section == "CP")

# Load miRNA gene models ############################
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

# Specify mirs to plot ################
mirs <- c("hsa-mir-92b", "hsa-mir-124-1", "hsa-mir-124-2", "hsa-mir-124-3", "hsa-mir-9-1", "hsa-mir-9-2", "hsa-mir-9-3")

# loop through mirs
for (mir in mirs) {
  cat("####################", "Working on", mir,"####################" , "\n")
  
  if (!mir %in% df_mirbase$id) {
    cat(mir, "not found. Skipping.\n")
    next
  }

  # Specify plotting region #############
  mir_index <- which(df_mirbase$id == mir)
  if (length(mir_index) > 1) {
    cat(mir, "identifier not unique. Skipping.\n")
    next
  }
  name <- df_mirbase$id[mir_index]
  chr <- as.character(df_mirbase$chromosome[mir_index])
  gen <- "hg38"
  start <- df_mirbase$start[mir_index]
  end <- df_mirbase$end[mir_index]
  
  upstream <- 25
  downstream <- 25
  
  gene.range <- GRanges(seqnames = chr,
                        ranges = IRanges(start = start - upstream, end = end + downstream))
  
  # Import bams ###################
  
  param <- ScanBamParam(what=c("pos","qwidth"),
                        which=GRanges(seqnames(gene.range),
                                      IRanges(start(gene.range),
                                              end(gene.range))),
                        flag=scanBamFlag(isUnmappedQuery=FALSE))
  
  areawidth <- width(gene.range)
  
  # empty values matrix, one each for gz and cp
  values_gz <- matrix(NA, nrow=length(gz_samples_index), ncol=areawidth)
  values_cp <- matrix(NA, nrow=length(cp_samples_index), ncol=areawidth)
  
  # loop over gz bams
  for (i in 1:length(gz_samples_index)) {
    cat("Loading GZ bam file:", i, "of", length(gz_samples_index), "\n")
    index <- gz_samples_index[i]
    
    # read in bam file
    bam <- readGAlignments(file.path(bam_dir, sampleData$bams[index]),
                           index=file.path(bam_dir, sampleData$bais[index]),
                           param=param)
    # convert to coverage
    bamcoverage <- coverage(bam)
    # index into bamcoverage object
    chr_index <- which(names(bamcoverage) == chr)
    # save coverage values in matrix
    values_gz[i,] <- as.numeric(bamcoverage[chr_index][[1]])[(start(gene.range)):(end(gene.range))]
  }
  
  # loop over cp bams
  for (i in 1:length(cp_samples_index)) {
    cat("Loading CP bam file:", i, "of", length(cp_samples_index), "\n")
    index <- cp_samples_index[i]
    
    # read in bam file
    bam <- readGAlignments(file.path(bam_dir, sampleData$bams[index]),
                           index=file.path(bam_dir, sampleData$bais[index]),
                           param=param)
    # convert to coverage
    bamcoverage <- coverage(bam)
    # index into bamcoverage object
    chr_index <- which(names(bamcoverage) == chr)
    # save coverage values in matrix
    values_cp[i,] <- as.numeric(bamcoverage[chr_index][[1]])[(start(gene.range)):(end(gene.range))]
  }
  
  rm(bam, bamcoverage, chr_index, param)
  
  # Scale coverage values #############
  # divide the values by the normalization factor
  scaledValues_gz <- apply(values_gz, 2, function (x) x / sampleData$sizeFactor[gz_samples_index])
  scaledValues_cp <- apply(values_cp, 2, function (x) x / sampleData$sizeFactor[cp_samples_index])
  # mean values (if only one sample, then scaledValues is already a vector)
  if(length(gz_samples_index) > 1) {
    meanValues_gz <- colMeans(scaledValues_gz[,])
  } else {
    meanValues_gz <- scaledValues_gz
  }
  if(length(cp_samples_index) > 1) {
    meanValues_cp <- colMeans(scaledValues_cp[,])
  } else {
    meanValues_cp <- scaledValues_cp
  }
  
  # max value for plotting both tracks on same scale
  maxValue <- max(meanValues_gz, meanValues_cp)
  
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
                                    collapseTranscripts = FALSE,
                                    just.group="below")
  
  # GZ coverage track
  cat("Adding GZ coverage track\n")
  
  fillcolor <- c("darkorange","darkorange")
  linecolor <- "darkorange"
  gztrack <- DataTrack(start = startseq,
                       end = endseq,
                       data = meanValues_gz,
                       chromosome = chr,
                       genome = gen,
                       type = "polygon",
                       ylim = c(0, maxValue),
                       name = paste("GZ Mean Counts\n(n=", length(gz_samples_index), ")", sep=""),
                       fill.mountain = fillcolor,
                       col = linecolor,
                       baseline = 0,
                       col.baseline = "black")
  
  # CP coverage track
  cat("Adding CP coverage track\n")
  
  fillcolor <- c("navy","navy")
  linecolor <- "navy"
  cptrack <- DataTrack(start = startseq,
                       end = endseq,
                       data = meanValues_cp,
                       chromosome = chr,
                       genome = gen,
                       type = "polygon",
                       ylim = c(0, maxValue),
                       name = paste("CP Mean Counts\n(n=", length(cp_samples_index), ")", sep=""),
                       fill.mountain = fillcolor,
                       col = linecolor,
                       baseline = 0,
                       col.baseline = "black")
  
  # miRNA track
  cat("Adding miRNA track\n")
  
  mirtrack <- AnnotationTrack(df_mirbase,
                              genome = gen,
                              chromosome = chr,
                              name = "miRBase V22",
                              featureAnnotation = "id",
                              showFeatureId = TRUE,
                              cex = .6,
                              fontcolor.item = "black")
  
  # Make Plot #################
  cat("Plotting...\n")
  
  plotname <- name
  
  pdf(paste(pdf_dir, plotname, ".pdf", sep=""), useDingbats = FALSE)
  
  plotTracks(trackList = list(itrack, gtrack, gztrack, cptrack, grtrack, mirtrack),
             chromosome = chr,
             from = start(gene.range),
             to = end(gene.range),
             main = plotname,
             background.title = "white",
             col.axis="black",
             col.title="black",
             cex.title=0.9,
             cex.axis=0.9)
  
  dev.off()
  
}
