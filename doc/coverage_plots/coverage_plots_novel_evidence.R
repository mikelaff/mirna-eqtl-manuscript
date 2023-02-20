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

# OUTPUT ##############################################################################################################
dir.pdf <- here("doc/coverage_plots/pdfs/novel_2/")
dir.png <- here("doc/coverage_plots/pngs/novel_2/")
dir.create(file.path(dir.png))
dir.create(file.path(dir.pdf))

# INPUT ###############################################################################################################
# directory with bam and bai index files
dir.bam <- here("results/small_rna_seq_bams/")
# mirna ranged summarized experiment
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GRanges for mirbase and novel mirna
mirbase.granges.rds <- here("data/gtf_and_granges/20190430_mirbase_v22_mirna_granges.rds")
fried.granges.rds <- here("data/gtf_and_granges/20190430_friedlander_mirna_granges.rds")
nowa.granges.rds <- here("data/gtf_and_granges/20190430_nowakowski_mirna_granges.rds")
mirdeep.granges.rds <- here("data/gtf_and_granges/20190430_mirdeep_novel_mirna_granges.rds")
mirge.granges.rds <- here("data/gtf_and_granges/20190430_mirge_novel_mirna_granges.rds")

# early and late gest. week values
early.gest <- 14L
late.gest <- 20L
# Load RSE ############################################################################################################
rse.mirna <- readRDS(rse.mirna.rds)

# get size factors
dds <- DESeqDataSet(rse.mirna, design = ~1)
dds <- estimateSizeFactors(dds)
sample.data <- as.data.frame(colData(dds))

rm(rse.mirna)

# Load GRanges Objects ################################################################################################
gr.mirbase <- readRDS(mirbase.granges.rds)
gr.fried <- readRDS(fried.granges.rds)
gr.nowa <- readRDS(nowa.granges.rds)
gr.mirdeep <- readRDS(mirdeep.granges.rds)
gr.mirge <- readRDS(mirge.granges.rds)

gr.nowa <- gr.nowa[!duplicated(gr.nowa$ID)]

# miRs to Plot ########################################################################################################
# threshold dds for novel mirna
dds.novel <- dds[mcols(dds)$source != "miRBase_v22",]
# threshold dds for at least 10 counts in at least 10 samples
dds.novel.thresh <- dds.novel[rowSums(counts(dds.novel) >= 10) >= 10,]

gr.novel <- rowRanges(dds.novel.thresh)

# # find overlaps against known mirnas
# hits <- findOverlaps(gr.novel, gr.mirbase, ignore.strand=TRUE)
#
# # granages, novel mirnas, no overlaps with mirbase mirnas
# gr.novel <- gr.novel[! 1:length(gr.novel) %in% queryHits(hits),]
#
# # # non overlapping ranges
# # hits <- findOverlaps(gr.novel, ignore.strand=TRUE,
# #                      drop.self=TRUE, drop.redundant=TRUE)
# # gr.novel <- gr.novel[! 1:length(gr.novel) %in% subjectHits(hits),]
#
# rm(hits, dds.novel, dds.novel.thresh)
#
# # add in miR-92b and miR-124-1
# gr.novel <- c(gr.mirbase[match("hsa-miR-92b-3p", gr.mirbase$Name)], gr.mirbase[match("hsa-miR-124-3p", gr.mirbase$Name)], gr.novel)

#############
tmp <- read_lines(here("doc/coverage_plots/tmp.txt"))

gr.novel <- gr.novel[which(gr.novel$Derives_from %in% tmp)]

#############

# Get BAM Filenames ###################################################################################################
# get bam and bai files, anything with "bam" at end of string ($)
bams <- list.files(dir.bam, ".bam$")
bais <- list.files(dir.bam, ".bai$")

# join files with samples
files <- data.frame(bams=bams, bais=bais, stringsAsFactors=FALSE)
files$rnaid <- sapply(strsplit(files$bams, "\\."), `[`, 1)
sample.data <- dplyr::left_join(sample.data, files, by = "rnaid")

rm(files)

# filter for samples with bams/bais
sample.data <- sample.data[!is.na(sample.data$bams),]

# filter for quick plotting
sample.data <- sample.data[sample.data$tissue_section == "CP" |
                             sample.data$tissue_section == "GZ" |
                             sample.data$gestation_week == early.gest |
                             sample.data$gestation_week == late.gest, ]

nsamples <- length(sample.data$rnaid)

rm(bams, bais)

# miRs Plot Loop ######################################################################################################
# plot hairpins instead of mature mirnas
hairpins <- gr.novel$Derives_from
hairpins <- hairpins[!duplicated(hairpins)]

# number of hairpins to plot
n <- length(hairpins)
#n <- 5L
# loop through novel hairpins
for (i in 1:n) {
  cat("####################", "Working on hairpin", i, "of", n, "####################" , "\n")

  # Specify Plotting Region #############
  name <- hairpins[i]
  gr.hairpin <- NULL
  if (name %in% gr.mirbase$Derives_from) {
    source <- "miRBase_v22"
    print(paste("Hairpin:", name, "from", source))
    gr.hairpin <- gr.mirbase[which(gr.mirbase$ID == name)]
  } else if (name %in% gr.mirdeep$Derives_from) {
    source <- "miRDeep2"
    print(paste("Hairpin:", name, "from", source))
    gr.hairpin <- gr.mirdeep[which(gr.mirdeep$Name == name)]
  } else if (name %in% gr.mirge$Derives_from) {
    source <- "miRge2.0"
    print(paste("Hairpin:", name, "from", source))
    gr.hairpin <- gr.mirge[which(gr.mirge$Name == name)]
  } else if (name %in% gr.fried$Derives_from) {
    source <- "Friedlander2014"
    print(paste("Hairpin:", name, "from", source))
    gr.hairpin <- gr.fried[which(gr.fried$Name == name)]
  } else if (name %in% gr.nowa$Derives_from) {
    source <- "Nowakowski2018"
    print(paste("Hairpin:", name, "from", source))
    gr.hairpin <- gr.nowa[which(gr.nowa$Name == name)]
  } else {
    stop("Hairpin name not in any granges object")
  }

  if (is.null(gr.hairpin)) {
    stop("No hairpin found")
  }
  if (length(gr.hairpin) != 1) {
    stop("Incorrect number of hairpins")
  }

  # add upstream and downstream basepairs to granges
  padding <- 20

  chr <- as.character(seqnames(gr.hairpin))
  gen <- "hg38"
  start <- start(gr.hairpin + padding)
  end <- end(gr.hairpin + padding)

  # Import Bams ###################

  param <- ScanBamParam(what=c("pos","qwidth"),
                        which=GRanges(chr,
                                      IRanges(start, end)),
                        flag=scanBamFlag(isUnmappedQuery=FALSE))

  areawidth <- end - start + 1

  # empty values matrix
  values <- matrix(NA, nrow=nsamples, ncol=areawidth)
  rownames(values) <- sample.data$rnaid

  # loop over bams
  for (k in 1:nsamples) {
    cat("Loading bam file:", k, "of", nsamples, "\n")

    # read in bam file
    bam <- readGAlignments(file.path(dir.bam, sample.data$bams[k]),
                           index=file.path(dir.bam, sample.data$bais[k]),
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
  if(nsamples > 1) {
    meanValues <- colMeans(scaledValues[,])
  } else {
    meanValues <- scaledValues
  }
  # total counts
  totalValues <- colSums(values[,])

  # gz and cp values
  n.gz <- sum(sample.data$tissue_section == "GZ")
  values.gz <- scaledValues[which(sample.data$tissue_section == "GZ"),]
  # mean values (if only one sample, then values.gz is already a vector)
  if(n.gz > 1) {
    values.gz <- colMeans(values.gz[,])
  }
  n.cp <- sum(sample.data$tissue_section == "CP")
  values.cp <- scaledValues[which(sample.data$tissue_section == "CP"),]
  # mean values (if only one sample, then values.cp is already a vector)
  if(n.cp > 1) {
    values.cp <- colMeans(values.cp[,])
  }
  maxValue.gzcp <- max(values.gz, values.cp)

  # early and late values
  n.early <- sum(sample.data$gestation_week == early.gest & sample.data$tissue_section == "CW")
  values.early <- scaledValues[which(sample.data$gestation_week == early.gest & sample.data$tissue_section == "CW"),]
  # mean values (if only one sample, then values.early is already a vector)
  if(n.early > 1) {
    values.early <- colMeans(values.early[,])
  }
  n.late <- sum(sample.data$gestation_week == late.gest & sample.data$tissue_section == "CW")
  values.late <- scaledValues[which(sample.data$gestation_week == late.gest & sample.data$tissue_section == "CW"),]
  # mean values (if only one sample, then values.late is already a vector)
  if(n.late > 1) {
    values.late <- colMeans(values.late[,])
  }
  maxValue.earlylate <- max(values.early, values.late)

  # Create tracks ##################
  # plotting sequences for data tracks
  startseq <- seq(start, end)
  endseq <- seq(start, end)

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
                                    start = start,
                                    end = end,
                                    chromosome = chr,
                                    showId = TRUE,
                                    geneSymbols = TRUE,
                                    rotation.title = 1,
                                    transcriptAnnotation = "symbol",
                                    name = "ENSEMBL",
                                    collapseTranscripts = TRUE,
                                    just.group="below")

  # GZ coverage track
  cat("Adding GZ coverage track\n")

  color <- "blue"

  fillcolor <- c(color,color)
  linecolor <- color
  track.gz <- DataTrack(start = startseq,
                        end = endseq,
                        data = values.gz,
                        chromosome = chr,
                        genome = gen,
                        type = "polygon",
                        ylim = c(0, maxValue.gzcp),
                        name = paste("GZ\nMean Counts\n(n=", n.gz, ")", sep=""),
                        fill.mountain = fillcolor,
                        col = linecolor,
                        baseline = 0,
                        col.baseline = "black")

  # CP coverage track
  cat("Adding CP coverage track\n")

  color <- "red"

  fillcolor <- c(color,color)
  linecolor <- color
  track.cp <- DataTrack(start = startseq,
                        end = endseq,
                        data = values.cp,
                        chromosome = chr,
                        genome = gen,
                        type = "polygon",
                        ylim = c(0, maxValue.gzcp),
                        name = paste("CP\nMean Counts\n(n=", n.cp, ")", sep=""),
                        fill.mountain = fillcolor,
                        col = linecolor,
                        baseline = 0,
                        col.baseline = "black")

  # early coverage track
  cat("Adding early coverage track\n")

  color <- "blue4"

  fillcolor <- c(color,color)
  linecolor <- color
  track.early <- DataTrack(start = startseq,
                           end = endseq,
                           data = values.early,
                           chromosome = chr,
                           genome = gen,
                           type = "polygon",
                           ylim = c(0, maxValue.earlylate),
                           name = paste("GW ", early.gest, "\nMean Counts\n(n=", n.early, ")", sep=""),
                           fill.mountain = fillcolor,
                           col = linecolor,
                           baseline = 0,
                           col.baseline = "black")

  # late coverage track
  cat("Adding late coverage track\n")

  color <- "red4"

  fillcolor <- c(color,color)
  linecolor <- color
  track.late <- DataTrack(start = startseq,
                          end = endseq,
                          data = values.late,
                          chromosome = chr,
                          genome = gen,
                          type = "polygon",
                          ylim = c(0, maxValue.earlylate),
                          name = paste("GW ", late.gest, "\nMean Counts\n(n=", n.late, ")", sep=""),
                          fill.mountain = fillcolor,
                          col = linecolor,
                          baseline = 0,
                          col.baseline = "black")

  # total coverage track
  cat("Adding mean coverage track\n")

  color <- "black"

  fillcolor <- c(color,color)
  linecolor <- color
  track.mean <- DataTrack(start = startseq,
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

  color <- "black"

  fillcolor <- c(color,color)
  linecolor <- color
  track.total <- DataTrack(start = startseq,
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

  track.mirbase <- AnnotationTrack(gr.mirbase,
                                   genome = gen,
                                   chromosome = chr,
                                   name = "miRBase_v22",
                                   rotation.title = 1,
                                   id = gr.mirbase$Name,
                                   showFeatureId = TRUE,
                                   cex = .6,
                                   fontcolor.item = "black",
                                   shape = "box")

  # friedlander track
  cat("Adding friedlander track\n")

  track.fried <- AnnotationTrack(gr.fried,
                                 genome = gen,
                                 chromosome = chr,
                                 name = "Friedlander2014",
                                 rotation.title = 1,
                                 id = gr.fried$Name,
                                 showFeatureId = TRUE,
                                 cex = .6,
                                 fontcolor.item = "black",
                                 shape = "box")

  # Nowakowski track
  cat("Adding nowakowski track\n")

  track.nowa <- AnnotationTrack(gr.nowa,
                                genome = gen,
                                chromosome = chr,
                                name = "Nowakowski2018",
                                rotation.title = 1,
                                id = gr.nowa$Name,
                                showFeatureId = TRUE,
                                cex = .6,
                                fontcolor.item = "black",
                                shape = "box")

  # mirdeep track
  cat("Adding mirdeep track\n")

  track.mirdeep <- AnnotationTrack(gr.mirdeep,
                                   genome = gen,
                                   chromosome = chr,
                                   name = "miRDeep2",
                                   rotation.title = 1,
                                   id = gr.mirdeep$Name,
                                   showFeatureId = TRUE,
                                   cex = .6,
                                   fontcolor.item = "black",
                                   shape = "box")

  # mirge track
  cat("Adding mirge track\n")

  track.mirge <- AnnotationTrack(gr.mirge,
                                 genome = gen,
                                 chromosome = chr,
                                 name = "miRge2.0",
                                 rotation.title = 1,
                                 id = gr.mirge$Name,
                                 showFeatureId = TRUE,
                                 cex = .6,
                                 fontcolor.item = "black",
                                 shape = "box")

  # Make Plot #################
  cat("Plotting...\n")

  plotname <- gr.hairpin$Name

  pdf(paste(dir.pdf, plotname, ".pdf", sep=""), useDingbats = FALSE)

  plotTracks(trackList = list(itrack, gtrack,
                              track.gz, track.cp,
                              track.early, track.late,
                              track.mean, track.total,
                              grtrack,
                              track.mirbase, track.fried, track.nowa, track.mirdeep, track.mirge),
             chromosome = chr,
             from = start,
             to = end,
             main = plotname,
             background.title = "white",
             col.axis="black",
             col.title="black",
             cex.title=0.5,
             cex.axis=0.5)

  dev.off()

}
