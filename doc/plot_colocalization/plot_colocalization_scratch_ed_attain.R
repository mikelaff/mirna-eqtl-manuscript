# gviz colocalization scratch

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
dir.pdf <- here("doc/plot_colocalization/pdfs/ed_attain")

dir.create(file.path(dir.pdf))

# INPUT ###############################################################################################################
# directory with bam and bai index files
dir.bam <- here("results/small_rna_seq_bams/")
# mirna ranged summarized experiment for rowranges
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GRanges for mirbase and novel mirna
mirbase.granges.rds <- here("data/gtf_and_granges/20190430_mirbase_v22_mirna_granges.rds")
fried.granges.rds <- here("data/gtf_and_granges/20190430_friedlander_mirna_granges.rds")
nowa.granges.rds <- here("data/gtf_and_granges/20190430_nowakowski_mirna_granges.rds")
mirdeep.granges.rds <- here("data/gtf_and_granges/20190430_mirdeep_novel_mirna_granges.rds")
mirge.granges.rds <- here("data/gtf_and_granges/20190430_mirge_novel_mirna_granges.rds")

# GRangesList containing SNP-miR associations within 1MB of each miR (only loaded if needed)
grangeslist.rds <- here("results/rdata_files/20190808_mirQTL.GRangesList.rds")
# Dataframe containing the GRangesList data (computation is faster on dataframe than the GRangesList)
dataframe.rds <- here("results/rdata_files/20190808_mirQTL.dataframe.rds")

# mirna ranged summarized experiment for rowranges
#rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# list of emiRs (eGenes)
emiRs.tsv <- here("results/emmax/association_results/20190808_mirQTL/emiRs.tsv")
# nominal p-value threshold
nom.pval.thresh.txt <- here("results/emmax/association_results/20190808_mirQTL/nominal.pval.threshold")

# directory of clumped results
clump.dir <- here("results/emmax/association_results/20190808_mirQTL/clumped/")

# Import SNP Associations #############################################################################################
# if the convertered dataframe is available, load it
# else import grangeslist and convert to dataframe
if (file.exists(dataframe.rds)) {
    df.snp.mir.assoc <- readRDS(dataframe.rds)
} else {
    grl <- readRDS(grangeslist.rds)

    # unlist GRangesList into GRanges
    gr.snps <- unlist(grl)
    # split name of each range into snpid and mir
    gr.snps$snpid <- sapply(strsplit(names(gr.snps), "\\."), `[`, 2)
    gr.snps$mir <- sapply(strsplit(names(gr.snps), "\\."), `[`, 1)
    # remove names from GRanges
    names(gr.snps) <- NULL

    # convert GRanges into dataframe
    df.snp.mir.assoc <- as.data.frame(gr.snps)

    # select and rename columns
    df.snp.mir.assoc %<>%
        select(chr = seqnames,
               coord.hg19 = start,
               snpid,
               mir, A1, A2.effect, A1.HOM.count, HET.count, A2.HOM.count, beta, pval)

    # save dataframe
    saveRDS(df.snp.mir.assoc, dataframe.rds)

    # clean up
    rm(grl, gr.snps)
}

# filter snp-mir associations to not be dependent on only 1 homozygous minor sample
df.snp.mir.assoc %<>% filter(A1.HOM.count != 1, HET.count >= 2)

nom.pval.thresh <- as.numeric(read_lines(nom.pval.thresh.txt))

# Import GWAS Data ####################################################################################################
df.gwas <- read_tsv(here("data/gwas_datasets/educational_attainment/GWAS_EA_excl23andMe.txt.gz"))

# Import miRNA GRanges ################################################################################################
# rowranges of rse used for GRanges of miRNAs
rse <- readRDS(rse.rds)
gr.mirs <- rowRanges(rse)

# subset for miRNAs with SNP associations
gr.mirs <- gr.mirs[unique(df.snp.mir.assoc$mir)]

# 864 expressed miRNAs, 7 had no SNPs within 1 Mb = 857 miRs remaining

rm(rse)

# Import emiRs ########################################################################################################
df.emirs <- read_tsv(emiRs.tsv)

# Import Index SNPs ###################################################################################################

# df.esnps <- data.frame(mir=character(), snpid=character(), stringsAsFactors = FALSE)
#
# # loop over each emir, get index snps
# for (i in 1:length(df.emirs$mir)) {
#
#     mir <- df.emirs$mir[i]
#     cat(mir, "\n")
#
#     # read in index snps
#     chr <- df.emirs$seqnames[i]
#     file <- file.path(clump.dir, paste(chr, "mirQTL", mir, "indexSNPs", sep="."))
#     esnps <- read_lines(file)
#
#     # fill in table
#     tmp <- data.frame(mir = mir, snpid = esnps, stringsAsFactors = FALSE)
#     df.esnps <- bind_rows(df.esnps, tmp)
#
#     rm(mir, chr, file, esnps, tmp)
# }
#
# df.esnps$key <- paste0(df.esnps$mir, df.esnps$snpid)
#
# df.snp.mir.assoc$key <- paste0(df.snp.mir.assoc$mir, df.snp.mir.assoc$snpid)
#
# sum(df.esnps$key %in% df.snp.mir.assoc$key)



# Plot ######

# loop over each emir

genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

for (i in 1:length(df.emirs$mir)) {
#for (i in 1:3) {
    cat(df.emirs$mir[i], "\n")

    mir <- df.emirs$mir[i]

    filename <- file.path(dir.pdf, paste0(mir, ".pdf"))

    MB.WINDOW <- 1e6

    startRange <- NULL
    endRange <- NULL

    chr <- df.emirs$seqnames[i]

    # genomic range for SNPs associated with this miRNA
    startRange <- start(gr.mirs[mir]) - MB.WINDOW
    endRange <- end(gr.mirs[mir]) + MB.WINDOW

    # miR-eQTL Data GRanges ##############################

    df.snp.mir.assoc %>%
        filter(mir == !!mir) -> df.snp.mir.sub

    gr.mir <- GRanges(seqnames = df.snp.mir.sub$chr,
                      ranges = IRanges(start = df.snp.mir.sub$coord.hg19,
                                       end = df.snp.mir.sub$coord.hg19,
                                       names = df.snp.mir.sub$snpid),
                      neglog10p = -log10(df.snp.mir.sub$pval))

    # SCZ Data GRanges #####################################

    df.gwas %>%
        filter(POS > startRange, POS < endRange) -> df.gwas.sub

    df.gwas.sub$CHR <- paste0("chr", df.gwas.sub$CHR)

    gr.gwas <- GRanges(seqnames = df.gwas.sub$CHR,
                       ranges = IRanges(start = df.gwas.sub$POS,
                                        end = df.gwas.sub$POS,
                                        names = df.gwas.sub$MarkerName),
                       neglog10p = -log10(df.gwas.sub$Pval))



    # Make Tracks

    dTrack.gwas <- DataTrack(gr.gwas,
                             name = "SCZ",
                             type = "p",
                             baseline=-log10(5e-8),
                             ylim=c(0,max(gr.gwas$neglog10p)),
                             col.baseline="black",
                             lty.baseline=2,
                             lwd.baseline=1,
                             legend=FALSE,
                             col="#192752")

    dTrack.mir <- DataTrack(gr.mir,
                            name = mir,
                            type = "p",
                            baseline=-log10(nom.pval.thresh),
                            ylim=c(0,max(gr.mir$neglog10p+1)),
                            col.baseline="black",
                            lty.baseline=2,
                            lwd.baseline=1,
                            legend=FALSE,
                            col="#192752")

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg19", chromosome = chr)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    # Gene region track
    cat("Adding gene region track\n")

    grtrack <- BiomartGeneRegionTrack(genome = "hg19",
                                      biomart = genemart,
                                      start = startRange,
                                      end = endRange,
                                      chromosome = chr,
                                      showId = TRUE,
                                      geneSymbols = TRUE,
                                      rotation.title = 90,
                                      transcriptAnnotation = "symbol",
                                      name = "ENSEMBL",
                                      collapseTranscripts = TRUE,
                                      just.group="below")

    mirTrack <- AnnotationTrack(gr.mirs,
                                   genome = "hg38",
                                   chromosome = chr,
                                   name = "miRNA",
                                   rotation.title = 90,
                                   id = gr.mirs$Name,
                                   showFeatureId = TRUE,
                                   cex = .6,
                                   fontcolor.item = "black",
                                   shape = "box")

    pdf(filename)

    plotTracks(list(itrack, gtrack, grtrack,dTrack.mir, dTrack.gwas, mirTrack),
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
               just.group="below")

    dev.off()
}


stop()
