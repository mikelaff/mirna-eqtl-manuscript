# Dot plots and Geno x Pheno plots for mirQTLor associations

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(DESeq2)
library(mikelaffr)

# OUTPUT ###############################################################################################################
# output folder for combined plots
dir.output <- here("doc/mirQTL_plots/pdfs/mirQTL_dot_emirs_emmax_v_mlma/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis from mlma
eqtls.dataframe.mlma.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor.eQTLs.dataFrame.rds")
eqtls.dataframe.mlma.npc.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor.eQTLs.dataFrame.rds")
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis from emmax
eqtls.dataframe.emmax.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor.eQTLs.dataFrame.rds")
eqtls.dataframe.emmax.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA from mlma
summarized.results.dataframe.mlma.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor_dataFrame.rds")
summarized.results.dataframe.mlma.npc.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA from emmax
summarized.results.dataframe.emmax.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_dataFrame.rds")
summarized.results.dataframe.emmax.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_dataFrame.rds")

# nominal p-value threshold from mlma
nom.p.val.mlma.txt <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor_nomPvalue.txt")
nom.p.val.mlma.npc.txt <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor_nomPvalue.txt")
# nominal p-value threshold from emmax
nom.p.val.emmax.fp.txt <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_nomPvalue.txt")
nom.p.val.emmax.pr.txt <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
# vsd.rds <- here("results/emmax/phenotype_files/20200103_mirQTLor_VST_miRNA_expression_preregress/20200103_mirQTLor_VST_miRNA_expression_preregress_rse.rds")

# Directory for sample genotypes at each index SNP
# dir.genotypes.fp <- here("results/emmax/association_results/20191204_mirQTLor/sample_genotypes/")
# dir.genotypes.pr <- here("results/emmax/association_results/20200103_mirQTLor_preregress/sample_genotypes/")
# dir.genotypes.mlma <- here("results/mlma/association_results/20200114_mirQTLor_mlma/sample_genotypes/")
# dir.genotypes.mlma.npc <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/sample_genotypes/")

# Directory for LD at each index SNP
# dir.ld.emmax.fp <- here("results/emmax/association_results/20191204_mirQTLor/ld/")
# dir.ld.emmax.pr <- here("results/emmax/association_results/20200103_mirQTLor_preregress/ld/")
# dir.ld.mlma <- here("results/mlma/association_results/20200114_mirQTLor_mlma/ld/")
# dir.ld.mlma.npc <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/ld/")

# mirnas used in this analysis
genes.hg38.txt <- here("results/mlma/phenotype_files/20200114_mirQTLor_VST_miRNA_expression/20200114_mirQTLor_genes_hg38.txt")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
# CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
#             "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################

# compiled eqtls
df.eqtls.emmax.fp <- as_tibble(readRDS(eqtls.dataframe.emmax.fp.rds))
df.eqtls.emmax.fp %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

df.eqtls.emmax.pr <- as_tibble(readRDS(eqtls.dataframe.emmax.pr.rds))
df.eqtls.emmax.pr %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

df.eqtls.mlma <- as_tibble(readRDS(eqtls.dataframe.mlma.rds))
df.eqtls.mlma %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

df.eqtls.mlma.npc <- as_tibble(readRDS(eqtls.dataframe.mlma.npc.rds))
df.eqtls.mlma.npc %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

df.eqtls <- bind_rows(df.eqtls.emmax.fp, df.eqtls.emmax.pr, df.eqtls.mlma, df.eqtls.mlma.npc)
df.eqtls %<>%
    filter(!duplicated(df.eqtls$eqtl))
df.eqtls %<>%
    mutate(sig.emmax.fp = eqtl %in% df.eqtls.emmax.fp$eqtl,
           sig.emmax.pr = eqtl %in% df.eqtls.emmax.pr$eqtl,
           sig.mlma = eqtl %in% df.eqtls.mlma$eqtl,
           sig.mlma.npc = eqtl %in% df.eqtls.mlma.npc$eqtl)

rm(df.eqtls.emmax.fp, df.eqtls.emmax.pr, df.eqtls.mlma, df.eqtls.mlma.npc)

# summarized results
df.results.emmax.fp <- as_tibble(readRDS(summarized.results.dataframe.emmax.fp.rds))
df.results.emmax.fp %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.emmax.pr <- as_tibble(readRDS(summarized.results.dataframe.emmax.pr.rds))
df.results.emmax.pr %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.mlma <- as_tibble(readRDS(summarized.results.dataframe.mlma.rds))
df.results.mlma %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.mlma.npc <- as_tibble(readRDS(summarized.results.dataframe.mlma.npc.rds))
df.results.mlma.npc %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# expression data
# vsd <- readRDS(vsd.rds)

# nom p-values
nom.p.val.emmax.fp <- as.numeric(read_lines(nom.p.val.emmax.fp.txt))
nom.p.val.emmax.pr <- as.numeric(read_lines(nom.p.val.emmax.pr.txt))
nom.p.val.mlma <- as.numeric(read_lines(nom.p.val.mlma.txt))
nom.p.val.mlma.npc <- as.numeric(read_lines(nom.p.val.mlma.npc.txt))

# mirs used in this analysis
df.mirs <- read_tsv(genes.hg38.txt, col_names = c("mir", "chr", "start", "end"))

# emirs
df.mirs %>%
    dplyr::filter(mir %in% df.eqtls$emir) -> df.emir

# Plot eQTLs ###########################################################################################################

#genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

trouble.ind <- NULL

# loop over each emir
for (i in 1:nrow(df.emir)) {
    #for (i in unique(trouble.ind)) {

    #for (i in 1:5) {
    printMessage(paste("Plotting", i, "of", nrow(df.emir), ":", df.emir$mir[i]))

    emir <- df.emir$mir[i]
    chr <- df.emir$chr[i]
    start.range <- df.emir$start[i] - 1e6
    stop.range <- df.emir$end[i] + 1e6

    if(start.range < 0) {
        start.range <- 0
    }

    # output plot file
    output.file <- paste0(dir.output, emir, ".pdf")

    # Dot Plots ########################################################################################################

    # subset results for this emiR
    df.results.emmax.fp.emir <- dplyr::filter(df.results.emmax.fp, UniName == emir)
    df.results.emmax.pr.emir <- dplyr::filter(df.results.emmax.pr, UniName == emir)
    df.results.mlma.emir <- dplyr::filter(df.results.mlma, UniName == emir)
    df.results.mlma.npc.emir <- dplyr::filter(df.results.mlma.npc, UniName == emir)

    # granges for this emiR
    gr.emir.emmax.fp <- GRanges(seqnames = df.results.emmax.fp.emir$CHR,
                                ranges = IRanges(start = df.results.emmax.fp.emir$BP.hg38,
                                                 end = df.results.emmax.fp.emir$BP.hg38,
                                                 names = df.results.emmax.fp.emir$SNP),
                                neglog10p = -log10(df.results.emmax.fp.emir$P))
    gr.emir.emmax.pr <- GRanges(seqnames = df.results.emmax.pr.emir$CHR,
                                ranges = IRanges(start = df.results.emmax.pr.emir$BP.hg38,
                                                 end = df.results.emmax.pr.emir$BP.hg38,
                                                 names = df.results.emmax.pr.emir$SNP),
                                neglog10p = -log10(df.results.emmax.pr.emir$P))
    gr.emir.mlma <- GRanges(seqnames = df.results.mlma.emir$CHR,
                            ranges = IRanges(start = df.results.mlma.emir$BP.hg38,
                                             end = df.results.mlma.emir$BP.hg38,
                                             names = df.results.mlma.emir$SNP),
                            neglog10p = -log10(df.results.mlma.emir$P))
    gr.emir.mlma.npc <- GRanges(seqnames = df.results.mlma.npc.emir$CHR,
                                ranges = IRanges(start = df.results.mlma.npc.emir$BP.hg38,
                                                 end = df.results.mlma.npc.emir$BP.hg38,
                                                 names = df.results.mlma.npc.emir$SNP),
                                neglog10p = -log10(df.results.mlma.npc.emir$P))


    # Ideogram, coords, gene annotations
    #startRange.emir <- min(df.results.emmax.fp.emir$BP.hg38)
    #endRange.emir <- max(df.results.emmax.fp.emir$BP.hg38)

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr, bands = ideo)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    # # Gene region track
    # cat("Adding gene region track\n")
    #
    # grtrack.emir <- BiomartGeneRegionTrack(genome = "hg38",
    #                                        biomart = genemart,
    #                                        start = startRange.emir,
    #                                        end = endRange.emir,
    #                                        chromosome = chr,
    #                                        showId = TRUE,
    #                                        geneSymbols = TRUE,
    #                                        rotation.title = 0,
    #                                        transcriptAnnotation = "symbol",
    #                                        name = "REFSEQ",
    #                                        collapseTranscripts = "meta",
    #                                        filter=list(with_refseq_mrna = TRUE),
    #                                        cex=0.1)


    # miRNA track
    cat("Adding miRNA track\n")

    miRTrack <- AnnotationTrack(gr.mirna.primary,
                                genome = "hg38",
                                chromosome = chr,
                                name = "miRNA",
                                rotation.title = 0,
                                #id = gr.mirna.primary$Name,
                                group = gr.mirna.primary$Name,
                                groupAnnotation = "group",
                                #featureAnnotation="id",
                                #showFeatureId = TRUE,
                                fontcolor.item = "black",
                                stacking = "dense",
                                shape = "box")



    # data track emmax.fp
    cat("Adding data track EMMAX\n")
    if (length(gr.emir.emmax.fp) == 0) {
        dTrack.emmax.fp.emir <- DataTrack()
    } else {
        dTrack.emmax.fp.emir <- DataTrack(gr.emir.emmax.fp,
                                          name = "EMMAX\n-log10(P)",
                                          type = "p",
                                          baseline=-log10(nom.p.val.emmax.fp),
                                          ylim=c(0,max(gr.emir.emmax.fp$neglog10p)),
                                          col.baseline="black",
                                          lty.baseline=2,
                                          lwd.baseline=1,
                                          legend=FALSE,
                                          col="black",
                                          cex=0.5)
    }
    # data track emmax.pr
    cat("Adding data track EMMAX.PR\n")
    if (length(gr.emir.emmax.pr) == 0) {
        dTrack.emmax.pr.emir <- DataTrack()
    } else {
    dTrack.emmax.pr.emir <- DataTrack(gr.emir.emmax.pr,
                                      name = "EMMAX.PR\n-log10(P)",
                                      type = "p",
                                      baseline=-log10(nom.p.val.emmax.pr),
                                      ylim=c(0,max(gr.emir.emmax.pr$neglog10p)),
                                      col.baseline="black",
                                      lty.baseline=2,
                                      lwd.baseline=1,
                                      legend=FALSE,
                                      col="black",
                                      cex=0.5)
    }
    # data track mlma
    if (length(gr.emir.mlma) == 0) {
        dTrack.mlma.emir <- DataTrack()
    } else {
    cat("Adding data track MLMA\n")
    dTrack.mlma.emir <- DataTrack(gr.emir.mlma,
                                  name = "MLMA\n-log10(P)",
                                  type = "p",
                                  baseline=-log10(nom.p.val.mlma),
                                  ylim=c(0,max(gr.emir.mlma$neglog10p)),
                                  col.baseline="black",
                                  lty.baseline=2,
                                  lwd.baseline=1,
                                  legend=FALSE,
                                  col="black",
                                  cex=0.5)
    }
    # data track mlma.npc
    if (length(gr.emir.mlma.npc) == 0) {
        dTrack.mlma.npc.emir <- DataTrack()
    } else {
    cat("Adding data track MLMA.NPC\n")
    dTrack.mlma.npc.emir <- DataTrack(gr.emir.mlma.npc,
                                      name = "MLMA.NPC\n-log10(P)",
                                      type = "p",
                                      baseline=-log10(nom.p.val.mlma.npc),
                                      ylim=c(0,max(gr.emir.mlma.npc$neglog10p)),
                                      col.baseline="black",
                                      lty.baseline=2,
                                      lwd.baseline=1,
                                      legend=FALSE,
                                      col="black",
                                      cex=0.5)
    }



    plot.title <- emir

    # new plot
    pdf(output.file, width = 8, height = 6)

    tryCatch(
        {
            plotTracks(list(itrack, gtrack, miRTrack, dTrack.emmax.fp.emir, dTrack.emmax.pr.emir, dTrack.mlma.emir, dTrack.mlma.npc.emir),
                       chromosome = chr,
                       from = start.range,
                       to = stop.range,
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
                       main = plot.title)
        },
        error = function(mess) {
            message(mess)
            trouble.ind <<- c(trouble.ind, i)
        }
    )

    dev.off()
}



