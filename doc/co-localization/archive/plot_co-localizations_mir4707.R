# co-localize mirQTL eSNPs with fetal brain mRNA-eQTLs

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
dir.pdf <- here("doc/co-localization/pdfs/")

# output folder for combined plots
dir.output <- paste0(dir.pdf, "mir4707/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GWAS munged datasets
gr.edu.rds <- here("data/gwas_datasets/educational_attainment/educational_attainment.hg38.GRanges.rds")
gr.sa.global.rds <- here("data/gwas_datasets/enigma3/ENIGMA3_Global/enigma3.surface_area.global.hg38.GRanges.rds")

# GWAS colocalizations
df.coloc.edu.rds <- here("results/external_data/educational_attainment/educational_attainment_mirQTL_colocalizations.rds")

# 1000 Genomes EUR plink files
dir.1kg.eur.plink <- here("data/1000genomes_phase3_hg38/EUR.plink/")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

GWAS.PVAL.THRESHOLD <- 5e-8

# FUNCTIONS ############################################################################################################
getRSIDfrom1kgEUR <- function(chr, pos) {

    # check for correct chr format
    if (! chr %in% CHROMS) {
        stop(paste("chromosome", chr, "not in correct format!"))
    }

    # load .bim plink file
    df.bim <- read_tsv(paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.bim"), col_names = c("chr", "rsid", "cm", "pos", "a1", "a2"))

    # return corresponding RSID
    rsid <- df.bim$rsid[match(pos, df.bim$pos)]

    rm(df.bim)

    return(rsid)
}

# Import ###############################################################################################################
# mirQTL eqtls
df.eqtls <- readRDS(eqtls.dataframe.rds)

# all variants mirQTL
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
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

# Edu. Attain. #########################################################################################################
# table of colocalized snps
df.coloc.edu <- readRDS(df.coloc.edu.rds)

# combination of eQTL to GWAS index SNP colocalizations
df.coloc.edu$eQTL_GWAS_INDEX <- paste(df.coloc.edu$eQTL, df.coloc.edu$SNP_A.gwas, sep = "_")

# unique eQTL to GWAS index SNP colocs
df.coloc.edu.unique <- dplyr::filter(df.coloc.edu, !duplicated(eQTL_GWAS_INDEX))

# import GWAS data
gr.gwas.edu <- readRDS(gr.edu.rds)
names(gr.gwas.edu) <- gr.gwas.edu$RSID

for (i in 1:nrow(df.eqtls.mir4707)) {

    esnp <- df.eqtls.mir4707$eSNP[i]
    emir <- df.eqtls.mir4707$emiR[i]
    eqtl <- df.eqtls.mir4707$eQTL[i]

    chr <- df.eqtls.mir4707$SNP.CHR[i]
    pos <- df.eqtls.mir4707$SNP.BP.hg38[i]

    gwas <- "Educational_Attainment"
    gwas.short <- "Edu.Attain."
    gwas.var <- getRSIDfrom1kgEUR(chr, pos)

    # output plot file
    output.file <- paste0(dir.output, gwas, "_", df.eqtls.mir4707$eQTL[i], "_", gwas.var, ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")
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

    # get LD for gwas index
    # 1000Genomes EUR bfile
    bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")
    outputPrefix <- paste0(dir.output, "tmp.", gwas.var)

    command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.2 --out %s",
                       bfile, gwas.var, outputPrefix)

    # call PLINK
    system(command = command)

    df.ld.gwas <- read_table(paste0(outputPrefix, ".ld"))
    df.ld.gwas %<>%
        dplyr::select(RSID = SNP_B, R2)

    # remove all plink tmp files
    system(command = paste0("rm ", dir.output, "tmp.*"))

    # subset GWAS data within 1MB of index SNP
    gr.gwas.window <- subsetByOverlaps(gr.gwas.edu, ranges = GRanges(seqnames = chr,
                                                                           ranges = IRanges(start = pos - 1e6,
                                                                                            end = pos + 1e6)))
    df.gwas.window <- as_tibble(gr.gwas.window)

    # add ld to gwas table
    df.gwas.window %<>%
        left_join(df.ld.gwas, by = "RSID")

    # dot color based on LD
    df.gwas.window$color <- "navy"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.8)] <- "red"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.6 & df.gwas.window$R2 < 0.8)] <- "orange"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.4 & df.gwas.window$R2 < 0.6)] <- "green"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.2 & df.gwas.window$R2 < 0.4)] <- "lightblue"
    df.gwas.window$color[which(df.gwas.window$RSID == gwas.var)] <- "purple"

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # plot window
    buffer.bp <- 1e5

    # always plot the same window, make sure all vars are in the window
    startRange.esnp <- pos - buffer.bp
    endRange.esnp <- pos + buffer.bp

    # if outside min/max coords, adjust
    if (startRange.esnp < min.bp) {
        startRange.esnp <- min.bp
    }
    if (endRange.esnp > max.bp) {
        endRange.esnp <- max.bp
    }

    # subset results dataframe and gwas dataframe
    df.mirQTL.plotting <- dplyr::filter(df.results.emir, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)
    df.gwas.plotting <- dplyr::filter(df.gwas.window, start >= startRange.esnp & end <= endRange.esnp)

    # granges for mirQTL
    gr.mirqtl <- GRanges(seqnames = df.mirQTL.plotting$CHR,
                         ranges = IRanges(start = df.mirQTL.plotting$BP.hg38,
                                          end = df.mirQTL.plotting$BP.hg38,
                                          names = df.mirQTL.plotting$SNP),
                         neglog10p = -log10(df.mirQTL.plotting$P))

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
    if (sum(df.mirQTL.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "red"],
                              name = "miRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.mirqtl$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.mirQTL.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "orange"],
                                 name = "miRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.mirqtl$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.mirQTL.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "green"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.mirqtl$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.mirQTL.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "lightblue"],
                                    name = "miRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.mirqtl$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "navy"],
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

    # granges for gwas
    gr.gwas <- GRanges(seqnames = df.gwas.plotting$seqnames,
                       ranges = IRanges(start = df.gwas.plotting$start,
                                        end = df.gwas.plotting$start,
                                        names = df.gwas.plotting$RSID),
                       neglog10p = -log10(df.gwas.plotting$Pval))

    # purple plot
    purpleTrack <- DataTrack(gr.gwas[gwas.var],
                             name = paste(gwas.short, "-log10(P)"),
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                             col="purple",
                             cex=1.25,
                             pch=18)

    # red plot
    if (sum(df.gwas.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "red"],
                              name = paste(gwas.short, "-log10(P)"),
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                              col="red")
    }

    # orange plot
    if (sum(df.gwas.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "orange"],
                                 name = paste(gwas.short, "-log10(P)"),
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                 col="orange")
    }

    # green plot
    if (sum(df.gwas.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "green"],
                                name = paste(gwas.short, "-log10(P)"),
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                col="green")
    }

    # lightblue plot
    if (sum(df.gwas.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "lightblue"],
                                    name = paste(gwas.short, "-log10(P)"),
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "navy"],
                           name = paste(gwas.short, "-log10(P)"),
                           type = "p",
                           baseline=-log10(GWAS.PVAL.THRESHOLD),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                           col="navy")

    overlayTrack.gwas <- OverlayTrack(trackList = list(navyTrack,
                                                       lightblueTrack,
                                                       greenTrack,
                                                       orangeTrack,
                                                       redTrack,
                                                       purpleTrack))

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

    pdf(output.file, width = 9, height = 5)

    plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.mirqtl, overlayTrack.gwas),
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
               main = paste(eqtl, paste(gwas, gwas.var), sep = "\n"),
               cex.main=1)

    dev.off()
}


# Glabal SA ############################################################################################################
# gwas data
gr.gwas.sa.global <- readRDS(gr.sa.global.rds)
names(gr.gwas.sa.global) <- gr.gwas.sa.global$RSID

df.eqtls.mir4707 <- filter(df.eqtls, emiR == "hsa-mir-4707_hsa-miR-4707-3p")

# loop
for (i in 1:nrow(df.eqtls.mir4707)) {

    esnp <- df.eqtls.mir4707$eSNP[i]
    emir <- df.eqtls.mir4707$emiR[i]
    eqtl <- df.eqtls.mir4707$eQTL[i]

    chr <- df.eqtls.mir4707$SNP.CHR[i]
    pos <- df.eqtls.mir4707$SNP.BP.hg38[i]

    gwas <- "Global_Surface_Area"
    gwas.short <- "Global.SA"
    gwas.var <- getRSIDfrom1kgEUR(chr, pos)

    # output plot file
    output.file <- paste0(dir.output, gwas, "_", df.eqtls.mir4707$eQTL[i], "_", gwas.var, ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")
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

    # get LD for gwas index
    # 1000Genomes EUR bfile
    bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")
    outputPrefix <- paste0(dir.output, "tmp.", gwas.var)

    command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.2 --out %s",
                       bfile, gwas.var, outputPrefix)

    # call PLINK
    system(command = command)

    df.ld.gwas <- read_table(paste0(outputPrefix, ".ld"))
    df.ld.gwas %<>%
        dplyr::select(RSID = SNP_B, R2)

    # remove all plink tmp files
    system(command = paste0("rm ", dir.output, "tmp.*"))

    # subset GWAS data within 1MB of index SNP
    gr.gwas.window <- subsetByOverlaps(gr.gwas.sa.global, ranges = GRanges(seqnames = chr,
                                                                     ranges = IRanges(start = pos - 1e6,
                                                                                      end = pos + 1e6)))
    df.gwas.window <- as_tibble(gr.gwas.window)

    # add ld to gwas table
    df.gwas.window %<>%
        left_join(df.ld.gwas, by = "RSID")

    # dot color based on LD
    df.gwas.window$color <- "navy"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.8)] <- "red"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.6 & df.gwas.window$R2 < 0.8)] <- "orange"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.4 & df.gwas.window$R2 < 0.6)] <- "green"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.2 & df.gwas.window$R2 < 0.4)] <- "lightblue"
    df.gwas.window$color[which(df.gwas.window$RSID == gwas.var)] <- "purple"

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # plot window
    buffer.bp <- 1e5

    # always plot the same window, make sure all vars are in the window
    startRange.esnp <- pos - buffer.bp
    endRange.esnp <- pos + buffer.bp

    # if outside min/max coords, adjust
    if (startRange.esnp < min.bp) {
        startRange.esnp <- min.bp
    }
    if (endRange.esnp > max.bp) {
        endRange.esnp <- max.bp
    }

    # subset results dataframe and gwas dataframe
    df.mirQTL.plotting <- dplyr::filter(df.results.emir, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)
    df.gwas.plotting <- dplyr::filter(df.gwas.window, start >= startRange.esnp & end <= endRange.esnp)

    # granges for mirQTL
    gr.mirqtl <- GRanges(seqnames = df.mirQTL.plotting$CHR,
                       ranges = IRanges(start = df.mirQTL.plotting$BP.hg38,
                                        end = df.mirQTL.plotting$BP.hg38,
                                        names = df.mirQTL.plotting$SNP),
                       neglog10p = -log10(df.mirQTL.plotting$P))

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
    if (sum(df.mirQTL.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "red"],
                              name = "miRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.mirqtl$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.mirQTL.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "orange"],
                                 name = "miRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.mirqtl$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.mirQTL.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "green"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.mirqtl$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.mirQTL.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "lightblue"],
                                    name = "miRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.mirqtl$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "navy"],
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

    # granges for gwas
    gr.gwas <- GRanges(seqnames = df.gwas.plotting$seqnames,
                         ranges = IRanges(start = df.gwas.plotting$start,
                                          end = df.gwas.plotting$start,
                                          names = df.gwas.plotting$RSID),
                         neglog10p = -log10(df.gwas.plotting$Pval))

    # purple plot
    purpleTrack <- DataTrack(gr.gwas[gwas.var],
                             name = paste(gwas.short, "-log10(P)"),
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                             col="purple",
                             cex=1.25,
                             pch=18)

    # red plot
    if (sum(df.gwas.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "red"],
                              name = paste(gwas.short, "-log10(P)"),
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                              col="red")
    }

    # orange plot
    if (sum(df.gwas.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "orange"],
                                 name = paste(gwas.short, "-log10(P)"),
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                 col="orange")
    }

    # green plot
    if (sum(df.gwas.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "green"],
                                name = paste(gwas.short, "-log10(P)"),
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                col="green")
    }

    # lightblue plot
    if (sum(df.gwas.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "lightblue"],
                                    name = paste(gwas.short, "-log10(P)"),
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "navy"],
                           name = paste(gwas.short, "-log10(P)"),
                           type = "p",
                           baseline=-log10(GWAS.PVAL.THRESHOLD),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                           col="navy")

    overlayTrack.gwas <- OverlayTrack(trackList = list(navyTrack,
                                                         lightblueTrack,
                                                         greenTrack,
                                                         orangeTrack,
                                                         redTrack,
                                                         purpleTrack))

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

    pdf(output.file, width = 9, height = 5)

    plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.mirqtl, overlayTrack.gwas),
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
               main = paste(eqtl, paste(gwas, gwas.var), sep = "\n"),
               cex.main=1)

    dev.off()
}


# Trans EQTLs #########################################################################################################
# table of colocalized snps
df.coloc.edu <- readRDS(df.coloc.edu.rds)

# combination of eQTL to GWAS index SNP colocalizations
df.coloc.edu$eQTL_GWAS_INDEX <- paste(df.coloc.edu$eQTL, df.coloc.edu$SNP_A.gwas, sep = "_")

# unique eQTL to GWAS index SNP colocs
df.coloc.edu.unique <- dplyr::filter(df.coloc.edu, !duplicated(eQTL_GWAS_INDEX))

# import GWAS data
gr.gwas.edu <- readRDS(gr.edu.rds)
names(gr.gwas.edu) <- gr.gwas.edu$RSID



    esnp <- df.eqtls$eSNP[44]
    emir <- df.eqtls$emiR[44]
    eqtl <- df.eqtls$eQTL[44]

    chr <- df.eqtls$SNP.CHR[44]
    pos <- df.eqtls$SNP.BP.hg38[44]


    # output plot file
    output.file <- paste0(dir.output, "trans_", eqtl, ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")
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


    trans1.ensg <- "ENSG00000110092"
    trans1.chr <- "chr11"
    df.results.trans1 <- read_table2(here("results/emmax_transTotalRNA/association_results/20200803_transTotalRNA_mir4707/raw/chr11.dosage.prefiltered.mir4707.ENSG00000110092.ps"))
    df.results.trans1 <- full_join(df.results.emir, df.results.trans1, by = "SNP", suffix = c(".emir", ".ensg"))
    df.ld.trans1 <- read_table(here("results/emmax_transTotalRNA/association_results/20200803_transTotalRNA_mir4707/clumped/ENSG00000110092.ld"), col_types = cols())
    df.ld.trans1 %<>%
        dplyr::select(SNP = SNP_B, R2)
    df.results.trans1 <- left_join(df.results.trans1, df.ld.trans1, by = "SNP", suffix = c(".emir", ".ensg"))
    df.results.trans1$color.ensg <- "navy"
    df.results.trans1$color.ensg[which(df.results.trans1$R2.ensg >= 0.8)] <- "red"
    df.results.trans1$color.ensg[which(df.results.trans1$R2.ensg >= 0.6 & df.results.trans1$R2.ensg < 0.8)] <- "orange"
    df.results.trans1$color.ensg[which(df.results.trans1$R2.ensg >= 0.4 & df.results.trans1$R2.ensg < 0.6)] <- "green"
    df.results.trans1$color.ensg[which(df.results.trans1$R2.ensg >= 0.2 & df.results.trans1$R2.ensg < 0.4)] <- "lightblue"
    df.results.trans1$color.ensg[which(df.results.trans1$R2.ensg == esnp)] <- "purple"

    trans2.ensg <- "ENSG00000065320"
    trans2.chr <- "chr17"
    df.results.trans2 <- read_table2(here("results/emmax_transTotalRNA/association_results/20200803_transTotalRNA_mir4707/raw/chr17.dosage.prefiltered.mir4707.ENSG00000065320.ps"))


    trans3.ensg <- "ENSG00000177663"
    trans3.chr <- "chr22"
    df.results.trans3 <- read_table2(here("results/emmax_transTotalRNA/association_results/20200803_transTotalRNA_mir4707/raw/chr22.dosage.prefiltered.mir4707.ENSG00000177663.ps"))


    # get LD for gwas index
    # 1000Genomes EUR bfile
    bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")
    outputPrefix <- paste0(dir.output, "tmp.", gwas.var)

    command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.2 --out %s",
                       bfile, gwas.var, outputPrefix)

    # call PLINK
    system(command = command)

    df.ld.gwas <- read_table(paste0(outputPrefix, ".ld"))
    df.ld.gwas %<>%
        dplyr::select(RSID = SNP_B, R2)

    # remove all plink tmp files
    system(command = paste0("rm ", dir.output, "tmp.*"))

    # subset GWAS data within 1MB of index SNP
    gr.gwas.window <- subsetByOverlaps(gr.gwas.edu, ranges = GRanges(seqnames = chr,
                                                                     ranges = IRanges(start = pos - 1e6,
                                                                                      end = pos + 1e6)))
    df.gwas.window <- as_tibble(gr.gwas.window)

    # add ld to gwas table
    df.gwas.window %<>%
        left_join(df.ld.gwas, by = "RSID")

    # dot color based on LD
    df.gwas.window$color <- "navy"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.8)] <- "red"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.6 & df.gwas.window$R2 < 0.8)] <- "orange"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.4 & df.gwas.window$R2 < 0.6)] <- "green"
    df.gwas.window$color[which(df.gwas.window$R2 >= 0.2 & df.gwas.window$R2 < 0.4)] <- "lightblue"
    df.gwas.window$color[which(df.gwas.window$RSID == gwas.var)] <- "purple"

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # plot window
    buffer.bp <- 1e5

    # always plot the same window, make sure all vars are in the window
    startRange.esnp <- pos - buffer.bp
    endRange.esnp <- pos + buffer.bp

    # if outside min/max coords, adjust
    if (startRange.esnp < min.bp) {
        startRange.esnp <- min.bp
    }
    if (endRange.esnp > max.bp) {
        endRange.esnp <- max.bp
    }

    # subset results dataframe and gwas dataframe
    df.mirQTL.plotting <- dplyr::filter(df.results.emir, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)
    df.gwas.plotting <- dplyr::filter(df.gwas.window, start >= startRange.esnp & end <= endRange.esnp)

    # granges for mirQTL
    gr.mirqtl <- GRanges(seqnames = df.mirQTL.plotting$CHR,
                         ranges = IRanges(start = df.mirQTL.plotting$BP.hg38,
                                          end = df.mirQTL.plotting$BP.hg38,
                                          names = df.mirQTL.plotting$SNP),
                         neglog10p = -log10(df.mirQTL.plotting$P))

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
    if (sum(df.mirQTL.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "red"],
                              name = "miRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.mirqtl$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.mirQTL.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "orange"],
                                 name = "miRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.mirqtl$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.mirQTL.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "green"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.mirqtl$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.mirQTL.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "lightblue"],
                                    name = "miRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.mirqtl$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.mirqtl[df.mirQTL.plotting$color == "navy"],
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

    # granges for gwas
    gr.gwas <- GRanges(seqnames = df.gwas.plotting$seqnames,
                       ranges = IRanges(start = df.gwas.plotting$start,
                                        end = df.gwas.plotting$start,
                                        names = df.gwas.plotting$RSID),
                       neglog10p = -log10(df.gwas.plotting$Pval))

    # purple plot
    purpleTrack <- DataTrack(gr.gwas[gwas.var],
                             name = paste(gwas.short, "-log10(P)"),
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                             col="purple",
                             cex=1.25,
                             pch=18)

    # red plot
    if (sum(df.gwas.plotting$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "red"],
                              name = paste(gwas.short, "-log10(P)"),
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                              col="red")
    }

    # orange plot
    if (sum(df.gwas.plotting$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "orange"],
                                 name = paste(gwas.short, "-log10(P)"),
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                 col="orange")
    }

    # green plot
    if (sum(df.gwas.plotting$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "green"],
                                name = paste(gwas.short, "-log10(P)"),
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                col="green")
    }

    # lightblue plot
    if (sum(df.gwas.plotting$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "lightblue"],
                                    name = paste(gwas.short, "-log10(P)"),
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.gwas[df.gwas.plotting$color == "navy"],
                           name = paste(gwas.short, "-log10(P)"),
                           type = "p",
                           baseline=-log10(GWAS.PVAL.THRESHOLD),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(c(gr.gwas$neglog10p, -log10(GWAS.PVAL.THRESHOLD)))),
                           col="navy")

    overlayTrack.gwas <- OverlayTrack(trackList = list(navyTrack,
                                                       lightblueTrack,
                                                       greenTrack,
                                                       orangeTrack,
                                                       redTrack,
                                                       purpleTrack))

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

    pdf(output.file, width = 9, height = 5)

    plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.mirqtl, overlayTrack.gwas),
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
               main = paste(eqtl, paste(gwas, gwas.var), sep = "\n"),
               cex.main=1)

    dev.off()


