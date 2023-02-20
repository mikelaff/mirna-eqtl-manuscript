
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(Gviz)
library(gridExtra)
library(biomaRt)

# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/co-localization/pdfs/")

# output folder for combined plots
dir.output <- paste0(dir.pdf)
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

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




# Co-localizations #####################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_eQTLs_dataFrame.rds")
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

# HAUS4 .ps file
haus4.ps <- here("results/emmax_transTotalRNA/association_results/20200928_transTotalRNA_mir4707_HAUS4_batchVars/raw/chr14.dosage.prefiltered.mir4707.ENSG00000092036.pcs.batchVars.ps")
# HAUS4 .ld file
haus4.ld <- here("results/emmax_transTotalRNA/association_results/20200928_transTotalRNA_mir4707_HAUS4_batchVars/raw/chr14.dosage.prefiltered.mir4707.ENSG00000092036.pcs.batchVars.ld")

# mir4707 conditional on edu attain index
edu.attain.cond.ps <- here("results/emmax/association_results/20201012_mirQTLor_4707_cond/ps/chr14.dosage.prefiltered.mirQTLor.hsa-mir-4707_hsa-miR-4707-3p.ps")

# 1000 Genomes EUR plink files
dir.1kg.eur.plink <- here("data/1000genomes_phase3_hg38/EUR.plink/")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# mirQTL eqtls
df.eqtls <- readRDS(eqtls.dataframe.rds)

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# ONLY looking at mir4707 region
eqtl <- df.eqtls$eQTL[14]
emir <- df.eqtls$emiR[14]
esnp <- df.eqtls$eSNP[14]

chr <- df.eqtls$SNP.CHR[14]
bp <- df.eqtls$SNP.BP.hg38[14]

# all variants mirQTL
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
# only for this emir
df.results %>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+")) %>%
    filter(UniName == emir) -> df.results.emir

rm(df.results)

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

rm(df.ld, ld.file)


# import GWAS data: EDU
gr.gwas.edu <- readRDS(gr.edu.rds)
names(gr.gwas.edu) <- gr.gwas.edu$RSID

# subset GWAS data within 1MB of index SNP
gr.gwas.edu.window <- subsetByOverlaps(gr.gwas.edu, ranges = GRanges(seqnames = chr,
                                                                     ranges = IRanges(start = bp - 1e6,
                                                                                      end = bp + 1e6)))
rm(gr.gwas.edu)

df.gwas.edu <- as_tibble(gr.gwas.edu.window)
rm(gr.gwas.edu.window)

#gwas.var <- getRSIDfrom1kgEUR(chr, bp)
gwas.var <- "rs1043209"


# get LD for gwas index
# 1000Genomes EUR bfile
# bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")
# outputPrefix <- paste0(dir.output, "tmp.", gwas.var)
#
# command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 5000 --ld-window-kb 10000 --ld-window-r2 0.2 --out %s",
#                    bfile, gwas.var, outputPrefix)
#
# # call PLINK
# system(command = command)

df.ld.gwas <- read_table(here("data/EUR.rs1043209.ld"))
df.ld.gwas %<>%
    dplyr::select(RSID = SNP_B, R2)

# remove all plink tmp files
#system(command = paste0("rm ", dir.output, "tmp.*"))

# add ld to gwas table
df.gwas.edu %<>%
    left_join(df.ld.gwas, by = "RSID")

# dot color based on LD
df.gwas.edu$color <- "navy"
df.gwas.edu$color[which(df.gwas.edu$R2 >= 0.8)] <- "red"
df.gwas.edu$color[which(df.gwas.edu$R2 >= 0.6 & df.gwas.edu$R2 < 0.8)] <- "orange"
df.gwas.edu$color[which(df.gwas.edu$R2 >= 0.4 & df.gwas.edu$R2 < 0.6)] <- "green"
df.gwas.edu$color[which(df.gwas.edu$R2 >= 0.2 & df.gwas.edu$R2 < 0.4)] <- "lightblue"
df.gwas.edu$color[which(df.gwas.edu$RSID == gwas.var)] <- "purple"

rm(df.ld.gwas)




# conditional data
df.cond <- read_tsv(edu.attain.cond.ps)

# combine tables
df.results.emir %>%
    dplyr::select(-BETA, -P) %>%
    left_join(df.cond, by = "SNP") -> df.cond


# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# plot window
buffer.bp <- 1e5

# always plot the same window, make sure all vars are in the window
startRange.esnp <- bp - buffer.bp
endRange.esnp <- bp + buffer.bp

y.axis.max <- -log10(min(df.results.emir$P, df.cond$P, df.gwas.edu$Pval))

# granges for mirQTL
gr.mirqtl <- GRanges(seqnames = df.results.emir$CHR,
                     ranges = IRanges(start = df.results.emir$BP.hg38,
                                      end = df.results.emir$BP.hg38,
                                      names = df.results.emir$SNP),
                     neglog10p = -log10(df.results.emir$P))

# purple plot
purpleTrack <- DataTrack(gr.mirqtl[esnp],
                         name = "mir4707 -log10(P)",
                         type = "p",
                         legend=FALSE,
                         ylim=c(0,y.axis.max+2),
                         col="purple",
                         cex=1.25,
                         pch=18)

# red plot
if (sum(df.results.emir$color == "red") == 0) {
    redTrack <- DataTrack()
} else {
    redTrack <- DataTrack(gr.mirqtl[df.results.emir$color == "red"],
                          name = "mir4707 -log10(P)",
                          type = "p",
                          legend=FALSE,
                          ylim=c(0,y.axis.max+2),
                          col="red")
}

# orange plot
if (sum(df.results.emir$color == "orange") == 0) {
    orangeTrack <- DataTrack()
} else {
    orangeTrack <- DataTrack(gr.mirqtl[df.results.emir$color == "orange"],
                             name = "mir4707 -log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,y.axis.max+2),
                             col="orange")
}

# green plot
if (sum(df.results.emir$color == "green") == 0) {
    greenTrack <- DataTrack()
} else {
    greenTrack <- DataTrack(gr.mirqtl[df.results.emir$color == "green"],
                            name = "mir4707 -log10(P)",
                            type = "p",
                            legend=FALSE,
                            ylim=c(0,y.axis.max+2),
                            col="green")
}

# lightblue plot
if (sum(df.results.emir$color == "lightblue") == 0) {
    lightblueTrack <- DataTrack()
} else {
    lightblueTrack <- DataTrack(gr.mirqtl[df.results.emir$color == "lightblue"],
                                name = "mir4707 -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,y.axis.max+2),
                                col="lightblue")
}

# navy plot
navyTrack <- DataTrack(gr.mirqtl[df.results.emir$color == "navy"],
                       name = "mir4707 -log10(P)",
                       type = "p",
                       baseline=-log10(nom.p.val),
                       col.baseline="black",
                       lty.baseline=2,
                       lwd.baseline=1,
                       legend=FALSE,
                       ylim=c(0,y.axis.max+2),
                       col="navy")

overlayTrack.mirqtl <- OverlayTrack(trackList = list(navyTrack,
                                                     lightblueTrack,
                                                     greenTrack,
                                                     orangeTrack,
                                                     redTrack,
                                                     purpleTrack))



# granges for edu
gr.edu <- GRanges(seqnames = df.gwas.edu$seqnames,
                  ranges = IRanges(start = df.gwas.edu$start,
                                   end = df.gwas.edu$start,
                                   names = df.gwas.edu$RSID),
                  neglog10p = -log10(df.gwas.edu$Pval))

# purple plot
purpleTrack <- DataTrack(gr.edu[gwas.var],
                         name = "Edu.Att. -log10(P)",
                         type = "p",
                         legend=FALSE,
                         ylim=c(0,y.axis.max+2),
                         col="purple",
                         cex=1.25,
                         pch=18)

# red plot
if (sum(df.gwas.edu$color == "red") == 0) {
    redTrack <- DataTrack()
} else {
    redTrack <- DataTrack(gr.edu[df.gwas.edu$color == "red"],
                          name = "Edu.Att. -log10(P)",
                          type = "p",
                          legend=FALSE,
                          ylim=c(0,y.axis.max+2),
                          col="red")
}

# orange plot
if (sum(df.gwas.edu$color == "orange") == 0) {
    orangeTrack <- DataTrack()
} else {
    orangeTrack <- DataTrack(gr.edu[df.gwas.edu$color == "orange"],
                             name = "Edu.Att. -log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,y.axis.max+2),
                             col="orange")
}

# green plot
if (sum(df.gwas.edu$color == "green") == 0) {
    greenTrack <- DataTrack()
} else {
    greenTrack <- DataTrack(gr.edu[df.gwas.edu$color == "green"],
                            name = "Edu.Att. -log10(P)",
                            type = "p",
                            legend=FALSE,
                            ylim=c(0,y.axis.max+2),
                            col="green")
}

# lightblue plot
if (sum(df.gwas.edu$color == "lightblue") == 0) {
    lightblueTrack <- DataTrack()
} else {
    lightblueTrack <- DataTrack(gr.edu[df.gwas.edu$color == "lightblue"],
                                name = "Edu.Att. -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,y.axis.max+2),
                                col="lightblue")
}

# navy plot
navyTrack <- DataTrack(gr.edu[df.gwas.edu$color == "navy"],
                       name = "Edu.Att. -log10(P)",
                       type = "p",
                       baseline=-log10(GWAS.PVAL.THRESHOLD),
                       col.baseline="black",
                       lty.baseline=2,
                       lwd.baseline=1,
                       legend=FALSE,
                       ylim=c(0,y.axis.max+2),
                       col="navy")

overlayTrack.edu <- OverlayTrack(trackList = list(navyTrack,
                                                  lightblueTrack,
                                                  greenTrack,
                                                  orangeTrack,
                                                  redTrack,
                                                  purpleTrack))



# granges for conditional
gr.cond <- GRanges(seqnames = df.cond$CHR,
                   ranges = IRanges(start = df.cond$BP.hg38,
                                    end = df.cond$BP.hg38,
                                    names = df.cond$SNP),
                   neglog10p = -log10(df.cond$P))

# purple plot
purpleTrack <- DataTrack(gr.cond[esnp],
                         name = "mir4707 Cond.\n-log10(P)",
                         type = "p",
                         legend=FALSE,
                         ylim=c(0,y.axis.max+2),
                         col="purple",
                         cex=1.25,
                         pch=18)

# red plot
if (sum(df.cond$color == "red") == 0) {
    redTrack <- DataTrack()
} else {
    redTrack <- DataTrack(gr.cond[df.cond$color == "red"],
                          name = "mir4707 Cond.\n-log10(P)",
                          type = "p",
                          legend=FALSE,
                          ylim=c(0,y.axis.max+2),
                          col="red")
}

# orange plot
if (sum(df.cond$color == "orange") == 0) {
    orangeTrack <- DataTrack()
} else {
    orangeTrack <- DataTrack(gr.cond[df.cond$color == "orange"],
                             name = "mir4707 Cond.\n-log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,y.axis.max+2),
                             col="orange")
}

# green plot
if (sum(df.cond$color == "green") == 0) {
    greenTrack <- DataTrack()
} else {
    greenTrack <- DataTrack(gr.cond[df.cond$color == "green"],
                            name = "mir4707 Cond.\n-log10(P)",
                            type = "p",
                            legend=FALSE,
                            ylim=c(0,y.axis.max+2),
                            col="green")
}

# lightblue plot
if (sum(df.cond$color == "lightblue") == 0) {
    lightblueTrack <- DataTrack()
} else {
    lightblueTrack <- DataTrack(gr.cond[df.cond$color == "lightblue"],
                                name = "mir4707 Cond.\n-log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,y.axis.max+2),
                                col="lightblue")
}

# navy plot
navyTrack <- DataTrack(gr.cond[df.cond$color == "navy"],
                       name = "mir4707 Cond.\n-log10(P)",
                       type = "p",
                       baseline=-log10(nom.p.val),
                       col.baseline="black",
                       lty.baseline=2,
                       lwd.baseline=1,
                       legend=FALSE,
                       ylim=c(0,y.axis.max+2),
                       col="navy")

overlayTrack.cond <- OverlayTrack(trackList = list(navyTrack,
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
                                       cex=1)
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

pdf("~/Desktop/mir4707_cond.pdf", width = 8, height = 8, useDingbats = FALSE)

plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.mirqtl, overlayTrack.edu, overlayTrack.cond),
           chromosome = chr,
           from = startRange.esnp,
           to = endRange.esnp,
           transcriptAnnotation="symbol",
           add53=TRUE,
           showBandID=TRUE,
           cex.bands=1,
           stackHeight=1,
           background.title = "white",
           col.axis="black",
           col.title="black",
           cex.title=.7,
           cex.axis=.7,
           add=FALSE,
           just.group = "right",
           cex.main=1,
           cex=0.8)

dev.off()
