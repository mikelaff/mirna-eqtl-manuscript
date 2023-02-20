# blood miRNA eQTLs that do and dont overlap with miRNA eqtls

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(DESeq2)
library(mikelaffr)


# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/mirQTL_plots/pdfs/")

# output folder for plots
dir.output <- paste0(dir.pdf, "blood_miRNA-eQTL_zoomIN_2/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor conditionally independent association analysis
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA: Primary results
rank1.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA: Secondary results
rank2.results.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/secondary/association_results/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_secondary_variants_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA: Tertiary results
rank3.results.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/tertiary/association_results/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_tertiary_variants_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA: Quarternary results
rank4.results.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/quarternary/association_results/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_quarternary_variants_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA: Quinary results
rank5.results.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/quinary/association_results/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_quinary_variants_dataFrame.rds")

# nominal p-value threshold, eigenMT 5% fdr
high.nom.p.val.txt <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_nomPvalue.txt")
# nominal p-value threshold, eigenMT 10% fdr
med.nom.p.val.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/20200120_mirQTLor_fdr10percent_conditional_nomPvalue.txt")
# nominal p-value threshold, 5% fdr
low.nom.p.val.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_conditional_nomPvalue.txt")


# blood miRNA-eQTLs
blood.eqtls.rds <- here("results/external_data/huan_2015_blood_miRNA-eQTL/huan_2015_blood_miRNA-eQTLs.rds")
# blood miRNA-eQTL data
blood.variants.rds <- here("results/external_data/huan_2015_blood_miRNA-eQTL/huan_2015_blood_miRNA-eQTL_hg38_variants_dataFrame.rds")

# Directory for LD at each mirQTL index SNP
dir.ld <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/ld/")
# LD for blood index snps
dir.blood.ld <- here("results/external_data/huan_2015_blood_miRNA-eQTL/ld/")

# summarized experiment with vst expression values used in this analysis
#vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/sample_genotypes/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Summarized Results ############################################################################################

# compiled eQTLs
df.eqtls <- as_tibble(readRDS(eqtls.dataframe.rds))

df.eqtls %<>%
    mutate(DEGREE = factor(DEGREE, levels = c("primary", "secondary", "tertiary", "quarternary", "quinary"), labels = c(1,2,3,4,5)))

# rank1 variants
df.rank1.results <- as_tibble(readRDS(rank1.results.dataframe.rds))
df.rank1.results %<>%
    filter(UniName %in% df.eqtls$emiR) %>%
    mutate(DEGREE = 1)
# rank2 variants
df.rank2.results <- as_tibble(readRDS(rank2.results.dataframe.rds))
df.rank2.results %<>%
    filter(UniName %in% df.eqtls$emiR) %>%
    mutate(DEGREE = 2)
# rank3 variants
df.rank3.results <- as_tibble(readRDS(rank3.results.dataframe.rds))
df.rank3.results %<>%
    filter(UniName %in% df.eqtls$emiR) %>%
    mutate(DEGREE = 3)
# rank4 variants
df.rank4.results <- as_tibble(readRDS(rank4.results.dataframe.rds))
df.rank4.results %<>%
    filter(UniName %in% df.eqtls$emiR) %>%
    mutate(DEGREE = 4)
# rank5 variants
df.rank5.results <- as_tibble(readRDS(rank5.results.dataframe.rds))
df.rank5.results %<>%
    filter(UniName %in% df.eqtls$emiR) %>%
    mutate(DEGREE = 5)

df.results <- bind_rows(df.rank1.results, df.rank2.results, df.rank3.results, df.rank4.results, df.rank5.results)

df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

rm(df.rank1.results, df.rank2.results, df.rank3.results, df.rank4.results, df.rank5.results)

# nom p value thresholds
high.p.val.thresh <- as.numeric(read_lines(high.nom.p.val.txt))
med.p.val.thresh <- as.numeric(read_lines(med.nom.p.val.txt))
low.p.val.thresh <- as.numeric(read_lines(low.nom.p.val.txt))

#import blood results
df.blood.results <- read_rds(blood.variants.rds)
df.blood.eqtls <- read_rds(blood.eqtls.rds)

df.blood.eqtls %<>%
    mutate(eqtl = paste(miRNA_FHS, snpID, sep = "+"))

# Find Overlaps ########################################################################################################

df.eqtls %>%
    filter(SIGNIFICANCE == "eigenMT_fdr5percent") -> df.eqtls.high

# brain emiRs that are also blood emiRs
df.eqtls.high %>%
    filter(NAME %in% df.blood.eqtls$hsa_miR_name) -> df.candidates

# join with blood eqtls
df.candidates %<>%
    left_join(df.blood.eqtls, by = c("NAME" = "hsa_miR_name"))

df.candidates %<>%
    mutate(combo = paste(eQTL, eqtl, sep = "_"))

df.candidates %<>%
    mutate(plot_name = paste(eQTL, eqtl, sep = "\n"))

# Track Data ###########################################################################################################
# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# Plot Colocalizations #################################################################################################


for (i in 1:nrow(df.candidates)) {
    printMessage(paste("Plotting", i, "of", nrow(df.candidates), ":", df.candidates$eQTL[i]))

    esnp <- df.candidates$eSNP[i]
    emir <- df.candidates$emiR[i]
    eqtl <- df.candidates$eQTL[i]

    chr <- df.candidates$SNP.CHR[i]
    esnp.bp <- df.candidates$SNP.BP.hg38[i]

    # output plot file
    output.file <- paste0(dir.output, df.candidates$combo[i], ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == df.candidates$emiR[i], DEGREE == df.candidates$DEGREE[i])

    # ld file: mirQTL.overlap.index.chr10:104260640:A:C.ld
    ld.file <- paste0(dir.ld, "conditional.mirQTLor.index.", df.candidates$eSNP[i], ".ld")
    # ld for this esnp
    suppressWarnings(
        df.ld <- read_table(ld.file, col_types = cols())
    )
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

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # eSNP window
    buffer.bp <- 3e5

    # always plot the same window, make sure all vars are in the window
    startRange.esnp <- esnp.bp - buffer.bp
    endRange.esnp <- esnp.bp + buffer.bp

    # if outside min/max coords, adjust
    if (startRange.esnp < min.bp) {
        startRange.esnp <- min.bp
    }
    if (endRange.esnp > max.bp) {
        endRange.esnp <- max.bp
    }

    # subset results dataframe
    df.results.esnp <- dplyr::filter(df.results.emir, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)



    # subset for this blood eqtl
    df.blood.results.emir <- dplyr::filter(df.blood.results, miRNA_FHS == df.candidates$miRNA_FHS[i], start >= startRange.esnp & start <= endRange.esnp)

    # ld file: mirQTL.overlap.index.chr10:104260640:A:C.ld
    ld.file <- paste0(dir.blood.ld, "blood.index.", df.candidates$snpID[i], ".ld")

    if (file.exists(ld.file)) {
        # ld for this esnp
        suppressWarnings(
            df.ld <- read_table(ld.file, col_types = cols())
        )
        df.ld %<>%
            dplyr::select(SNP = SNP_B, R2)

        df.blood.results.emir %<>%
            left_join(df.ld, by = c("snpID" = "SNP"))

        # dot color based on LD
        df.blood.results.emir$color <- "navy"
        df.blood.results.emir$color[which(df.blood.results.emir$R2 >= 0.8)] <- "red"
        df.blood.results.emir$color[which(df.blood.results.emir$R2 >= 0.6 & df.blood.results.emir$R2 < 0.8)] <- "orange"
        df.blood.results.emir$color[which(df.blood.results.emir$R2 >= 0.4 & df.blood.results.emir$R2 < 0.6)] <- "green"
        df.blood.results.emir$color[which(df.blood.results.emir$R2 >= 0.2 & df.blood.results.emir$R2 < 0.4)] <- "lightblue"
        df.blood.results.emir$color[which(df.blood.results.emir$snpID == df.candidates$snpID[i])] <- "purple"

    } else {
        df.blood.results.emir$color <- "navy"
    }

    cat("Adding miRNA track\n")


    # granges for this esnp
    gr.esnp <- GRanges(seqnames = df.results.esnp$CHR,
                       ranges = IRanges(start = df.results.esnp$BP.hg38,
                                        end = df.results.esnp$BP.hg38,
                                        names = df.results.esnp$SNP),
                       neglog10p = -log10(df.results.esnp$P))

    # purple plot
    purpleTrack <- DataTrack(gr.esnp[esnp],
                             name = "miRNA -log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(gr.esnp$neglog10p)),
                             col="purple",
                             cex=1.25,
                             pch=18)

    # red plot
    if (sum(df.results.esnp$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.esnp[df.results.esnp$color == "red"],
                              name = "miRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.esnp$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.results.esnp$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.esnp[df.results.esnp$color == "orange"],
                                 name = "miRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.esnp$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.results.esnp$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.esnp[df.results.esnp$color == "green"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.esnp$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.results.esnp$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.esnp[df.results.esnp$color == "lightblue"],
                                    name = "miRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.esnp$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.esnp[df.results.esnp$color == "navy"],
                           name = "miRNA -log10(P)",
                           type = "p",
                           baseline=-log10(high.p.val.thresh),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(gr.esnp$neglog10p)),
                           col="navy")

    overlayTrack.esnp <- OverlayTrack(trackList = list(navyTrack,
                                                       lightblueTrack,
                                                       greenTrack,
                                                       orangeTrack,
                                                       redTrack,
                                                       purpleTrack))



    cat("Adding mRNA track\n")

    # granges for this ensg overlap
    gr.ensg <- GRanges(seqnames = df.blood.results.emir$seqnames,
                       ranges = IRanges(start = df.blood.results.emir$start,
                                        end = df.blood.results.emir$start,
                                        names = df.blood.results.emir$snpID),
                       neglog10p = -log10(df.blood.results.emir$Pval))

    # purple plot
    if (sum(df.blood.results.emir$color == "purple") == 0) {
        purpleTrack <- DataTrack()
    } else {
        purpleTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "purple"],
                                 name = "mRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.ensg$neglog10p)),
                                 col="purple",
                                 cex=1.25,
                                 pch=18)
    }

    # red plot
    if (sum(df.blood.results.emir$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "red"],
                              name = "mRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.ensg$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.blood.results.emir$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "orange"],
                                 name = "mRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.ensg$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.blood.results.emir$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "green"],
                                name = "mRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.ensg$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.blood.results.emir$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "lightblue"],
                                    name = "mRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.ensg$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    if (sum(df.blood.results.emir$color == "navy") == 0) {
        navyTrack <- DataTrack(name = "blood -log10(P)")
    } else {
        navyTrack <- DataTrack(gr.ensg[df.blood.results.emir$color == "navy"],
                               name = "blood -log10(P)",
                               type = "p",
                               legend=FALSE,
                               ylim=c(0,max(gr.ensg$neglog10p)),
                               col="navy")
    }

    overlayTrack.ensg <- OverlayTrack(trackList = list(navyTrack,
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

    pdf(output.file, width = 12, height = 5)
    plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack.esnp, overlayTrack.ensg),
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
               main = paste(df.candidates$plot_name[i]),
               cex.main=1)
    dev.off()




}


