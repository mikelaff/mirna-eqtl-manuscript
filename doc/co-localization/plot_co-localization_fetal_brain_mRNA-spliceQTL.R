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
library(DESeq2)
library(mikelaffr)


# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/co-localization/pdfs/")

# output folder for combined plots
dir.output <- paste0(dir.pdf, "fetal_brain_mrna-spliceqtl_colocalizations/2e5/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
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

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
#nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# Directory for LD at each index SNP
mirQTL.dir.ld <- here("results/co-localization/fetal_brain_sQTL/ld/mirQTL/")
mQTL.dir.ld <- here("results/co-localization/fetal_brain_sQTL/ld/mQTL/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# fetal brain raw data by chr dir
dir.mrna.eqtl.raw <- "/proj/steinlab/projects/R00/eQTLanalysis/Raw_data/splicing/bulk/"

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# table of overlaps
overlap.rds <- here("results/co-localization/fetal_brain_sQTL/fetal_brain_sQTL_mirQTL_overlaps_r2at0.8.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import ###############################################################################################################
# overlaps
df.overlaps <- readRDS(overlap.rds)
df.overlaps %<>%
    mutate(overlap = paste(eQTL.mirQTL, intron.mQTL, sep = "+"),
           DEGREE.mirQTL = factor(DEGREE.mirQTL, levels = c("primary", "secondary", "tertiary", "quarternary", "quinary"), labels = c(1,2,3,4,5)))

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



# # expression data
# vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
#rse.gene <- readRDS(rse.gene.rds)

# # subset samples to those in vsd.mirna
# rse.gene <- rse.gene[,colnames(vsd.mirna)]

# # set seqnames for gene expression to UCSC
# genome(rse.gene) <- "hg38"
# # remove non-standard chroms
# seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS
# seqlevelsStyle(rse.gene) <- "UCSC"

# remove biotype = mirna
#rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# # normalize expression (VST)
# dds.gene <- DESeqDataSet(rse.gene, design = ~1)
# vsd.gene <- vst(dds.gene)

# nominal p-value for plotting
#nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# # host genes
# df.hosts <- readRDS(hosts.rds)

#gene_names <- rowData(rse.gene)$gene_name
#names(gene_names) <- rowData(rse.gene)$gene_id

# Plot Colocalizations #################################################################################################


for (i in 1:nrow(df.overlaps)) {
    printMessage(paste("Plotting", i, "of", nrow(df.overlaps), ":", df.overlaps$eQTL.mirQTL[i]))

    esnp <- df.overlaps$eSNP.mirQTL[i]
    emir <- df.overlaps$emiR.mirQTL[i]
    eqtl <- df.overlaps$eQTL.mirQTL[i]

    chr <- df.overlaps$SNP.CHR.mirQTL[i]
    esnp.bp <- df.overlaps$SNP.BP.hg38.mirQTL[i]

    ensg <- df.overlaps$ensemblID.mQTL[i]
    mrna.esnp <- df.overlaps$snp.mQTL[i]

    gene_name <- df.overlaps$gene.mQTL[i]

    nom.p.val <- as.numeric(df.overlaps$NOM.P.VALUE.THRESH.mirQTL[i])
    rank <- df.overlaps$DEGREE.mirQTL[i]

    # output plot file
    output.file <- paste0(dir.output, df.overlaps$overlap[i], ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir, DEGREE == rank)

    # ld file: mirQTL.overlap.index.chr10:104260640:A:C.ld
    ld.file <- paste0(mirQTL.dir.ld, "mirQTL.overlap.index.", esnp, ".ld")
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

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # eSNP window
    buffer.bp <- 2e5

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


    # load mrna splice qtl data
    df.results.mrna <- read_tsv(paste0(dir.mrna.eqtl.raw, chr, "/", df.overlaps$intron.mQTL[i], ".ps"), col_names = c("SNP", "BETA", "P"))

    df.results.mrna %<>%
        mutate(CHR = sapply(strsplit(SNP, ":"), `[`, 1),
               BP = as.integer(sapply(strsplit(SNP, ":"), `[`, 2)))

    # subset for this gene
    df.results.ensg <- dplyr::filter(df.results.mrna, BP >= startRange.esnp & BP <= endRange.esnp)

    # ld file: mirQTL.overlap.index.chr10:104260640:A:C.ld
    ld.file <- paste0(mQTL.dir.ld, "mQTL.overlap.index.", esnp, ".ld")
    # ld for this esnp
    df.ld <- read_table(ld.file, col_types = cols())
    df.ld %<>%
        dplyr::select(SNP = SNP_B, R2)

    df.results.ensg %<>%
        left_join(df.ld, by = "SNP")

    # dot color based on LD
    df.results.ensg$color <- "navy"
    df.results.ensg$color[which(df.results.ensg$R2 >= 0.8)] <- "red"
    df.results.ensg$color[which(df.results.ensg$R2 >= 0.6 & df.results.ensg$R2 < 0.8)] <- "orange"
    df.results.ensg$color[which(df.results.ensg$R2 >= 0.4 & df.results.ensg$R2 < 0.6)] <- "green"
    df.results.ensg$color[which(df.results.ensg$R2 >= 0.2 & df.results.ensg$R2 < 0.4)] <- "lightblue"
    df.results.ensg$color[which(df.results.ensg$SNP == mrna.esnp)] <- "purple"


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
                           baseline=-log10(nom.p.val),
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
    gr.ensg <- GRanges(seqnames = df.results.ensg$CHR,
                       ranges = IRanges(start = df.results.ensg$BP,
                                        end = df.results.ensg$BP,
                                        names = df.results.ensg$SNP),
                       neglog10p = -log10(df.results.ensg$P))

    # purple plot
    if (sum(df.results.ensg$color == "purple") == 0) {
        purpleTrack <- DataTrack()
    } else {
        purpleTrack <- DataTrack(gr.ensg[df.results.ensg$color == "purple"],
                                 name = "mRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.ensg$neglog10p)),
                                 col="red")
    }

    # red plot
    if (sum(df.results.ensg$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.ensg[df.results.ensg$color == "red"],
                              name = "mRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.ensg$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.results.ensg$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.ensg[df.results.ensg$color == "orange"],
                                 name = "mRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.ensg$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.results.ensg$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.ensg[df.results.ensg$color == "green"],
                                name = "mRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.ensg$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.results.ensg$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.ensg[df.results.ensg$color == "lightblue"],
                                    name = "mRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.ensg$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.ensg[df.results.ensg$color == "navy"],
                           name = "mRNA -log10(P)",
                           type = "p",
                           baseline=-log10(0.0002957020645),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(gr.ensg$neglog10p)),
                           col="navy")

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
               main = paste(df.overlaps$overlap[i], gene_name),
               cex.main=1)
    dev.off()




}

# Expression Correlation ###############################################################################################

df.coloc$corr <- NA

for (k in 1:nrow(df.coloc)) {

    emiR.expression <- assay(vsd.mirna)[df.coloc$emiR[k],]

    if (!df.coloc$ENSG[k] %in% rownames(vsd.gene)) {
        df.coloc$corr[k] <- 0
        next
    } else {
        ensg.expression <- assay(vsd.gene)[df.coloc$ENSG[k],]
    }
    stopifnot(all(names(ensg.expression) == names(emiR.expression)))

    df.coloc$corr[k] <- cor(emiR.expression, ensg.expression)

}

df.coloc %>%
    ggplot(aes(x = corr)) +
    geom_histogram(bins = 10) +
    labs(title = "emiR-mRNA Expression Correlation for Co-Localized eQTLs")

ggsave(paste0(dir.pdf, "emir-mrna_expression_correlation_for_coloc_eqtls.pdf"))

# Scratch ############

df.tmp <- df.coloc[!duplicated(df.coloc$COLOC),]

#df.tmp <- left_join(df.tmp, df.hosts, by = c("emiR.mirQTL" = "uniName"))

df.coloc.with.host <- dplyr::filter(df.tmp, emiR.mirQTL %in% df.hosts$uniName)
df.coloc.no.host <- dplyr::filter(df.tmp, !emiR.mirQTL %in% df.hosts$uniName)

df.tmp$host <- df.tmp$emiR.mirQTL %in% df.hosts$uniName


df.tmp %>%
    ggplot(aes(x = emiR_eGene_Corr, fill = host)) +
    geom_histogram() +
    facet_wrap(~host) +
    labs(title = "Colocalizations",
         fill = "emiR Has Host") +
    scale_fill_manual(values = c("navy", "darkorange"))

ggsave("~/Desktop/colocs_emiR_eGene_Corr_by_host.pdf")


df.tmp2 <- df.coloc[!duplicated(df.coloc$emiR.mirQTL),]
df.tmp2$host <- df.tmp2$emiR.mirQTL %in% df.hosts$uniName

sum(df.tmp2$host)

df.tmp2 %>%
    ggplot(aes(x = emiR_eGene_Corr, fill = host)) +
    geom_histogram() +
    facet_wrap(~host) +
    labs(title = "emiRs",
         fill = "emiR Has Host") +
    scale_fill_manual(values = c("navy", "darkorange"))

ggsave("~/Desktop/emirs_emiR_eGene_Corr_by_host.pdf")

df.tmp %>%
    left_join(df.hosts, by = c("emiR.mirQTL" = "uniName")) -> tmp3

tmp3 <- tmp3[!is.na(tmp3$gene_id),]
sum(tmp3$ENSG.mQTL == tmp3$gene_id)

tmp3 <- tmp3[!is.na(tmp3$ENSG.mQTL),]

tmp3 %>%
    dplyr::filter(!is.na())
