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
dir.output <- paste0(dir.pdf, "fetal_brain_colocalizations_updated/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# mirQTL eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
esnps.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# mirQTL emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# fetal brain mRNA-eQTL clumped with ld buddies
clumped.ldbuddies.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs_and_LDbuddies.rds")
index.snps.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_index_SNPs.rds")

# fetal brain raw data by chr dir
dir.mrna.eqtl.raw <- here("results/external_data/fetal_brain_mRNA-eQTL/raw/")

# summarized experiment with mirna vst expression values used in this analysis
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# table of colocalizations
coloc.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalizations.rds")

# table of hosts
hosts.rds <- here("results/rdata_files/tmp_20200224_miRNA_host_expression_corr.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import ###############################################################################################################
# mirQTL eqtls
df.eqtls <- readRDS(eqtls.dataframe.rds)

# mrna eqtls
#df.mrna.eqtls <- readRDS(index.snps.rds)

# mrna eqtl ld buddies
df.mrna.ldbuddies <- readRDS(clumped.ldbuddies.rds)

# r2 > 0.8 ld buddies
df.mrna.high.ldbuddies <- dplyr::filter(df.mrna.ldbuddies, R2 > 0.8)

# all variants mirQTL
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# expression data
vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd.mirna)]

# set seqnames for gene expression to UCSC
seqlevelsStyle(rse.gene) <- "UCSC"
genome(rse.gene) <- "hg38"
# remove non-standard chroms
seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

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

# host genes
df.hosts <- readRDS(hosts.rds)

# Colocalizations ######################################################################################################

df.coloc <- readRDS(coloc.rds)

# coloc label
df.coloc$COLOC <- paste(df.coloc$eQTL.mirQTL, df.coloc$SNP_ENSG.mQTL, sep = "-")

length(unique(df.coloc$COLOC))

unique.colocs <- unique(df.coloc$COLOC)


gene_names <- rowData(vsd.gene)$gene_name
names(gene_names) <- rowData(vsd.gene)$gene_id

# Plot Colocalizations #################################################################################################


for (i in 1:length(unique.colocs)) {
    printMessage(paste("Plotting", i, "of", nrow(df.coloc), ":", df.coloc$eQTL[i]))

    df.this.coloc <- dplyr::filter(df.coloc, COLOC == unique.colocs[i])

    esnp <- df.coloc$eSNP[i]
    emir <- df.coloc$emiR[i]
    eqtl <- df.coloc$eQTL[i]

    chr <- df.coloc$SNP.CHR[i]
    esnp.bp <- df.coloc$SNP.BP.hg38[i]

    ensg <- df.coloc$ENSG[i]
    mrna.esnp <- df.coloc$SNP_A[i]

    gene_name <- gene_names[ensg]

    # output plot file
    output.file <- paste0(dir.output, df.coloc$COLOC[i], ".pdf")

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

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # eSNP window
    buffer.bp <- 1e6

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


    # load mrna eqtl data
    df.results.mrna <- readRDS(paste0(dir.mrna.eqtl.raw, chr, "_fetal_brain_mRNA-eQTL_raw.rds"))

    # subset for this gene
    df.results.ensg <- dplyr::filter(df.results.mrna, ENSG == ensg, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)
    df.mrna.ldbuddies.ensg <- dplyr::filter(df.mrna.ldbuddies, ENSG == ensg, BP.hg38_B >= startRange.esnp & BP.hg38_B <= endRange.esnp)

    df.results.ensg <- left_join(df.results.ensg, df.mrna.ldbuddies.ensg, by = c("SNP" = "SNP_B"))

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
                       ranges = IRanges(start = df.results.ensg$BP.hg38,
                                        end = df.results.ensg$BP.hg38,
                                        names = df.results.ensg$SNP),
                       neglog10p = -log10(df.results.ensg$P.x))

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

    pdf(output.file, width = 7, height = 5)
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
               main = paste(eqtl, paste(df.coloc$SNP_ENSG[i], gene_name), sep = "\n"),
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
