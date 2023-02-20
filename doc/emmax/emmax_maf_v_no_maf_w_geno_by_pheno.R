# compare 1kg v topmed emmax association test

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(lattice)
library(reshape2)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/maf_v_no_maf_geno_pheno")
#dir.pngs <- here("doc/emmax/pngs/")

# INPUT FILES ##########################################################################################################
# GRangesList containing SNP-miR associations within 1MB of each miR
grl.hardcall.rds <- here("results/emmax/association_results/20191101_mirQTL_hardcall/compiled/20191101_mirQTL_topmed_hardcall_GRangesList.rds")
grl.dosage.rds <- here("results/emmax/association_results/20191101_mirQTL_dosage/compiled/20191101_mirQTL_topmed_dosage_GRangesList.rds")
# Dataframe containing the association results
df.hardcall.rds <- here("results/emmax/association_results/20191101_mirQTL_hardcall/compiled/20191101_mirQTL_topmed_hardcall_dataFrame.rds")
df.dosage.rds <- here("results/emmax/association_results/20191101_mirQTL_dosage/compiled/20191101_mirQTL_topmed_dosage_dataFrame.rds")

# genes (mirs) used in this analysis
genes.topmed.hg38.txt <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/20191101_mirQTL_genes_hg38.txt")

# phenotype files directory
pheno.dir <- here("results/emmax/phenotype_files/20191101_mirQTL_VST_miRNA_expression/")

# transposed genotype data directory
hardcall.tfile.dir <- here("results/emmax/tfiles/hardcall/")
dosage.tfile.dir <- here("results/emmax/tfiles/dosage/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")

FDR.THRESHOLD <- 0.05

# Import SNP Associations ##############################################################################################
df.hardcall <- readRDS(df.hardcall.rds)
df.dosage <- readRDS(df.dosage.rds)

# Import Genes #########################################################################################################
df.genes <- read_tsv(genes.topmed.hg38.txt, col_names = c("uniqename", "chr", "start", "end"))

# Format Tables ########################################################################################################

df.hardcall %<>% as_tibble()
df.hardcall$A1.count <- (df.hardcall$A1.HOM.count * 2) + df.hardcall$HET.count
df.hardcall$A2.count <- (df.hardcall$A2.HOM.count * 2) + df.hardcall$HET.count
df.hardcall$A1.freq <- df.hardcall$A1.count / (df.hardcall$A1.count + df.hardcall$A2.count)

df.dosage %<>% as_tibble()
df.dosage$ALT.freq <- df.dosage$ALT_CTS / df.dosage$OBS_CT
df.dosage$REF.freq <- 1 - df.dosage$ALT.freq
df.dosage$MAF <- ifelse(df.dosage$ALT.freq < df.dosage$REF.freq, df.dosage$ALT.freq, df.dosage$REF.freq)

# Find emiRs ###########################################################################################################

# MAF > 0.01 and (HOM.MINOR != 1 & HET > 1)
# use hardcall data to filter variants in the dosage data
df.hardcall %>%
    dplyr::filter(A1.freq > 0.01) %>%
    dplyr::filter(A1.HOM.count != 1 & HET.count > 1) -> df.hardcall.filt

df.dosage %>%
    dplyr::filter(SNP %in% df.hardcall.filt$SNP) -> df.dosage.maf

# NO MAF: HOM.MINOR != 1 & HET > 1
# use hardcall data to filter variants in the dosage data
df.hardcall %>%
    dplyr::filter(A1.HOM.count != 1 & HET.count > 1) -> df.hardcall.filt

df.dosage %>%
    dplyr::filter(SNP %in% df.hardcall.filt$SNP) -> df.dosage.nomaf

rm(df.hardcall.filt)

# adjusted p-value
df.dosage.maf$P.adj <- p.adjust(df.dosage.maf$P, method = "fdr")
df.dosage.nomaf$P.adj <- p.adjust(df.dosage.nomaf$P, method = "fdr")

# get emirs
df.dosage.maf %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(P.adj), min.pval = min(P)) %>%
    dplyr::filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.maf

df.dosage.nomaf %>%
    group_by(UniName) %>%
    summarise(min.pval.adj = min(P.adj), min.pval = min(P)) %>%
    dplyr::filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(UniName) -> emiRs.nomaf

# get nominal p-value threshold
df.dosage.maf %>%
    dplyr::filter(P.adj < FDR.THRESHOLD) %>%
    top_n(1, P) %>%
    pull(P) -> p.value.nominal.maf

df.dosage.nomaf %>%
    dplyr::filter(P.adj < FDR.THRESHOLD) %>%
    top_n(1, P) %>%
    pull(P) -> p.value.nominal.nomaf

p.value.nominal.nomaf <- p.value.nominal.nomaf[1]

# Plot ##########

emiRs.maf.only <- emiRs.maf[!emiRs.maf %in% emiRs.nomaf]
emiRs.nomaf.only <- emiRs.nomaf[!emiRs.nomaf %in% emiRs.maf]

emiRs.combined <- unique(c(emiRs.maf, emiRs.nomaf))

df.emirs <- dplyr::filter(df.genes, uniqename %in% emiRs.combined)

# loop over each emir

#genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

for (i in 1:length(df.emirs$uniqename)) {
    #for (i in 294:306) {
    cat(df.emirs$uniqename[i], "\n")

    mir <- df.emirs$uniqename[i]
    chr <- df.emirs$chr[i]

    # pdf directory
    if (mir %in% emiRs.maf.only) {
        filename <- file.path(dir.pdfs, paste0("/maf_only/", mir, ".pdf"))
        filename.box.maf <- file.path(dir.pdfs, paste0("/maf_only/", mir, "_boxplot_maf.pdf"))
        filename.box.nomaf <- file.path(dir.pdfs, paste0("/maf_only/", mir, "_boxplot_nomaf.pdf"))
    } else if (mir %in% emiRs.nomaf.only) {
        filename <- file.path(dir.pdfs, paste0("/nomaf_only/", mir, ".pdf"))
        filename.box.maf <- file.path(dir.pdfs, paste0("/nomaf_only/", mir, "_boxplot_maf.pdf"))
        filename.box.nomaf <- file.path(dir.pdfs, paste0("/nomaf_only/", mir, "_boxplot_nomaf.pdf"))
    } else {
        filename <- file.path(dir.pdfs, paste0("/both/", mir, ".pdf"))
        filename.box.maf <- file.path(dir.pdfs, paste0("/both/", mir, "_boxplot_maf.pdf"))
        filename.box.nomaf <- file.path(dir.pdfs, paste0("/both/", mir, "_boxplot_nomaf.pdf"))
    }

    # phenotype file
    pheno.file <- paste0(pheno.dir, chr, "/", chr, ".mirQTL.", mir, ".pheno")
    # hardcall genotype file
    geno.hardcall.file <- paste0(hardcall.tfile.dir, chr, "/", chr, ".topmed.hardcall.r2g03.mirQTL.", mir, ".traw")
    geno.dosage.file <- paste0(dosage.tfile.dir, chr, "/", chr, ".topmed.dosage.r2g03.mirQTL.", mir, ".traw")

    # import phenotype data
    df.pheno <- read_delim(pheno.file, delim = " ", col_names = c("DonorID", "DNAID", "Exprs"))
    df.pheno$ID <- paste(df.pheno$DonorID, df.pheno$DNAID, sep = "_")

    # import hardcall genotypes
    df.geno.hardcall <- read_tsv(geno.hardcall.file)
    df.geno.dosage <- read_tsv(geno.dosage.file)

    all(df.pheno$ID == colnames(df.geno.hardcall[7:223]))
    all(df.pheno$ID == colnames(df.geno.dosage[7:223]))

    # subset variant association data
    df.dosage.maf %>%
        dplyr::filter(UniName == !!mir) -> df.dosage.maf.sub

    df.dosage.nomaf %>%
        dplyr::filter(UniName == !!mir) -> df.dosage.nomaf.sub

    # select top 5 variants
    df.dosage.maf.sub %>%
        dplyr::filter(P.adj < FDR.THRESHOLD) %>%
        top_n(5, -P) %>%
        pull(SNP) -> top.five.vars.maf

    df.dosage.nomaf.sub %>%
        dplyr::filter(P.adj < FDR.THRESHOLD) %>%
        top_n(5, -P) %>%
        pull(SNP) -> top.five.vars.nomaf

    top.five.vars.maf <- top.five.vars.maf[1:5]
    top.five.vars.nomaf <- top.five.vars.nomaf[1:5]

    # plot top 5 with maf filter
    pdf(filename.box.maf, height = 4, width = 8)

    df.geno.hardcall %>%
        dplyr::filter(SNP %in% top.five.vars.maf) %>%
        dplyr::select(SNP, starts_with("D")) %>%
        melt() %>%
        left_join(dplyr::select(df.pheno, ID, Exprs), by = c("variable" = "ID")) %>%
        mutate(value = factor(value, levels = c(0,1,2))) %>%
        ggplot(aes(x = value, y = Exprs)) +
        geom_boxplot(outlier.color = NA) +
        geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
        facet_wrap(~SNP, nrow = 1) +
        scale_x_discrete(drop = FALSE) +
        labs(title = paste(mir, "top 5 snps, filter: maf > 0.01"))

    df.geno.dosage %>%
        dplyr::filter(SNP %in% top.five.vars.maf) %>%
        dplyr::select(SNP, starts_with("D")) %>%
        melt() %>%
        left_join(dplyr::select(df.pheno, ID, Exprs), by = c("variable" = "ID")) %>%
        ggplot(aes(x = value, y = Exprs)) +
        geom_point(size = 1, alpha = 0.5) +
        scale_x_continuous(limits = c(0,2)) +
        facet_wrap(~SNP, nrow = 1) +
        labs(title = paste(mir, "top 5 snps, filter: maf > 0.01"))

    dev.off()

    # plot top 5 with no maf filter
    pdf(filename.box.nomaf, height = 4, width = 8)

    df.geno.hardcall %>%
        dplyr::filter(SNP %in% top.five.vars.nomaf) %>%
        dplyr::select(SNP, starts_with("D")) %>%
        melt() %>%
        left_join(dplyr::select(df.pheno, ID, Exprs), by = c("variable" = "ID")) %>%
        mutate(value = factor(value, levels = c(0,1,2))) %>%
        ggplot(aes(x = value, y = Exprs)) +
        geom_boxplot(outlier.color = NA) +
        geom_jitter(width = 0.1, size = 1, alpha = 0.5) +
        facet_wrap(~SNP, nrow = 1) +
        scale_x_discrete(drop = FALSE) +
        labs(title = paste(mir, "top 5 snps, no maf filter"))

    df.geno.dosage %>%
        dplyr::filter(SNP %in% top.five.vars.nomaf) %>%
        dplyr::select(SNP, starts_with("D")) %>%
        melt() %>%
        left_join(dplyr::select(df.pheno, ID, Exprs), by = c("variable" = "ID")) %>%
        ggplot(aes(x = value, y = Exprs)) +
        geom_point(size = 1, alpha = 0.5) +
        scale_x_continuous(limits = c(0,2)) +
        facet_wrap(~SNP, nrow = 1) +
        labs(title = paste(mir, "top 5 snps, no maf filter1"))

    dev.off()

    MB.WINDOW <- 1e6

    startRange <- NULL
    endRange <- NULL

    # genomic range for SNPs associated with this miRNA
    startRange <- df.emirs$start[i] - MB.WINDOW
    endRange <- df.emirs$end[i] + MB.WINDOW

    if (startRange < 0) {
        startRange <- 0
    }

    # MAF Data ##############################

    gr.maf <- GRanges(seqnames = df.dosage.maf.sub$CHR,
                      ranges = IRanges(start = df.dosage.maf.sub$BP.hg38,
                                       end = df.dosage.maf.sub$BP.hg38,
                                       names = df.dosage.maf.sub$SNP),
                      neglog10p = -log10(df.dosage.maf.sub$P))

    # No MAF Data #####################################

    gr.nomaf <- GRanges(seqnames = df.dosage.nomaf.sub$CHR,
                        ranges = IRanges(start = df.dosage.nomaf.sub$BP.hg38,
                                         end = df.dosage.nomaf.sub$BP.hg38,
                                         names = df.dosage.nomaf.sub$SNP),
                        neglog10p = -log10(df.dosage.nomaf.sub$P))



    # Make Tracks

    col.maf <- "black"
    if (mir %in% emiRs.maf) {
        col.maf <- "darkgreen"
    } else {
        col.maf <- "darkred"
    }

    if (length(gr.maf) > 0) {
        dTrack.maf <- DataTrack(gr.maf,
                                name = "MAF Filter",
                                type = "p",
                                baseline=-log10(p.value.nominal.maf),
                                ylim=c(0,max(c(gr.maf$neglog10p, p.value.nominal.maf))),
                                col.baseline="black",
                                lty.baseline=2,
                                lwd.baseline=1,
                                legend=FALSE,
                                col=col.maf)
    } else {
        dTrack.maf <- DataTrack()
    }

    col.nomaf <- "black"
    if (mir %in% emiRs.nomaf) {
        col.nomaf <- "darkgreen"
    } else {
        col.nomaf <- "darkred"
    }

    dTrack.nomaf <- DataTrack(gr.nomaf,
                              name = "No MAF Filter",
                              type = "p",
                              baseline=-log10(p.value.nominal.nomaf),
                              ylim=c(0,max(c(gr.nomaf$neglog10p, p.value.nominal.nomaf))),
                              col.baseline="black",
                              lty.baseline=2,
                              lwd.baseline=1,
                              legend=FALSE,
                              col=col.nomaf)

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    # # Gene region track
    # cat("Adding gene region track\n")
    #
    # grtrack <- BiomartGeneRegionTrack(genome = "hg19",
    #                                   biomart = genemart,
    #                                   start = startRange,
    #                                   end = endRange,
    #                                   chromosome = chr,
    #                                   showId = TRUE,
    #                                   geneSymbols = TRUE,
    #                                   rotation.title = 90,
    #                                   transcriptAnnotation = "symbol",
    #                                   name = "ENSEMBL",
    #                                   collapseTranscripts = TRUE,
    #                                   just.group="below")
    #
    # mirTrack <- AnnotationTrack(gr.mirs,
    #                             genome = "hg38",
    #                             chromosome = chr,
    #                             name = "miRNA",
    #                             rotation.title = 90,
    #                             id = gr.mirs$Name,
    #                             showFeatureId = TRUE,
    #                             cex = .6,
    #                             fontcolor.item = "black",
    #                             shape = "box")

    pdf(filename)

    plotTracks(list(itrack, gtrack, dTrack.maf, dTrack.nomaf),
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
               just.group="below",
               main = mir)

    dev.off()
}

