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
#dir.pdfs <- here("doc/mirQTL_plots/pdfs/")

# output folder for combined plots
dir.output <- here("doc/mirQTL_plots/pdfs/mirQTL_dot_geno_by_pheno_conditional_compiled/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
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

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/sample_genotypes/")

# Directory for LD at each index SNP
dir.ld <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Covariate matrix autosomes
#cov.mat.auto.file <- here("results/emmax/covariates/20191204_chrAutosomes.mirQTLor.columnLabeled.cov")
# Covariate matrix x chromosome
#cov.mat.x.file <- here("results/emmax/covariates/20191204_chrX.mirQTLor.columnLabeled.cov")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

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

# expression data
vsd <- readRDS(vsd.rds)

# Plot eQTLs ###########################################################################################################

genemart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl")

ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

trouble.ind <- NULL

# loop over each eQTL
for (i in 1:nrow(df.eqtls)) {
    #for (i in unique(trouble.ind)) {

    #for (i in 1:5) {
    printMessage(paste("Plotting", i, "of", nrow(df.eqtls), ":", df.eqtls$eQTL[i]))

    esnp <- df.eqtls$eSNP[i]
    emir <- df.eqtls$emiR[i]
    eqtl <- df.eqtls$eQTL[i]

    rank <- df.eqtls$DEGREE[i]
    sig <- df.eqtls$SIGNIFICANCE[i]

    chr <- df.eqtls$SNP.CHR[i]
    esnp.bp <- df.eqtls$SNP.BP.hg38[i]

    # reference allele
    allele.ref <- df.results$REF[match(eqtl, df.results$UniName_SNP)]
    # alt allele (effect allele)
    allele.alt <- df.results$ALT[match(eqtl, df.results$UniName_SNP)]
    # minor allele
    allele.minor <- df.results$A1[match(eqtl, df.results$UniName_SNP)]
    # major allele
    allele.major <- df.results$A2[match(eqtl, df.results$UniName_SNP)]
    # effect allele (for dosage based associations)
    effect.allele<- df.results$EFFECT.ALLELE[match(eqtl, df.results$UniName_SNP)]

    if (effect.allele != allele.ref) {
        stop("Effect allele != reference allele!")
    }


    # check for flipped genotypes: plink1 stores genotypes as Minor/Major (coded as number of minor alleles)

    # ld file: conditional.mirQTLor.index.chrX:84834596:G:T.ld
    ld.file <- paste0(dir.ld, "conditional.mirQTLor.index.", esnp, ".ld")

    # genotype files: chrX.hardcall.prefiltered.mirQTLor.chrX:109505804:TA:T.raw
    geno.hardcall.file <- paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".raw")
    geno.dosage.file <- paste0(dir.genotypes, chr, ".dosage.prefiltered.mirQTLor.", esnp, ".raw")

    # output plot file
    output.file <- paste0(dir.output, eqtl, "_rank", rank, "_", sig, ".pdf")

    # Genotype x Phenotype Plot ########################################################################################

    # subset expression data
    df.expression <- data.frame(expression = assays(vsd[emir,])[["VST"]][1,],
                                expression.corrected = assays(vsd[emir,])[["VST.RESIDUAL"]][1,],
                                DonorID = paste0("D", vsd$donor_id),
                                stringsAsFactors = FALSE)
    df.expression$RNAID <- rownames(df.expression)

    # load genotype data
    df.hardcall.genotypes <- read_delim(geno.hardcall.file,
                                        delim = " ",
                                        col_types = "cccciiii")
    counted.allele.hardcall <- strsplit(colnames(df.hardcall.genotypes)[7], "_")[[1]][2]

    df.hardcall.genotypes %<>%
        dplyr::select(DonorID = FID, DNAID = IID, Sex = SEX, Genotype_Hard = 7)

    # hardcall genotype labels
    if (counted.allele.hardcall == effect.allele) {
        hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                paste(allele.alt, allele.ref, sep = "/"),
                                paste(allele.ref, allele.ref, sep = "/"))
        rev.flag <- FALSE
    } else if (counted.allele.hardcall == allele.alt) {
        hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                paste(allele.alt, allele.ref, sep = "/"),
                                paste(allele.ref, allele.ref, sep = "/"))
        rev.flag <- TRUE
    } else {
        stop("Something is wrong with genotype labels!")
    }


    df.dosage.genotypes <- read_tsv(geno.dosage.file,
                                    col_types = "cccciinn")
    counted.allele.dosage <- strsplit(colnames(df.dosage.genotypes)[7], "_")[[1]][2]
    # df.dosage.genotypes <- read_tsv(geno.dosage.file,
    #                                 skip = 1,
    #                                 col_names = c("DonorID", "DNAID", "PAT", "MAT", "Sex", "PHE", "Genotype_Dose", "Geno_Dose_Het"),
    #                                 col_types = "cccciinn")
    df.dosage.genotypes %<>%
        dplyr::select(DonorID = FID, DNAID = IID, Sex = SEX, Genotype_Dose = 7)

    # dosage genotype check
    if (counted.allele.dosage == effect.allele) {
        dosage.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                              paste(allele.alt, allele.ref, sep = "/"),
                              paste(allele.ref, allele.ref, sep = "/"))
    } else {
        stop("Something is wrong with genotype labels!")
    }

    df.hardcall.genotypes %>%
        dplyr::select(DonorID, Genotype_Hard) %>%
        dplyr::left_join(df.dosage.genotypes, by = "DonorID") %>%
        dplyr::select(DonorID, DNAID, Sex, Genotype_Hard, Genotype_Dose)-> df.genotypes

    rm(df.dosage.genotypes, df.hardcall.genotypes)

    # format genotype hardcalls as factor
    if (rev.flag) {
        df.genotypes$Genotype_Hard <- factor(df.genotypes$Genotype_Hard, levels = rev(c(0, 1, 2)), ordered = TRUE, labels = hardcall.genotypes)
    } else {
        df.genotypes$Genotype_Hard <- factor(df.genotypes$Genotype_Hard, levels = c(0, 1, 2), ordered = TRUE, labels = hardcall.genotypes)
    }
    # format sex as factor
    df.genotypes$Sex <- factor(df.genotypes$Sex, levels = c(1, 2), labels = c("M", "F"))

    # combine expression and genotypes
    df.geno.pheno <- dplyr::left_join(df.genotypes, df.expression, by = "DonorID")

    # plot geno x pheno
    df.geno.pheno %>%
        ggplot(aes(x = Genotype_Hard, y = expression)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 1) +
        #geom_jitter(mapping = aes(color = Sex), width = 0.1, size = 1) +
        scale_x_discrete(drop = FALSE) +
        labs(y = "VST Expression",
             x = paste0("REF(effect):", allele.ref, " ALT:", allele.alt)) -> p1

    df.geno.pheno %>%
        ggplot(aes(x = Genotype_Dose, y = expression)) +
        geom_point(size = 1) +
        scale_x_continuous(limits = c(0,2)) +
        labs(y = "VST Expression",
             x = "REF Dosage") -> p2

    df.geno.pheno %>%
        ggplot(aes(x = Genotype_Hard, y = expression.corrected)) +
        geom_boxplot(outlier.shape = NA) +
        geom_jitter(width = 0.1, size = 1) +
        #geom_jitter(mapping = aes(color = Sex), width = 0.1, size = 1) +
        scale_x_discrete(drop = FALSE) +
        labs(y = "VST Expression (residual)",
             x = paste0("REF(effect):", allele.ref, " ALT:", allele.alt)) -> p3

    df.geno.pheno %>%
        ggplot(aes(x = Genotype_Dose, y = expression.corrected)) +
        geom_point(size = 1) +
        geom_smooth(method = "lm") +
        scale_x_continuous(limits = c(0,2)) +
        labs(y = "VST Expression (residual)",
             x = "REF Dosage") -> p4


    # Dot Plots ########################################################################################################

    # subset results for this miRNA
    df.results %>%
        dplyr::filter(UniName == emir) %>%
        dplyr::filter(DEGREE == rank) -> df.results.emir

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # ld for this esnp
    suppressWarnings(
        df.ld <- read_table(ld.file))
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

    # Dot Plot: emiR ###################################################################################################
    # plotting window
    #emir.WINDOW <- 1e6

    startRange.emir <- min.bp
    endRange.emir <- max.bp

    # # genomic range for SNPs associated with this miRNA
    # startRange <- start(rowRanges(vsd)[emir]) - emir.WINDOW
    # endRange <- end(rowRanges(vsd)[emir]) + emir.WINDOW

    # granges for this plot
    gr.emir <- GRanges(seqnames = df.results.emir$CHR,
                       ranges = IRanges(start = df.results.emir$BP.hg38,
                                        end = df.results.emir$BP.hg38,
                                        names = df.results.emir$SNP),
                       neglog10p = -log10(df.results.emir$P))

    # data track for this plot
    cat("Adding data track\n")
    dTrack <- DataTrack(gr.emir,
                        name = "-log10(P)",
                        type = "p",
                        baseline=c(-log10(high.p.val.thresh), -log10(med.p.val.thresh), -log10(low.p.val.thresh)),
                        ylim=c(0,max(gr.emir$neglog10p)),
                        col.baseline="black",
                        lty.baseline=2,
                        lwd.baseline=1,
                        legend=FALSE,
                        col="black",
                        cex=0.5)

    # Ideogram track
    cat("Adding ideogram track\n")

    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr, bands = ideo)

    # Genome axis track
    cat("Adding genome axis track\n")

    gtrack <- GenomeAxisTrack(exponent = 0)

    # Gene region track
    cat("Adding gene region track\n")

    grtrack.emir <- BiomartGeneRegionTrack(genome = "hg38",
                                           biomart = genemart,
                                           start = startRange.emir,
                                           end = endRange.emir,
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
    # miRTrack <- BiomartGeneRegionTrack(genome = "hg38",
    #                                    biomart = genemart,
    #                                    start = startRange.emir,
    #                                    end = endRange.emir,
    #                                    chromosome = chr,
    #                                    showId = TRUE,
    #                                    geneSymbols = TRUE,
    #                                    rotation.title = 0,
    #                                    transcriptAnnotation = "symbol",
    #                                    name = "miRBase",
    #                                    collapseTranscripts = "meta",
    #                                    filter=list(with_mirbase = TRUE),
    #                                    cex=0.1)

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


    # highlight track
    highlightTrack <- HighlightTrack(trackList = list(gtrack, grtrack.emir, miRTrack, dTrack),
                                     start = esnp.bp - 1000,
                                     end = esnp.bp + 1000,
                                     chromosome = chr)

    # plot title
    plot.title <- paste0(eqtl,
                         " / BETA: ", signif(df.results.emir$BETA[match(eqtl, df.results.emir$UniName_SNP)], 3),
                         " / P: ", signif(df.results.emir$P[match(eqtl, df.results.emir$UniName_SNP)], 3),
                         " / RANK: ", rank,
                         " / ", sig)

    # new multi plot grid
    pdf(output.file, width = 12, height = 6)

    grid.newpage()

    pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 5, heights = unit(c(0.5, 5, 5), "null"))))

    grid.text(plot.title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:5))

    pushViewport(viewport(layout.pos.row = 2, layout.pos.col = 1:3))

    tryCatch(
        {
            plotTracks(list(itrack, highlightTrack),
                       chromosome = chr,
                       from = startRange.emir,
                       to = endRange.emir,
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
                       just.group = "right")
            popViewport()
        },
        error = function(mess) {
            message(mess)
            trouble.ind <<- c(trouble.ind, i)
        }
    )

    # Dot Plot: eSNP ###################################################################################################

    # window: eSNP to emiR + 1e4
    emir.bp <- start(rowRanges(vsd)[emir])
    eqtl.WINDOW <- 5e5
    # if (esnp.bp < emir.bp) {
    #     # esnp upstream emir
    #     startRange.esnp <- esnp.bp - eqtl.WINDOW
    #     endRange.esnp <- emir.bp + eqtl.WINDOW
    # } else {
    #     # esnp downstream emir
    #     startRange.esnp <- emir.bp - eqtl.WINDOW
    #     endRange.esnp <- esnp.bp + eqtl.WINDOW
    # }

    startRange.esnp <- esnp.bp - eqtl.WINDOW
    endRange.esnp <- esnp.bp + eqtl.WINDOW

    # if outside min/max coords, adjust
    if (startRange.esnp < min.bp) {
        startRange.esnp <- min.bp
    }
    if (endRange.esnp > max.bp) {
        endRange.esnp <- max.bp
    }

    # subset results dataframe
    df.results.esnp <- dplyr::filter(df.results.emir, BP.hg38 >= startRange.esnp & BP.hg38 <= endRange.esnp)

    # granges for this plot
    gr.esnp <- GRanges(seqnames = df.results.esnp$CHR,
                       ranges = IRanges(start = df.results.esnp$BP.hg38,
                                        end = df.results.esnp$BP.hg38,
                                        names = df.results.esnp$SNP),
                       neglog10p = -log10(df.results.esnp$P))

    # purple plot
    purpleTrack <- DataTrack(gr.esnp[esnp],
                             name = "-log10(P)",
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
                              name = "-log10(P)",
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
                                 name = "-log10(P)",
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
                                name = "-log10(P)",
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
                                    name = "-log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.esnp$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.esnp[df.results.esnp$color == "navy"],
                           name = "-log10(P)",
                           type = "p",
                           baseline=c(-log10(high.p.val.thresh), -log10(med.p.val.thresh), -log10(low.p.val.thresh)),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(gr.esnp$neglog10p)),
                           col="navy")

    overlayTrack <- OverlayTrack(trackList = list(navyTrack,
                                                  lightblueTrack,
                                                  greenTrack,
                                                  orangeTrack,
                                                  redTrack,
                                                  purpleTrack))

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
    # miRTrack.esnp <- BiomartGeneRegionTrack(genome = "hg38",
    #                                         biomart = genemart,
    #                                         start = startRange.esnp,
    #                                         end = endRange.esnp,
    #                                         chromosome = chr,
    #                                         showId = TRUE,
    #                                         geneSymbols = TRUE,
    #                                         rotation.title = 0,
    #                                         transcriptAnnotation = "symbol",
    #                                         name = "miRBase",
    #                                         collapseTranscripts = "meta",
    #                                         filter=list(with_mirbase = TRUE),
    #                                         cex=0.1)

    miRTrack.esnp <- AnnotationTrack(gr.mirna.primary,
                                     genome = "hg38",
                                     chromosome = chr,
                                     name = "miRNA",
                                     rotation.title = 0,
                                     #id = gr.mirna.primary$Name,
                                     group = gr.mirna.primary$Name,
                                     #groupAnnotation = "group",
                                     stacking = "dense",
                                     featureAnnotation="group",
                                     #showID = TRUE,
                                     showFeatureId = TRUE,
                                     fontcolor.item = "black",
                                     shape = "box",
                                     cex = 0.7,
                                     just.id = "left")

    pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:3))

    tryCatch(
        {
            plotTracks(list(gtrack, overlayTrack, grtrack.esnp, miRTrack.esnp),
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
                       just.group="right",
                       add=TRUE)
            popViewport()
        },
        error = function(mess) {
            message(mess)
            trouble.ind <<- c(trouble.ind, i)
        }
    )


    print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
    print(p2, vp = viewport(layout.pos.row = 3, layout.pos.col = 4))

    print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
    print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 5))

    dev.off()
}



