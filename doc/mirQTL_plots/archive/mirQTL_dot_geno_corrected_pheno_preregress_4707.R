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
dir.pdfs <- here("doc/mirQTL_plots/pdfs/")
dir.pngs <- here("doc/mirQTL_plots/pngs/")

# output folder for combined plots
dir.output <- here("doc/mirQTL_plots/pdfs/mirQTL_dot_geno_corrected_pheno_preregress_4707/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200103_mirQTLor_preregress/ld/")

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


# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# expression data
vsd <- readRDS(vsd.rds)

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# # autosome covariates file
# cov.mat.auto <- read_delim(cov.mat.auto.file, delim = " ")
# cov.mat.auto %<>%
#     dplyr::rename(DonorID = X1, DNAID = X2)
#
# # x chromosome covariates file
# cov.mat.x <- read_delim(cov.mat.x.file, delim = " ")
# cov.mat.x %<>%
#     dplyr::rename(DonorID = X1, DNAID = X2)
#
# # covariates
# batchVars.auto <- c(paste0("PC", 1:10, ".genotype"),
#                     paste0("PC", 1:10, ".expression"),
#                     "PoolPool2", "PoolPool3", "PoolPool4",
#                     "PoolPool5", "PoolPool6", "PoolPool7",
#                     "PoolPool8","PurificationMethodmiRNeasy",
#                     "PurificationMethodmiRNeasy_mini", "SexM", "RIN", "GestationWeek")
# batchVars.x <- c(paste0("PC", 1:10, ".genotype"),
#                  paste0("PC", 1:10, ".expression"),
#                  "PoolPool2", "PoolPool3", "PoolPool4",
#                  "PoolPool5", "PoolPool6", "PoolPool7",
#                  "PoolPool8","PurificationMethodmiRNeasy",
#                  "PurificationMethodmiRNeasy_mini", "RIN", "GestationWeek")
#
# # formula strings
# formula.string.auto <- paste("expression ~", paste(batchVars.auto, collapse = " + "))
# formula.string.x <- paste("expression ~", paste(batchVars.x, collapse = " + "))

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
    printMessage(paste("Plotting", i, "of", nrow(df.eqtls), ":", df.eqtls$eqtl[i]))

    esnp <- df.eqtls$esnp[i]
    emir <- df.eqtls$emir[i]
    eqtl <- df.eqtls$eqtl[i]

    chr <- strsplit(esnp, ":")[[1]][1]
    esnp.bp <- as.integer(strsplit(esnp, ":")[[1]][2])

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

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")

    # genotype files: chr1.hardcall.prefiltered.mirQTLor.chr1:172411341:C:T.raw
    geno.hardcall.file <- paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".raw")
    geno.dosage.file <- paste0(dir.genotypes, chr, ".dosage.prefiltered.mirQTLor.", esnp, ".raw")

    # output plot file
    output.file <- paste0(dir.output, eqtl, ".pdf")

    # Genotype x Phenotype Plot ########################################################################################

    # subset expression data
    df.expression <- data.frame(expression = assays(vsd[emir,])[["VST"]][1,],
                                expression.corrected = assays(vsd[emir,])[["VST.RESIDUAL"]][1,],
                                DonorID = paste0("D", vsd$donor_id),
                                stringsAsFactors = FALSE)
    df.expression$RNAID <- rownames(df.expression)

    # # join expression data with covariate matrix
    # if (chr == "chrX") {
    #     df.expression %<>%
    #         left_join(cov.mat.x, by = "DonorID")
    # } else {
    #     df.expression %<>%
    #         left_join(cov.mat.auto, by = "DonorID")
    # }

    # # set rownames
    # rownames(df.expression) <- df.expression$DonorID
    # # fit model
    # if (chr == "chrX") {
    #     fit <- lm(formula = as.formula(formula.string.x), data = df.expression)
    # } else {
    #     fit <- lm(formula = as.formula(formula.string.auto), data = df.expression)
    # }
    #
    # # corrected expression from residuals?
    # df.corrected <- data.frame(expression.corrected = resid(fit),
    #                            DonorID = names(resid(fit)),
    #                            stringsAsFactors = FALSE)
    #
    # df.expression %<>%
    #     left_join(df.corrected, by = "DonorID")

    # # genotype labels
    # labels.genotypes <- c(paste(allele.ref, allele.ref, sep = "/"),
    #                                paste(allele.ref, allele.alt.effect, sep = "/"),
    #                                paste(allele.alt.effect, allele.alt.effect, sep = "/"))

    # load genotype data
    df.hardcall.genotypes <- read_delim(geno.hardcall.file,
                                        delim = " ",
                                        col_types = "cccciiii")
    counted.allele.hardcall <- strsplit(colnames(df.hardcall.genotypes)[7], "_")[[1]][2]
    # df.hardcall.genotypes <- read_delim(geno.hardcall.file,
    #                                     delim = " ",
    #                                     skip = 1,
    #                                     col_names = c("DonorID", "DNAID", "PAT", "MAT", "Sex", "PHE", "Genotype_Hard", "Geno_Hard_Het"),
    #                                     col_types = "cccciiii")
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
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # ld for this esnp
    df.ld <- read_table(ld.file)
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
                        baseline=-log10(nom.p.val),
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
    if (eqtl %in% preregress.only) {
        plot.title <- paste0(eqtl,
                             " / BETA: ", signif(df.results.emir$BETA[match(eqtl, df.results.emir$UniName_SNP)], 3),
                             " / P: ", signif(df.results.emir$P[match(eqtl, df.results.emir$UniName_SNP)], 3),
                             " PREREGRESSED COV. ONLY")
    } else if (eqtl %in% both) {
        plot.title <- paste0(eqtl,
                             " / BETA: ", signif(df.results.emir$BETA[match(eqtl, df.results.emir$UniName_SNP)], 3),
                             " / P: ", signif(df.results.emir$P[match(eqtl, df.results.emir$UniName_SNP)], 3),
                             " PREREGRESSED COV. (BOTH)")
    } else {
        stop("eQTL not in list of pregress only or both")
    }

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
    eqtl.WINDOW <- 1e5
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
                           baseline=-log10(nom.p.val),
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
                                     groupAnnotation = "group",
                                     stacking = "dense",
                                     #featureAnnotation="id",
                                     #showFeatureId = TRUE,
                                     fontcolor.item = "black",
                                     shape = "box")

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



