# investigate anomalies with miR-101

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

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_plots/pdfs/")
dir.pngs <- here("doc/mirQTL_plots/pngs/")

# output folder for combined plots
dir.output <- here("doc/mirQTL_plots/pdfs/mirQTL_ld_snps/")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20191204_mirQTLor_VST_miRNA_expression/20191204_mirQTLor_VST_miRNA_expression_rse.rds")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/emmax/association_results/20191204_mirQTLor/sample_genotypes/")

# Root directory for prefiltered hardcall genotypes per mirna
dir.root.hardcalls <- here("results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/")
# Root directory for prefiltered dosage genotyeps per mirna
dir.root.dosages <- here("results/emmax/tfiles/prefiltered_mirQTLor_dosages/")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20191204_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Covariate matrix autosomes
cov.mat.auto.file <- here("results/emmax/covariates/20191204_chrAutosomes.mirQTLor.columnLabeled.cov")
# Covariate matrix x chromosome
cov.mat.x.file <- here("results/emmax/covariates/20191204_chrX.mirQTLor.columnLabeled.cov")

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
    mutate(eqtl = paste(emir, esnp, sep = "+"))
#filter for only hsa-mir-101-1 and hsa-mir-101-2 results
#df.eqtls %<>%
#    dplyr::filter(!grepl("hsa-mir-101-", emir))

# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))
#filter for only hsa-mir-101-1 and hsa-mir-101-2 results
#df.results %<>%
#    dplyr::filter(!grepl("hsa-mir-101-", UniName))

# expression data
vsd <- readRDS(vsd.rds)

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# autosome covariates file
cov.mat.auto <- read_delim(cov.mat.auto.file, delim = " ")
cov.mat.auto %<>%
    dplyr::rename(DonorID = X1, DNAID = X2)

# covariates
batchVars.auto <- c(paste0("PC", 1:10, ".genotype"),
                    paste0("PC", 1:10, ".expression"),
                    "PoolPool2", "PoolPool3", "PoolPool4",
                    "PoolPool5", "PoolPool6", "PoolPool7",
                    "PoolPool8","PurificationMethodmiRNeasy",
                    "PurificationMethodmiRNeasy_mini", "SexM", "RIN", "GestationWeek")

# formula strings
formula.string.auto <- paste("expression ~", paste(batchVars.auto, collapse = " + "))

genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# Investigation Plots ##################################################################################################
# investigate variant associations at: hsa-mir-101-1_hsa-miR-101-3p+chr1:64125648:G:A

#eqtl.ind <- which(df.eqtls$eqtl == "hsa-mir-101-1_hsa-miR-101-3p+chr1:64125648:G:A")

trouble.ind.eqtls <- NULL

for (k in 1:nrow(df.eqtls)) {

    printMessage(paste("Plotting", k, "of", nrow(df.eqtls), ":", df.eqtls$eqtl[k]))

    eqtl.ind <- k

    esnp <- df.eqtls$esnp[eqtl.ind]
    emir <- df.eqtls$emir[eqtl.ind]
    eqtl <- df.eqtls$eqtl[eqtl.ind]

    chr <- strsplit(esnp, ":")[[1]][1]
    esnp.bp <- as.integer(strsplit(esnp, ":")[[1]][2])

    # genotype hardcall file: chr1.hardcall.prefiltered.mirQTLor.hsa-mir-101-1_hsa-miR-101-3p.traw
    geno.hardcall.file <- paste0(dir.root.hardcalls, chr, "/", chr, ".hardcall.prefiltered.mirQTLor.", emir, ".traw")
    df.genotype.hardcalls <- read_tsv(geno.hardcall.file, col_types = cols())

    # genotype dosage file: chr1.dosage.prefiltered.mirQTLor.hsa-mir-101-1_hsa-miR-101-3p.traw
    geno.dosage.file <- paste0(dir.root.dosages, chr, "/", chr, ".dosage.prefiltered.mirQTLor.", emir, ".traw")
    df.genotype.dosages <- read_tsv(geno.dosage.file, col_types = cols())

    # output plot file
    output.file <- paste0(dir.output, eqtl, ".pdf")

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # max and min base pair which can be plotted
    min.bp <- min(df.results.emir$BP.hg38)
    max.bp <- max(df.results.emir$BP.hg38)

    # subset expression data
    df.expression <- data.frame(expression = assay(vsd[emir,])[1,],
                                DonorID = paste0("D", vsd$donor_id),
                                stringsAsFactors = FALSE)
    df.expression$RNAID <- rownames(df.expression)

    # join expression data with covariate matrix
    df.expression %<>%
        left_join(cov.mat.auto, by = "DonorID")

    # set rownames
    rownames(df.expression) <- df.expression$DonorID
    # fit model
    fit <- lm(formula = as.formula(formula.string.auto), data = df.expression)

    # corrected expression from residuals?
    df.corrected <- data.frame(expression.corrected = resid(fit),
                               DonorID = names(resid(fit)),
                               stringsAsFactors = FALSE)

    df.expression %<>%
        left_join(df.corrected, by = "DonorID") %>%
        mutate(DonorID_DNAID = paste(DonorID, DNAID, sep = "_"))

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

    # loop over variants in high LD with index (red variants)
    df.results.emir %>%
        dplyr::filter(!is.na(R2)) %>%
        arrange(desc(R2)) -> df.vars

    if (nrow(df.vars) > 30) {
        df.vars <- df.vars[1:30,]
    }

    pdf(output.file, width = 12, height = 6)

    for (i in 1:nrow(df.vars)) {

        # this variant
        this.var <- df.vars$SNP[i]

        printMessage(paste("Plotting", i, "of", nrow(df.vars), ":", this.var),
                     fillChar = "-")



        # reference allele
        allele.ref <- df.vars$REF[match(this.var, df.vars$SNP)]
        # alt allele (effect allele)
        allele.alt <- df.vars$ALT[match(this.var, df.vars$SNP)]
        # minor allele
        allele.minor <- df.vars$A1[match(this.var, df.vars$SNP)]
        # major allele
        allele.major <- df.vars$A2[match(this.var, df.vars$SNP)]
        # effect allele (for dosage based associations)
        effect.allele<- df.vars$EFFECT.ALLELE[match(this.var, df.vars$SNP)]

        if (effect.allele != allele.ref) {
            stop("Effect allele != reference allele!")
        }


        # hardcall genotype counted allele
        counted.hardcall.allele <- df.genotype.hardcalls$COUNTED[which(df.genotype.hardcalls$SNP == this.var)]
        # genotypes for this variant
        df.genotype.hardcalls %>%
            dplyr::filter(SNP == this.var) %>%
            dplyr::select(starts_with("D")) %>%
            pivot_longer(cols = everything(),
                         names_to = "DonorID_DNAID",
                         values_to = "Genotype_Hard") -> df.hardcalls

        # hardcall genotype labels
        if (counted.hardcall.allele == effect.allele) {
            hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                    paste(allele.alt, allele.ref, sep = "/"),
                                    paste(allele.ref, allele.ref, sep = "/"))
            df.hardcalls$Genotype_Hard <- factor(df.hardcalls$Genotype_Hard, levels = c(0, 1, 2), ordered = TRUE, labels = hardcall.genotypes)
        } else if (counted.hardcall.allele == allele.alt) {
            hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                    paste(allele.alt, allele.ref, sep = "/"),
                                    paste(allele.ref, allele.ref, sep = "/"))
            # format genotype hardcalls as factor
            df.hardcalls$Genotype_Hard <- factor(df.hardcalls$Genotype_Hard, levels = rev(c(0, 1, 2)), ordered = TRUE, labels = hardcall.genotypes)
        } else {
            stop("Something is wrong with genotype labels!")
        }

        # dosage genotype counted allele
        counted.dosage.allele <- df.genotype.dosages$COUNTED[which(df.genotype.dosages$SNP == this.var)]
        # genotypes for this variant
        df.genotype.dosages %>%
            dplyr::filter(SNP == this.var) %>%
            dplyr::select(starts_with("D")) %>%
            pivot_longer(cols = everything(),
                         names_to = "DonorID_DNAID",
                         values_to = "Genotype_Dose") -> df.dosages

        # dosage genotype check
        if (counted.dosage.allele == effect.allele) {
            dosage.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                  paste(allele.alt, allele.ref, sep = "/"),
                                  paste(allele.ref, allele.ref, sep = "/"))
        } else {
            stop("Something is wrong with genotype labels!")
        }

        # combine into single genotype dataframe
        df.genotypes <- left_join(df.dosages, df.hardcalls, by = "DonorID_DNAID")

        rm(df.dosages, df.hardcalls)

        # join genotypes with expression data
        df.expression %>%
            left_join(df.genotypes, by = "DonorID_DNAID") -> df.expression.this.var

        # linear model: corrected_expression ~ genotype_dose
        model.fit <- lm(formula = expression.corrected ~ Genotype_Dose, data = df.expression.this.var)

        lm.intercept <- model.fit$coefficients["(Intercept)"]
        lm.beta <- signif(model.fit$coefficients["Genotype_Dose"], 3)
        lm.p <- signif(summary(model.fit)$coefficients[,4]["Genotype_Dose"], 3)

        #model.fit2 <- lm(formula = as.formula(paste(formula.string.auto, "Genotype_Dose", sep = " + ")), data = df.expression)

        # EMMAX beta and p value
        #emmax.beta <- signif(df.vars$BETA[i], 3)
        #emmax.p <- signif(df.vars$P[i], 3)

        # Geno x Pheno Plots
        # plot geno x pheno
        df.expression.this.var %>%
            ggplot(aes(x = Genotype_Hard, y = expression)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.1, size = 1) +
            #geom_jitter(mapping = aes(color = Sex), width = 0.1, size = 1) +
            scale_x_discrete(drop = FALSE) +
            labs(y = "VST Expression",
                 x = paste0("REF(effect):", allele.ref, " ALT:", allele.alt)) -> p1

        df.expression.this.var %>%
            ggplot(aes(x = Genotype_Dose, y = expression)) +
            geom_point(size = 1) +
            scale_x_continuous(limits = c(0,2)) +
            labs(y = "VST Expression",
                 x = "REF Dosage") -> p2

        df.expression.this.var %>%
            ggplot(aes(x = Genotype_Hard, y = expression.corrected)) +
            geom_boxplot(outlier.shape = NA) +
            geom_jitter(width = 0.1, size = 1) +
            #geom_jitter(mapping = aes(color = Sex), width = 0.1, size = 1) +
            scale_x_discrete(drop = FALSE) +
            labs(y = "VST Expression (residual)",
                 x = paste0("REF(effect):", allele.ref, " ALT:", allele.alt)) -> p3

        df.expression.this.var %>%
            ggplot(aes(x = Genotype_Dose, y = expression.corrected)) +
            geom_point(size = 1) +
            #geom_smooth(method = "lm") +
            geom_abline(slope = lm.beta, intercept = lm.intercept, color = df.vars$color[i], size = 1) +
            scale_x_continuous(limits = c(0,2)) +
            labs(y = "VST Expression (residual)",
                 x = paste("BETA:", lm.beta, "P:", lm.p)) -> p4

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
        #cat("Adding data track\n")
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
        #cat("Adding ideogram track\n")

        itrack <- IdeogramTrack(genome = "hg38", chromosome = chr, bands = ideo)

        # Genome axis track
        #cat("Adding genome axis track\n")

        gtrack <- GenomeAxisTrack(exponent = 0)

        # Gene region track
        #cat("Adding gene region track\n")

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
        #cat("Adding miRNA track\n")

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
                                         chromosome = chr,
                                         col = "purple")

        # plot title
        plot.title <- paste0(df.vars$UniName_SNP[i],
                             " / EMMAX BETA: ", signif(df.vars$BETA[i], 3),
                             " / EMMAX P: ", signif(df.vars$P[i], 3),
                             " / R2 To Index: ", signif(df.vars$R2[i], 2))

        # new multi plot grid
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
                trouble.ind.eqtls <<- c(trouble.ind.eqtls, k)
            }
        )

        # Dot Plot: eSNP ###################################################################################################

        # window: eSNP to emiR + 1e4
        #emir.bp <- start(rowRanges(vsd)[emir])
        buffer.bp <- 5e4

        # always plot the same window, make sure all vars are in the window
        startRange.esnp <- min(df.vars$BP.hg38) - buffer.bp
        endRange.esnp <- max(df.vars$BP.hg38) + buffer.bp

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
        #cat("Adding gene region track\n")

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
        #cat("Adding miRNA track\n")

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

        # highlight track
        highlightTrack.snp <- HighlightTrack(trackList = list(gtrack, overlayTrack, grtrack.esnp, miRTrack.esnp),
                                             start = df.vars$BP.hg38[i] - 5,
                                             end = df.vars$BP.hg38[i] + 5,
                                             chromosome = chr,
                                             col = df.vars$color[i])

        pushViewport(viewport(layout.pos.row = 3, layout.pos.col = 1:3))

        tryCatch(
            {
        plotTracks(list(highlightTrack.snp),
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
                trouble.ind.eqtls <<- c(trouble.ind.eqtls, k)
            }
        )


        print(p1, vp = viewport(layout.pos.row = 2, layout.pos.col = 4))
        print(p2, vp = viewport(layout.pos.row = 3, layout.pos.col = 4))

        print(p3, vp = viewport(layout.pos.row = 2, layout.pos.col = 5))
        print(p4, vp = viewport(layout.pos.row = 3, layout.pos.col = 5))

    }

    dev.off()

}

write_lines(trouble.ind.eqtls, "trouble-ind-eqtls.txt")


# scratch #########

# stop()
#
# library(pheatmap)
#
# df.mir101.expr <- data.frame(mir_101_3p = assay(vsd)["hsa-mir-101-1_hsa-miR-101-3p",], rnaid = names(assay(vsd)["hsa-mir-101-1_hsa-miR-101-3p",]), stringsAsFactors = FALSE)
#
# df.mir101.expr <- dplyr::select(left_join(df.mir101.expr, as.data.frame(colData(vsd)), by = "rnaid"), mir_101_3p, DonorID = donor_id)
# df.mir101.expr$DonorID <- paste0("D", df.mir101.expr$DonorID)
#
# df.cov <- left_join(cov.mat.auto, df.mir101.expr, by = "DonorID")
#
# mat.cov <- as.matrix(df.cov[,4:36])
#
# cormat <- cor(mat.cov, use = "pair")
#
# pdf(paste0(dir.pdfs, "cov_mat_to_mir101_expr_heatmap.pdf"))
# pheatmap(cormat,
#          cluster_rows = FALSE,
#          cluster_cols = FALSE)
# dev.off()


