
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(AnnotationHub)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(psych)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mir4707_trans/pdfs/")

dir.output <- here("doc/mir4707_trans/pdfs/top_prospects/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# summary of trans results by gene
gene.summary.rds <- here("results/emmax_transTotalRNA/association_results/20201019_transTotalRNA_mir4707/compiled/20201019_transTotalRNA_mir4707_genes_dataFrame.rds")

# target predictions by miRNAtap databases
predictions.2sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_2sources.csv")
predictions.3sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_3sources.csv")
predictions.4sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_4sources.csv")

# miRDB 2019 predictions, known mirnas
predictions.mirdb.txt.gz <- here("data/target_predictions/miRDB_v6.0_prediction_result.txt.gz")

# miRWalk 4707-3p predictions
predictions.miRWalk.4707.csv <- here("data/target_predictions/miRWalk_miR-4707-3p_Targets.csv")

# TargetScan7.2 4707-3p predictions
predictions.targetScan.4707.txt <- here("data/target_predictions/TargetScan7.2_miR-4707-3p.predicted_targets.txt")

# miRDB 4707-3p ALT sequence predictions
predictions.mirdb.4707alt.txt <- here("data/target_predictions/miRDB_v6.0_miR-4707-3p_ALT.txt")

# dir with select .ps files for plotting
dir.select.ps <- here("results/emmax_transTotalRNA/association_results/20201019_transTotalRNA_mir4707/compiled/select_raw_ps/")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# cytoband file for ideogram track
file.cyto <- here("data/cytoBandIdeo.txt.gz")

# samples file
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# ranged summarized experiment for expression values
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# ranged summarized experiment for expression values
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")


# GLOBALS ##############################################################################################################
FDR.THRESHOLD <- 0.05 / (1 * 14225)

# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Import Summary Data ##################################################################################################

df.gene.summary <- readRDS(gene.summary.rds)
df.gene.summary <- as_tibble(df.gene.summary)

df.gene.summary %<>%
    filter(!is.na(smallest_P))

# Import Target Predictions ############################################################################################
# miRDB predictions
df.mirdb.predictions <- read_tsv(predictions.mirdb.txt.gz, col_names = c("Name", "REFSEQ", "TargetScore", "ENTREZID"))
df.mirdb.predictions %<>%
    filter(Name == "hsa-miR-4707-3p")

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.predictions$REFSEQ),
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

# combine with predictions df
df.mirdb.predictions %<>%
    dplyr::left_join(rowdata, by = "REFSEQ")

df.mirdb.predictions %<>%
    dplyr::filter(!is.na(ENSEMBL)) %>%
    dplyr::filter(!duplicated(ENSEMBL))

# TargetScan Predictions
df.targetScan.predictions <- read_tsv(predictions.targetScan.4707.txt, na = c("N/A"))

rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.targetScan.predictions$`Target gene`),
                                 keytype = "SYMBOL",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

rowdata %<>%
    dplyr::filter(!is.na(ENSEMBL))

df.targetScan.predictions %<>%
    dplyr::left_join(rowdata, by = c("Target gene" = "SYMBOL"))

df.targetScan.predictions %<>%
    dplyr::filter(!is.na(ENSEMBL)) %>%
    dplyr::filter(!duplicated(ENSEMBL))

# miRWalk Predictions
df.miRWalk.predictions <- read_csv(predictions.miRWalk.4707.csv, guess_max = 5000)

# filter predictions for binding probability of 1?
df.miRWalk.predictions %<>%
    dplyr::filter(bindingp == 1)

rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.miRWalk.predictions$refseqid),
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

df.miRWalk.predictions %<>%
    dplyr::left_join(rowdata, by = c("refseqid" = "REFSEQ"))

df.miRWalk.predictions %<>%
    dplyr::filter(!is.na(ENSEMBL)) %>%
    dplyr::filter(!duplicated(ENSEMBL))

# miRDB 4707 ALT predictions
df.mirdb.4707alt.predictions <- read_tsv(predictions.mirdb.4707alt.txt)

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.4707alt.predictions$`Gene Symbol`),
                                 keytype = "SYMBOL",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

# combine with predictions df
df.mirdb.4707alt.predictions %<>%
    dplyr::left_join(rowdata, by = c("Gene Symbol" = "SYMBOL"))

df.mirdb.4707alt.predictions %<>%
    dplyr::filter(!is.na(ENSEMBL)) %>%
    dplyr::filter(!duplicated(ENSEMBL))


# combine with gene summary
df.gene.summary %<>%
    mutate(predicted_mir4707_target_miRDB = uniqueName %in% df.mirdb.predictions$ENSEMBL,
           predicted_mir4707_target_targetScan = uniqueName %in% df.targetScan.predictions$ENSEMBL,
           predicted_mir4707_target_miRWalk = uniqueName %in% df.miRWalk.predictions$ENSEMBL,
           predicted_mir4707_target_miRDB_ALT = uniqueName %in% df.mirdb.4707alt.predictions$ENSEMBL)

df.gene.summary %<>%
    mutate(predicted_mir4707_target_miRDB_OR_ALT = predicted_mir4707_target_miRDB | predicted_mir4707_target_miRDB_ALT,
           predicted_mir4707_target_OR_ALL = predicted_mir4707_target_miRDB | predicted_mir4707_target_targetScan | predicted_mir4707_target_miRWalk | predicted_mir4707_target_miRDB_ALT)

df.gene.summary %<>%
    rowwise() %>%
    mutate(predicted_mir4707_target_num_sources = sum(predicted_mir4707_target_miRDB + predicted_mir4707_target_targetScan + predicted_mir4707_target_miRWalk + predicted_mir4707_target_miRDB_ALT))

# df.gene.summary %>%
#     ggplot(aes(x = -log10(smallest_P), color = predicted_mir4707_target)) +
#     geom_density() +
#     geom_vline(xintercept = -log10(FDR.THRESHOLD)) +
#     labs(x = "-log10(P)",
#          title = "trans mRNA-eQTLs near miR-4707") +
#     scale_color_manual(values = c("navy", "orange"))
#
# #ggsave("~/Desktop/trans_mRNA-eQTLs_near_4707_P.pdf", height = 5, width = 7)
#
# df.gene.summary %>%
#     ggplot(aes(x = abs(beta_at_smallest_P), color = predicted_mir4707_target)) +
#     geom_density() +
#     labs(x = "abs(BETA)",
#          title = "trans mRNA-eQTLs near miR-4707") +
#     scale_color_manual(values = c("navy", "orange"))
#
# #ggsave("~/Desktop/trans_mRNA-eQTLs_near_4707_beta.pdf", height = 5, width = 7)
#
# df.gene.summary %>%
#     ggplot(aes(x = abs(beta_at_smallest_P), color = predicted_mir4707_target)) +
#     geom_density() +
#     facet_wrap(~sig_P) +
#     labs(x = "abs(BETA)",
#          title = "trans mRNA-eQTLs near miR-4707") +
#     scale_color_manual(values = c("navy", "orange"))
#
# #ggsave("~/Desktop/trans_mRNA-eQTLs_near_4707_beta_by_significant.pdf", height = 5, width = 7)
#
# pdf(paste0(dir.pdfs, "trans_mRNA_eQTLs_at_miR-4707-3p_eQTL_by_Target_Predictions.pdf"))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_miRDB)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "miRDB target predictions (35 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_targetScan)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "TargetScan target predictions (1250 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_miRWalk)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "miRWalk target predictions (3077 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_miRDB_ALT)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "miRDB 4707-3p ALT sequence target predictions (290 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_miRDB_OR_ALT)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "miRDB 4707-3p OR ALT sequence target predictions (324 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = predicted_mir4707_target_OR_ALL)) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "ANY target prediction (4138 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = factor(predicted_mir4707_target_num_sources))) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = "Target sources 0 (10,087), 1 (3649), 2 (464), 3 (25)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("grey40", "navy", "orange", "red")) +
#     scale_x_continuous(limits = c(-0.05, 0.05))
#
# df.gene.summary %>%
#     ggplot(aes(x = beta_at_eqtl, color = factor(predicted_mir4707_target_num_sources))) +
#     geom_density() +
#     labs(x = "Beta",
#          title = "trans mRNA-eQTLs at miR-4707-3p eQTL",
#          caption = ,
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("grey40", "navy", "orange", "red")) +
#     scale_x_continuous(limits = c(-0.05, 0.05)) +
#     facet_wrap(~predicted_mir4707_target_miRDB_ALT)
#
# dev.off()
#

# Expression correlation ##########


# Import Samples
samples <- read_delim(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# Import RSE
rse <- readRDS(rse.gene.rds)

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 100) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd.gene <- vst(dds)

# subset samples
vsd.gene <- vsd.gene[,samples$RNAID]

rm(rse, dds)

rse <- readRDS(rse.mirna.rds)

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd.mirna <- varianceStabilizingTransformation(dds)

# subset samples
vsd.mirna <- vsd.mirna[,samples$RNAID]

rm(rse, dds)


mat.mir <- assay(vsd.mirna["hsa-miR-4707-3p", ])

mat.gene <- assay(vsd.gene)

all(colnames(mat.mir) == colnames(mat.gene))

mat.corr <- cor(x=t(mat.mir), y=t(mat.gene), method="pearson")

# get correlations and significance values
corr <- corr.test(x = t(mat.mir),
                  y = t(mat.gene),
                  method = "pearson",
                  adjust = "fdr",
                  ci=FALSE)
# extract correlation
mat.corr.d <- corr$r
# extract corrected p-values
mat.pval <- corr$p

df.corr <- data.frame(ensg = colnames(mat.corr.d), corr.pearson = mat.corr.d[1,], corr.pval = mat.pval[1,])


df.gene.summary %<>%
    dplyr::left_join(df.corr, by = c("uniqueName" = "ensg"))

# pdf(paste0(dir.pdfs, "expression_corr_by_Target_Predictions.pdf"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_miRDB)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "miRDB target predictions (35 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_targetScan)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "TargetScan target predictions (1250 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_miRWalk)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "miRWalk target predictions (3077 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_miRDB_ALT)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "miRDB 4707-3p ALT sequence target predictions (290 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_miRDB_OR_ALT)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "miRDB 4707-3p OR ALT sequence target predictions (324 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = predicted_mir4707_target_OR_ALL)) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "ANY target prediction (4138 targets)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("navy", "orange"))
#
# df.gene.summary %>%
#     ggplot(aes(x = corr.pearson, color = factor(predicted_mir4707_target_num_sources))) +
#     geom_density() +
#     labs(x = "Pearson Correlation",
#          title = "Expression Correlation (Pearson)",
#          caption = "Target sources 0 (10,087), 1 (3649), 2 (464), 3 (25)",
#          color = "Predicted\nTarget") +
#     scale_color_manual(values = c("grey40", "navy", "orange", "red"))
#
# dev.off()


# Plot specific trans associations #########################

# genes to plot
# all 3 source genes
df.gene.summary %>%
    dplyr::filter(predicted_mir4707_target_num_sources == 3) -> df.toPlot

# top 20 by smallest P
df.gene.summary %>%
    as.data.frame() %>%
    dplyr::top_n(n = 20, wt = P_at_eqtl) %>%
    as_tibble() %>%
    bind_rows(df.toPlot) -> df.toPlot

# top 20 with lowest correlation
df.gene.summary %>%
    as.data.frame() %>%
    dplyr::top_n(n = 20, wt = -corr.pearson) %>%
    as_tibble() %>%
    bind_rows(df.toPlot) -> df.toPlot

# top 20 with lowest beta
df.gene.summary %>%
    as.data.frame() %>%
    dplyr::top_n(n = 20, wt = -beta_at_eqtl) %>%
    as_tibble() %>%
    bind_rows(df.toPlot) -> df.toPlot

# top 20 with highest beta
df.gene.summary %>%
    as.data.frame() %>%
    dplyr::top_n(n = 20, wt = beta_at_eqtl) %>%
    as_tibble() %>%
    bind_rows(df.toPlot) -> df.toPlot


df.names <- AnnotationDbi::select(orgdb,
                                  keys = df.toPlot$uniqueName,
                                  keytype = "ENSEMBL",
                                  columns = c("SYMBOL"))
df.names %>%
    dplyr::distinct() -> df.names

df.toPlot %<>%
    dplyr::left_join(df.names, by = c("uniqueName" = "ENSEMBL"))

# files to plot
#ps.files <- list.files(dir.select.ps)

#dir.output <- paste0(dir.pdfs, "select_raw_ps/")

# ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
ld.file <- paste0(dir.ld, "chr14.hardcall.prefiltered.mirQTLor.chr14:22953244:A:G.ld")
# ld for this esnp
df.ld <- read_table(ld.file, col_types = cols())
df.ld %<>%
    dplyr::select(SNP = SNP_B, R2)

esnp <- "chr14:22953244:A:G"
pos <- 22953244
chr <- "chr14"

# biomart genemart
genemart <- useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")

# ideogram file
ideo <- read.table(file.cyto, header = FALSE, sep = "\t")
colnames(ideo) <- c("chrom", "chromStart", "chromEnd", "name", "gieStain")

# mirna
gr.mirna <- readRDS(gr.mirna.rds)
gr.mirna.primary <- gr.mirna[gr.mirna$type == "miRNA_primary_transcript" | gr.mirna$type == "miRNA_putative_precursor"]

# loop over files for plotting
for (i in 1:nrow(df.toPlot)) {

    chr.gene <- df.toPlot$chr[i]
    ensg <- df.toPlot$uniqueName[i]

    file <- paste(chr.gene, "dosage", "prefiltered", "mir4707", ensg, "ps", sep = ".")

    dir.file <- paste0(here("results/emmax_transTotalRNA/association_results/20201019_transTotalRNA_mir4707/"), chr.gene, "/ps/")

    print(file)

    gene <- as.character(df.toPlot$SYMBOL[i])

    # ps file to plot
    ps.file <- paste0(dir.file, file)
    # output pdf file
    pdf.output <- paste0(dir.output, paste(ensg, gene, "trans.pdf", sep = "_"))

    # read results file
    df.result <- read_tsv(ps.file)
    df.result %<>%
        dplyr::left_join(df.ld, by = "SNP")
    df.result$CHR <- chr
    df.result$BP.hg38 <- as.numeric(sapply(strsplit(df.result$SNP, ":"), `[`, 2))

    # dot color based on LD
    df.result$color <- "navy"
    df.result$color[which(df.result$R2 >= 0.8)] <- "red"
    df.result$color[which(df.result$R2 >= 0.6 & df.result$R2 < 0.8)] <- "orange"
    df.result$color[which(df.result$R2 >= 0.4 & df.result$R2 < 0.6)] <- "green"
    df.result$color[which(df.result$R2 >= 0.2 & df.result$R2 < 0.4)] <- "lightblue"
    df.result$color[which(df.result$SNP == esnp)] <- "purple"

    # max and min base pair which can be plotted
    min.bp <- min(df.result$BP.hg38)
    max.bp <- max(df.result$BP.hg38)

    # plot window
    buffer.bp <- 1e5

    # always plot the same window, make sure all vars are in the window
    startRange.esnp <- pos - buffer.bp
    endRange.esnp <- pos + buffer.bp

    # granges for mirQTL
    gr.result <- GRanges(seqnames = df.result$CHR,
                         ranges = IRanges(start = df.result$BP.hg38,
                                          end = df.result$BP.hg38,
                                          names = df.result$SNP),
                         neglog10p = -log10(df.result$P))

    # purple plot
    purpleTrack <- DataTrack(gr.result[esnp],
                             name = "miRNA -log10(P)",
                             type = "p",
                             legend=FALSE,
                             ylim=c(0,max(gr.result$neglog10p)),
                             col="purple",
                             cex=1.25,
                             pch=18)

    # red plot
    if (sum(df.result$color == "red") == 0) {
        redTrack <- DataTrack()
    } else {
        redTrack <- DataTrack(gr.result[df.result$color == "red"],
                              name = "miRNA -log10(P)",
                              type = "p",
                              legend=FALSE,
                              ylim=c(0,max(gr.result$neglog10p)),
                              col="red")
    }

    # orange plot
    if (sum(df.result$color == "orange") == 0) {
        orangeTrack <- DataTrack()
    } else {
        orangeTrack <- DataTrack(gr.result[df.result$color == "orange"],
                                 name = "miRNA -log10(P)",
                                 type = "p",
                                 legend=FALSE,
                                 ylim=c(0,max(gr.result$neglog10p)),
                                 col="orange")
    }

    # green plot
    if (sum(df.result$color == "green") == 0) {
        greenTrack <- DataTrack()
    } else {
        greenTrack <- DataTrack(gr.result[df.result$color == "green"],
                                name = "miRNA -log10(P)",
                                type = "p",
                                legend=FALSE,
                                ylim=c(0,max(gr.result$neglog10p)),
                                col="green")
    }

    # lightblue plot
    if (sum(df.result$color == "lightblue") == 0) {
        lightblueTrack <- DataTrack()
    } else {
        lightblueTrack <- DataTrack(gr.result[df.result$color == "lightblue"],
                                    name = "miRNA -log10(P)",
                                    type = "p",
                                    legend=FALSE,
                                    ylim=c(0,max(gr.result$neglog10p)),
                                    col="lightblue")
    }

    # navy plot
    navyTrack <- DataTrack(gr.result[df.result$color == "navy"],
                           name = "miRNA -log10(P)",
                           type = "p",
                           baseline=-log10(FDR.THRESHOLD),
                           col.baseline="black",
                           lty.baseline=2,
                           lwd.baseline=1,
                           legend=FALSE,
                           ylim=c(0,max(gr.result$neglog10p)),
                           col="navy")

    overlayTrack <- OverlayTrack(trackList = list(navyTrack,
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

    pdf(pdf.output, width = 12, height = 4)

    plotTracks(list(itrack, gtrack, grtrack.esnp, miRTrack.esnp, overlayTrack),
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
               main = paste("mir4707_trans", ensg, gene, sep = " - "),
               cex.main=1)

    dev.off()









}







