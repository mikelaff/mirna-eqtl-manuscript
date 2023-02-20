# look eqtls and expressed mirnas in relation to spliced introns

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/eqtl_annotation/pdfs/")
dir.create(dir.pdfs, showWarnings = FALSE, recursive = TRUE)

# summarized results of miRNAs and host expression
#df.output.rds <- here("results/rdata_files/2020")

# INPUT ################################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")
# mRNA-eQTLs
mrna.eqtls.dataframe.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTLs.rds")
# sQTLs
sqtls.dataframe.rds <- here("results/external_data/fetal_brain_sQTL/fetal_brain_sQTLs.rds")


# fetal brain mRNA-eQTL co-localizations
mrna.colocs.rds <- here("results/co-localization/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_overlaps_r2at0.8.rds")
# fetal brain sQTL co-localizations
sqtl.colocs.rds <- here("results/co-localization/fetal_brain_sQTL/fetal_brain_sQTL_mirQTL_overlaps_r2at0.8.rds")

# summarized experiment with mirna vst expression values used in this analysis
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# splice ratios by leafcutter
rse.splicing.rds <- here("results/rdata_files/20211217_rse_leafcutter_splicing.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import eQTL Data #####################################################################################################
df.eqtl <- readRDS(eqtls.dataframe.rds)

df.eqtl.gene <- readRDS(mrna.eqtls.dataframe.rds)

df.sqtl.gene <- readRDS(sqtls.dataframe.rds)

# Import Co-localizations ##############################################################################################
df.mrna.colocs <- readRDS(mrna.colocs.rds)

df.sqtl.colocs <- readRDS(sqtl.colocs.rds)

# Import Expression Data ###############################################################################################

# mirna expression
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

# threshold for expression
rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

rm(rse.gene, dds.gene)

# Import Splice Data ###################################################################################################
rse.splicing <- read_rds(rse.splicing.rds)

# subset samples to those in vsd.mirna
rse.splicing <- rse.splicing[,colnames(vsd.mirna)]

# set seqnames for gene expression to UCSC
seqlevelsStyle(rse.splicing) <- "UCSC"
genome(rse.splicing) <- "hg38"
# remove non-standard chroms
seqlevels(rse.splicing, pruning.mode = "coarse") <- CHROMS

# Import miRNA GRanges #################################################################################################

# known and novel mirnas
gr.mirna <- readRDS(gr.mirna.rds)

# expressed primary mirnas
gr.expressed.pri.mirna <-  gr.mirna[ gr.mirna$ID %in% rowRanges(vsd.mirna)$Derives_from ]

# remove non-standard chroms
seqlevels(gr.expressed.pri.mirna, pruning.mode = "coarse") <- CHROMS

# Distance to miRNA ####################################################################################################
df.eqtl$distance.eqtl.to.mirna <- NA

for (i in 1:nrow(df.eqtl)) {
    if (df.eqtl$miR.STRAND[i] == "-") {
        df.eqtl$distance.eqtl.to.mirna[i] <- df.eqtl$miR.END.hg38[i] - df.eqtl$SNP.BP.hg38[i]
    } else {
        df.eqtl$distance.eqtl.to.mirna[i] <- df.eqtl$miR.START.hg38[i] - df.eqtl$SNP.BP.hg38[i]
    }
}

# Find Host Genes ######################################################################################################

# hits.expr.pri.mirna.TO.genes <- findOverlaps(gr.expressed.pri.mirna, vsd.gene)
#
# #length(unique(queryHits(hits.expr.pri.mirna.TO.genes)))
#
# df.expr.pri.mirna.with.host <- as_tibble(gr.expressed.pri.mirna[queryHits(hits.expr.pri.mirna.TO.genes)])
# df.expr.pri.mirna.with.host %<>%
#     dplyr::rename(seqnames.primir = seqnames,
#            start.primir = start,
#            end.primir = end,
#            width.primir = width,
#            strand.primir = strand)
#
# df.host.genes <- as_tibble(rowRanges(vsd.gene)[subjectHits(hits.expr.pri.mirna.TO.genes)])
# df.host.genes %<>%
#     dplyr::rename(seqnames.gene = seqnames,
#            start.gene = start,
#            end.gene = end,
#            width.gene = width,
#            strand.gene = strand)
#
# df.expr.pri.mirna.TO.hosts <- bind_cols(df.expr.pri.mirna.with.host, df.host.genes)
#
# # how many possible host genes per pri-miR
# df.expr.pri.mirna.TO.hosts %>%
#     group_by(UniqueName) %>%
#     summarise(count = as.integer(n())) %>%
#     ggplot(aes(x=count)) +
#     geom_histogram()
#
# df.gene.expr <- tibble(gene.expr = rowMeans( assay(vsd.gene) ), gene_id = names(rowMeans( assay(vsd.gene) )))
#
# df.expr.pri.mirna.TO.hosts <- left_join(df.expr.pri.mirna.TO.hosts, df.gene.expr, by = "gene_id")
#
# df.dups <- df.expr.pri.mirna.TO.hosts[duplicated(df.expr.pri.mirna.TO.hosts$ID) | duplicated(df.expr.pri.mirna.TO.hosts$ID, fromLast = TRUE), ]

# Find Host Introns ####################################################################################################


hits.expr.pri.mirna.TO.introns <- findOverlaps(gr.expressed.pri.mirna, rse.splicing)

#length(unique(queryHits(hits.expr.pri.mirna.TO.genes)))

df.expr.pri.mirna.with.host.intron <- as_tibble(gr.expressed.pri.mirna[queryHits(hits.expr.pri.mirna.TO.introns)])
df.expr.pri.mirna.with.host.intron %<>%
    dplyr::rename(seqnames.primir = seqnames,
                  start.primir = start,
                  end.primir = end,
                  width.primir = width,
                  strand.primir = strand)

df.host.introns <- as_tibble(rowRanges(rse.splicing)[subjectHits(hits.expr.pri.mirna.TO.introns)])
df.host.introns$intronID <- names(rowRanges(rse.splicing)[subjectHits(hits.expr.pri.mirna.TO.introns)])
df.host.introns %<>%
    dplyr::rename(seqnames.intron = seqnames,
                  start.intron = start,
                  end.intron = end,
                  width.intron = width,
                  strand.intron = strand)

df.expr.pri.mirna.TO.introns <- bind_cols(df.expr.pri.mirna.with.host.intron, df.host.introns)

# how many possible host introns per pri-miR
df.expr.pri.mirna.TO.introns %>%
    group_by(UniqueName) %>%
    summarise(count = as.integer(n())) %>%
    ggplot(aes(x=count)) +
    geom_histogram()

#df.gene.expr <- tibble(gene.expr = rowMeans( assay(vsd.gene) ), gene_id = names(rowMeans( assay(vsd.gene) )))

#df.expr.pri.mirna.TO.hosts <- left_join(df.expr.pri.mirna.TO.hosts, df.gene.expr, by = "gene_id")

#df.dups <- df.expr.pri.mirna.TO.introns[duplicated(df.expr.pri.mirna.TO.introns$ID) | duplicated(df.expr.pri.mirna.TO.introns$ID, fromLast = TRUE), ]



# Expression Correlation ###############################################################################################

# loop over host genes and look at correlation to miRNA expression

# build table of mature mirnas with hosts
df.mirna <- as_tibble(rowRanges(vsd.mirna))

df.expr.pri.mirna.TO.introns %<>%
    select(ends_with(".primir"),
           ID.primir = ID,
           Alias.primir = Alias,
           Name.primir = Name,
           sequence.primir = sequence,
           UniqueName.primir = UniqueName,
           ends_with(".intron"),
           intronID)

# for each expression miRNA, combine with primary/host gene information
df.mirna %<>%
    left_join(df.expr.pri.mirna.TO.introns, by = c("Derives_from" = "ID.primir"))

# remove expression miRNAs that don't have hosts
#df.mirna <- df.mirna[!is.na(df.mirna$gene_id),]

# get average expression of the miRNA
#df.mirna.expr <- tibble(expression.mirna = rowMeans( assay(vsd.mirna) ), uniqueName = names(rowMeans( assay(vsd.mirna) )))
#df.mirna %<>%
#    left_join(df.mirna.expr, by = "uniqueName")


# confirm expression matrix has the same columns in miRNA and gene expression
all(colnames(vsd.mirna) == colnames(rse.splicing))

# only mirna with hosts
df.mirna %>%
    filter(! is.na(intronID)) -> df.mirna.with.host.intron

# expression matrix for miRNA and intron: raw
mat.mirna.expr <- assay(vsd.mirna)[df.mirna.with.host.intron$uniqueName, ]
mat.splicing.usage.raw <- assay(rse.splicing, 1)[df.mirna.with.host.intron$intronID, ]

# get correlation
df.mirna.with.host.intron$corr.raw <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.splicing.usage.raw[i,]))

# expression matrix for miRNA and intron: normalized
mat.mirna.expr <- assay(vsd.mirna)[df.mirna.with.host.intron$uniqueName, ]
mat.splicing.usage.norm <- assay(rse.splicing, 2)[df.mirna.with.host.intron$intronID, ]

# get correlation
df.mirna.with.host.intron$corr.norm <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.splicing.usage.norm[i,]))


# do they have eqtls
df.mirna.with.host.intron$have.mirna.eqtl <- df.mirna.with.host.intron$uniqueName %in% df.eqtl$emiR
df.mirna.with.host.intron$have.mirna.eqtl.high <- df.mirna.with.host.intron$uniqueName %in% df.eqtl$emiR[df.eqtl$SIGNIFICANCE == "eigenMT_fdr5percent"]
df.mirna.with.host.intron$have.gene.sqtl <- df.mirna.with.host.intron$intronID %in% df.sqtl.gene$intron

# is there a coloc
df.mirna.with.host.intron$have.mirna.sqtl.coloc <- df.mirna.with.host.intron$uniqueName %in% df.sqtl.colocs$emiR.mirQTL & df.mirna.with.host.intron$intronID %in% df.sqtl.colocs$intron.mQTL

colscale <- scale_fill_manual(values = c("darkblue", "darkorange"))

df.mirna.with.host.intron %>%
    ggplot(aes(x = corr.norm, fill = have.gene.sqtl)) +
    geom_histogram(bins = 11) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mirna_host_corr_eqtl.pdf", height = 2, width = 2.5)

df.mirna.with.host %>%
    ggplot(aes(x = corr, fill = have.mirna.eqtl.high)) +
    geom_histogram(bins = 15) +
    scale_fill_manual(values = c("darkblue", "darkorange"))

df.mirna.with.host %>%
    ggplot(aes(x = corr, fill = have.gene.eqtl)) +
    geom_histogram(bins = 15) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mirna_host_corr_mrna-eqtl.pdf", height = 2, width = 2.5)

df.mirna.with.host %>%
    ggplot(aes(x = corr, fill = have.gene.sqtl)) +
    geom_histogram(bins = 15) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mirna_host_corr_sqtl.pdf", height = 2, width = 2.5)

df.mirna.with.host %>%
    ggplot(aes(x = corr, fill = have.mirna.gene.coloc)) +
    geom_histogram(bins = 15) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mirna_host_corr_mrna-eqtl_coloc.pdf", height = 2, width = 2.5)

df.mirna.with.host %>%
    ggplot(aes(x = corr, fill = have.mirna.sqtl.coloc)) +
    geom_histogram(bins = 15) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mirna_host_corr_sqtl_coloc.pdf", height = 2, width = 2.5)

# get correlation for mrna eqtl colocs
mat.mirna.expr <- assay(vsd.mirna)[df.mrna.colocs$emiR.mirQTL, ]
mat.gene.expr <- assay(vsd.gene)[df.mrna.colocs$ENSG.mQTL, ]

df.mrna.colocs$corr <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.gene.expr[i,]))

# are they in a host
df.mrna.colocs$emir.in.egene <- df.mrna.colocs$emiR.mirQTL %in% df.mirna.with.host$uniqueName & df.mrna.colocs$ENSG.mQTL %in% df.mirna.with.host$gene_id

df.mrna.colocs %>%
    filter(! duplicated(corr)) %>%
    ggplot(aes(x = corr, fill = emir.in.egene)) +
    geom_histogram(bins = 5) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mrna-eqtl_colocs_host_corr.pdf", height = 2, width = 2.5)


# beta direction
df.mrna.colocs$BETA.mQTL.corrected <- NA
for (i in 1:nrow(df.mrna.colocs)) {
    if (df.mrna.colocs$EFFECT.ALLELE.mirQTL[i] == df.mrna.colocs$ALLELE_MAJOR_EFFECT.mQTL[i]) {
        df.mrna.colocs$BETA.mQTL.corrected[i] <- df.mrna.colocs$BETA.mQTL[i]
    } else {
        df.mrna.colocs$BETA.mQTL.corrected[i] <- -df.mrna.colocs$BETA.mQTL[i]
    }
}

df.mrna.colocs %>%
    filter(! duplicated(corr)) %>%
    ggplot(aes(x = BETA.mirQTL, y = BETA.mQTL.corrected, color = emir.in.egene)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(limits = c(-1,3)) +
    scale_x_continuous(limits = c(-1,3)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/mrna-eqtl_colocs_betas.pdf", height = 3, width = 3)


# get correlation for sqtl colocs
mat.mirna.expr <- assay(vsd.mirna)[df.sqtl.colocs$emiR.mirQTL, ]
mat.gene.expr <- assay(vsd.gene)[df.sqtl.colocs$ensemblID.mQTL, ]

df.sqtl.colocs$corr.expression <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.gene.expr[i,]))

mat.splice.raw <- as.matrix(df.splice.raw[, colnames(mat.mirna.expr)])
rownames(mat.splice.raw) <- df.splice.raw$ID
mat.splice.raw <- mat.splice.raw[df.sqtl.colocs$intron.mQTL, ]

mat.splice.norm <- as.matrix(df.splice.norm[, colnames(mat.mirna.expr)])
rownames(mat.splice.norm) <- df.splice.norm$ID
mat.splice.norm <- mat.splice.norm[df.sqtl.colocs$intron.mQTL, ]

all(colnames(mat.mirna.expr) == colnames(mat.splice.raw))
all(colnames(mat.mirna.expr) == colnames(mat.splice.norm))

df.sqtl.colocs$corr.splice.raw <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.splice.raw[i,]))
df.sqtl.colocs$corr.splice.norm <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.splice.norm[i,]))



# are they in a host
df.sqtl.colocs$emir.in.sgene <- df.sqtl.colocs$emiR.mirQTL %in% df.mirna.with.host$uniqueName & df.sqtl.colocs$ensemblID.mQTL %in% df.mirna.with.host$gene_id

df.sqtl.colocs %>%
    filter(! duplicated(corr.expression)) %>%
    ggplot(aes(x = corr, fill = emir.in.sgene)) +
    geom_histogram(bins = 5) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/sqtl_colocs_host_corr.pdf", height = 2, width = 2.5)

# beta direction
df.sqtl.colocs$BETA.mQTL.corrected <- NA
for (i in 1:nrow(df.sqtl.colocs)) {
    if (df.sqtl.colocs$EFFECT.ALLELE.mirQTL[i] == df.sqtl.colocs$A1.mQTL[i]) {
        df.sqtl.colocs$BETA.mQTL.corrected[i] <- df.sqtl.colocs$beta.mQTL[i]
    } else {
        df.sqtl.colocs$BETA.mQTL.corrected[i] <- -df.sqtl.colocs$beta.mQTL[i]
    }
}

df.sqtl.colocs %>%
    filter(! duplicated(corr.expression)) %>%
    ggplot(aes(x = BETA.mirQTL, y = BETA.mQTL.corrected, color = emir.in.sgene)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    scale_y_continuous(limits = c(-2,2)) +
    scale_x_continuous(limits = c(-2,2)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/sqtl_colocs_betas.pdf", height = 3, width = 3)


# distances
df.eqtl$mirna.has.host <- df.eqtl$emiR %in% df.mirna.with.host$uniqueName
df.eqtl$eqtl.has.mrna.coloc <- df.eqtl$eQTL %in% df.mrna.colocs$eQTL.mirQTL
df.eqtl$eqtl.has.sqtl.coloc <- df.eqtl$eQTL %in% df.sqtl.colocs$eQTL.mirQTL


df.eqtl %>%
    ggplot(aes(x = distance.eqtl.to.mirna)) +
    geom_histogram(bins = 11)

df.eqtl %>%
    ggplot(aes(x = distance.eqtl.to.mirna, fill = mirna.has.host)) +
    geom_histogram(bins = 11) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    plotTheme("figure")

#ggsave("~/Desktop/eqtl_distance_hosts.pdf", height = 2, width = 2.5)


df.mirna.with.host$category <- "neither"
df.mirna.with.host$category[df.mirna.with.host$have.mirna.gene.coloc] <- "mRNA-eQTL-coloc"
df.mirna.with.host$category[df.mirna.with.host$have.mirna.sqtl.coloc] <- "sQTL-coloc"

df.mirna.with.host$category <- factor(df.mirna.with.host$category, levels = c("neither", "mRNA-eQTL-coloc", "sQTL-coloc"), ordered = TRUE)

df.mirna.with.host %>%
    filter(have.mirna.eqtl) %>%
    ggplot(aes(x = corr, fill = category)) +
    geom_histogram(bins = 11) +
    scale_fill_manual(values = c("gray", "darkblue", "darkorange")) +
    scale_y_continuous(expand = expansion(mult = c(0,.1))) +
    scale_x_continuous(breaks = c(-0.5,0,0.5)) +
    plotTheme("figure")







rm()





















# # for each miRNA with a host gene, does its expression correlate with the host gene
#
# # mature mirnas with hosts
#
# df.mirna <- data.frame(ID = rowRanges(vsd.mirna)$ID, uniName = rownames(vsd.mirna), stringsAsFactors = FALSE)
#
# df.mirna$primir <- gr.mirna$Derives_from[match(df.mirna$ID, gr.mirna$ID)]
#
# df.mirna <- left_join(df.mirna, df.expr.pri.mirna.TO.hosts, by = c("primir" = "ID"))
#
# df.mirna <- df.mirna[!is.na(df.mirna$gene_id),]
#
# df.mirna.expr <- data.frame(mirna.expr = rowMeans( assay(vsd.mirna) ), uniName = names(rowMeans( assay(vsd.mirna) )), stringsAsFactors = FALSE)
# df.mirna <- left_join(df.mirna, df.mirna.expr, by = "uniName")
#
# all(colnames(vsd.mirna) == colnames(vsd.gene))
#
# mat.mirna.expr <- assay(vsd.mirna)[df.mirna$uniName,]
# mat.gene.expr <- assay(vsd.gene)[df.mirna$gene_id,]
#
# df.mirna$corr <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.gene.expr[i,]))
#
# df.mirna$eqtl <- df.mirna$uniName %in% df.eqtl$emiR
#
# df.mirna %>%
#     ggplot(mapping = aes(x = gene.expr, y = mirna.expr, color = corr)) +
#     geom_point() +
#     #geom_line(aes(group = uniName)) +
#     labs(y = "Mean miRNA Expression (VST)",
#          x = "Mean Host Gene Expression (VST)",
#          title = "miRNA-Host Expression Correlation",
#          color = "Pearson Corr.") +
#     scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
#
# ggsave(paste0(dir.pdfs, "mirna_host_expression_correlation.pdf"))
#
# df.mirna %>%
#     dplyr::filter(eqtl) %>%
#     ggplot(mapping = aes(x = gene.expr, y = mirna.expr, color = corr)) +
#     geom_point() +
#     geom_line(aes(group = uniName)) +
#     labs(y = "Mean emiR Expression (VST)",
#          x = "Mean Host Gene Expression (VST)",
#          title = "emiR-Host Expression Correlation",
#          color = "Pearson Corr.") +
#     scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)
#
# ggsave(paste0(dir.pdfs, "emir_host_expression_correlation.pdf"))
#
# # 196 emiRs: 146 have host gene
#
#
#
#
#
# # Scratch ###############
#
#
# df.coloc <- readRDS(here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalizations.rds"))
# df.coloc$MIR_ENSG <- paste(df.coloc$emiR.mirQTL, df.coloc$ENSG.mQTL, sep = "_")
#
# length(unique(df.coloc$MIR_ENSG))
#
# df.coloc$COLOC <- paste(df.coloc$eQTL.mirQTL, df.coloc$SNP_ENSG.mQTL, sep = "_")
# length(unique(df.coloc$COLOC))
#
# df.host.mir.pairs <- df.mirna
#
# df.host.mir.pairs$MIR_HOST <- paste(df.host.mir.pairs$uniName, df.host.mir.pairs$gene_id, sep = "_")
# sum(duplicated(df.host.mir.pairs$MIR_HOST))
#
# df.coloc.unique.mir.gene.pairs <- df.coloc[!duplicated(df.coloc$MIR_ENSG),]
#
# sum(df.coloc.unique.mir.gene.pairs$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST)
#
#
# df.coloc.unique.mir.gene.pairs$emir_hosted_by_egene <- df.coloc.unique.mir.gene.pairs$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST
#
#
# df.coloc.unique.mir.gene.pairs %<>%
#     left_join(select(df.host.mir.pairs, MIR_HOST, gene.expr, mirna.expr, corr, eqtl), by = c("MIR_ENSG" = "MIR_HOST"))
#
#
# df.coloc.unique.mir.gene.pairs %>%
#     group_by(emir_hosted_by_egene) %>%
#     summarise(mean_corr = mean(emiR_eGene_Corr)) -> df.mean
#
# df.coloc.unique.mir.gene.pairs %>%
#     ggplot(aes(x = emiR_eGene_Corr,
#                fill = emir_hosted_by_egene)) +
#     geom_histogram() +
#     facet_wrap(~emir_hosted_by_egene) +
#     geom_vline(data = df.mean, aes(xintercept = mean_corr)) +
#     geom_text(data = df.mean, aes(x = mean_corr, y = Inf, label = round(mean_corr, digits = 3)), vjust = 1.1, hjust = -0.1, size = 7) +
#     scale_fill_manual(values = cbPalette) +
#     labs(fill = "emiR\nHosted By\neGene",
#          x = "emiR / eGene Expression Correlation") +
#     theme(axis.text = element_text(size = 14),
#           axis.title = element_text(size = 16),
#           legend.text = element_text(size = 14),
#           legend.title = element_text(size = 16))
#
# ggsave("~/Desktop/hosted_emir_egene_coloc.pdf", width = 8, height = 5)
#
#
# df.coloc.unique.esnp <- filter(df.coloc, !duplicated(eSNP.mirQTL))
#
# df.coloc.unique.esnp$emir_hosted_by_egene <- df.coloc.unique.esnp$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST
#
# df.coloc.unique.esnp %>%
#     ggplot(aes(x = BETA.mirQTL, y = BETA.mQTL, color = emir_hosted_by_egene)) +
#     geom_point(size = 3) +
#     scale_color_manual(values = cbPalette) +
#     labs(color = "emiR\nHosted By\neGene",
#          x = "miRNA-eQTL BETA",
#          y = "mRNA-eQTL BETA") +
#     theme(axis.text = element_text(size = 14),
#           axis.title = element_text(size = 16),
#           legend.text = element_text(size = 14),
#           legend.title = element_text(size = 16))
#
# ggsave("~/Desktop/hosted_emir_egene_beta.pdf", width = 8, height = 5)








