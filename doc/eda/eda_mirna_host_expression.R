# look at expression of miRNAs and their host genes

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(psych)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/eda/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# summarized experiment with mirna vst expression values used in this analysis
# using only the QCed samples from the miRNA-eQTL analysis
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import and Normalize #################################################################################################
# mirna expression
vsd.mirna <- readRDS(vsd.mirna.rds)

# gene expression
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd.mirna)]

# threshold for expression
rse.gene <- rse.gene[rowSums(assay(rse.gene) >= 10) >= 10, ]

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

rm(rse.gene, dds.gene)

# Import miRNA GRanges #################################################################################################
# known and novel mirnas
gr.mirna <- readRDS(gr.mirna.rds)

# expressed primary mirnas
gr.expr.pri.mirna <-  gr.mirna[ gr.mirna$ID %in% rowRanges(vsd.mirna)$Derives_from ]

# remove non-standard chroms
seqlevels(gr.expr.pri.mirna, pruning.mode = "coarse") <- CHROMS

# Find Host Genes ######################################################################################################

hits.expr.pri.mirna.TO.genes <- findOverlaps(gr.expr.pri.mirna, vsd.gene)

#length(unique(queryHits(hits.expr.pri.mirna.TO.genes)))

df.expr.pri.mirna.with.host <- as_tibble(gr.expr.pri.mirna[queryHits(hits.expr.pri.mirna.TO.genes)])
df.expr.pri.mirna.with.host %<>%
    dplyr::rename(seqnames.primir = seqnames,
                  start.primir = start,
                  end.primir = end,
                  width.primir = width,
                  strand.primir = strand)

df.host.genes <- as_tibble(rowRanges(vsd.gene)[subjectHits(hits.expr.pri.mirna.TO.genes)])
df.host.genes %<>%
    dplyr::rename(seqnames.gene = seqnames,
                  start.gene = start,
                  end.gene = end,
                  width.gene = width,
                  strand.gene = strand)

df.expr.pri.mirna.TO.hosts <- bind_cols(df.expr.pri.mirna.with.host, df.host.genes)

# how many possible host genes per pri-miR
df.expr.pri.mirna.TO.hosts %>%
    group_by(UniqueName) %>%
    summarise(count = as.integer(n())) %>%
    ggplot(aes(x=count)) +
    geom_histogram()

df.gene.expr <- data.frame(gene.expr = rowMeans( assay(vsd.gene) ), gene_id = names(rowMeans( assay(vsd.gene) )), stringsAsFactors = FALSE)

df.expr.pri.mirna.TO.hosts <- left_join(df.expr.pri.mirna.TO.hosts, df.gene.expr, by = "gene_id")

df.dups <- df.expr.pri.mirna.TO.hosts[duplicated(df.expr.pri.mirna.TO.hosts$ID) | duplicated(df.expr.pri.mirna.TO.hosts$ID, fromLast = TRUE), ]


# Expression Correlation ###############################################################################################

# loop over host genes and look at correlation to miRNA expression

# mature mirnas with hosts
df.mirna <- tibble(ID = rowRanges(vsd.mirna)$ID, uniqueName = rownames(vsd.mirna))

# get primary miRNA ID
df.mirna$Derives_from <- gr.mirna$Derives_from[match(df.mirna$ID, gr.mirna$ID)]

# for each expression miRNA, combine with primary/host gene information
df.mirna <- left_join(df.mirna, df.expr.pri.mirna.TO.hosts, by = c("Derives_from" = "ID"))

# remove expression miRNAs that don't have hosts
df.mirna <- df.mirna[!is.na(df.mirna$gene_id),]

# get average expression of the miRNA
df.mirna.expr <- tibble(mirna.expr = rowMeans( assay(vsd.mirna) ), uniqueName = names(rowMeans( assay(vsd.mirna) )))
df.mirna <- left_join(df.mirna, df.mirna.expr, by = "uniqueName")


# confirm expression matrix has the same columns in miRNA and gene expression
all(colnames(vsd.mirna) == colnames(vsd.gene))

# expression matrix for miRNA and gene expression
mat.mirna.expr <- assay(vsd.mirna)[df.mirna$uniqueName,]
mat.gene.expr <- assay(vsd.gene)[df.mirna$gene_id,]

# get correlation
df.mirna$corr <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.gene.expr[i,]))



# all(colnames(mat.mirna.expr) == colnames(mat.gene.expr))
#
# mat.corr <- cor(x=t(mat.mirna.expr), y=t(mat.gene.expr), method="pearson")
#
# # get correlations and significance values
# corr <- corr.test(x = t(mat.mir),
#                   y = t(mat.gene),
#                   method = "pearson",
#                   adjust = "fdr",
#                   ci=FALSE)
# # extract correlation
# mat.corr.d <- corr$r
# # extract corrected p-values
# mat.pval <- corr$p
#
# df.corr <- data.frame(ensg = colnames(mat.corr.d), corr.pearson = mat.corr.d[1,], corr.pval = mat.pval[1,])
#
#
# df.gene.summary %<>%
#     dplyr::left_join(df.corr, by = c("uniqueName" = "ensg"))

df.mirna$lm_beta <- NA
df.mirna$lm_p <- NA


#pdf(paste0(dir.pdfs, "host_gene_expression_corr.pdf"))

for (i in 1:nrow(df.mirna)) {
#for (i in 1:5) {

    print(i)

    mirna <- df.mirna$uniqueName[i]
    host.gene <- df.mirna$gene_id[i]
    host.gene.symbol <- df.mirna$symbol[i]

    # get mirna and host gene expression
    df.expression <- tibble(mirna_expression = assay(vsd.mirna)[mirna, ], host_gene_expression = assay(vsd.gene)[host.gene, ], rnaid_mirna = colnames(vsd.mirna), rnaid_gene = colnames(vsd.gene))

    stopifnot(all(df.expression$rnaid_mirna == df.expression$rnaid_gene))

    # linear model: miRNA expression ~ host gene expression
    model.fit <- lm(formula = host_gene_expression ~ mirna_expression, data = df.expression)

    lm.intercept <- model.fit$coefficients["(Intercept)"]
    lm.beta <- signif(model.fit$coefficients["mirna_expression"], 3)
    lm.p <- signif(summary(model.fit)$coefficients[,4]["mirna_expression"], 3)

    df.mirna$lm_beta[i] <- lm.beta
    df.mirna$lm_p[i] <- lm.p

    line.color <- "grey"
    if (lm.p < 0.05) {
        if (lm.beta > 0) {
            line.color <- "red"
        }
        if (lm.beta < 0) {
            line.color <- "blue"
        }
    }

    df.expression %>%
        ggplot(aes(x = mirna_expression, y = host_gene_expression)) +
        geom_point(size = 1) +
        #geom_smooth(method = "lm") +
        geom_abline(slope = lm.beta, intercept = lm.intercept, color = line.color, size = 1) +
        labs(y = paste(host.gene.symbol, host.gene, "Expression (vst normalized)"),
             x = paste(mirna, "Expression (vst normalized)"),
             title = "Host Gene Expression Correlation",
             caption = paste("beta:", lm.beta, "p:", lm.p)) -> p


    #print(p)

}

#dev.off()


#df.mirna$eqtl <- df.mirna$uniName %in% df.eqtl$emiR

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

# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_eQTLs_dataFrame.rds")

# mirQTL eqtls
df.eqtls <- readRDS(eqtls.dataframe.rds)

dir.root.hardcalls <- here("results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/")

emir <- "hsa-mir-4707_hsa-miR-4707-3p"
chr <- 14

geno.hardcall.file <- here("results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/chr14/chr14.hardcall.prefiltered.mirQTLor.hsa-mir-4707_hsa-miR-4707-3p.traw")

df.geno <- df.genotype.hardcalls <- read_tsv(geno.hardcall.file, col_types = cols())





