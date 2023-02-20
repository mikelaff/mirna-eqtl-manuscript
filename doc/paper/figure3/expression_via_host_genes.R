# miRNA eQTLs expressed from host gene

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)

# OUTPUT ###############################################################################################################
dir.pdf <- here("doc/paper/figure3/pdf/")
dir.create(dir.pdf, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_eQTLs_dataFrame.rds")

# summarized experiment with mirna vst expression values used in this analysis
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

# Import eQTL Data #####################################################################################################
df.eqtl <- readRDS(eqtls.dataframe.rds)

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

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# normalize expression (VST)
dds.gene <- DESeqDataSet(rse.gene, design = ~1)
vsd.gene <- vst(dds.gene)

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

df.host.genes <- as_tibble(rowRanges(rse.gene)[subjectHits(hits.expr.pri.mirna.TO.genes)])
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


# for each miRNA with a host gene, does its expression correlate with the host gene

# mature mirnas with hosts

df.mirna <- data.frame(ID = rowRanges(vsd.mirna)$ID, uniName = rownames(vsd.mirna), stringsAsFactors = FALSE)

df.mirna$primir <- gr.mirna$Derives_from[match(df.mirna$ID, gr.mirna$ID)]

df.mirna <- left_join(df.mirna, df.expr.pri.mirna.TO.hosts, by = c("primir" = "ID"))

df.mirna <- df.mirna[!is.na(df.mirna$gene_id),]

df.mirna.expr <- data.frame(mirna.expr = rowMeans( assay(vsd.mirna) ), uniName = names(rowMeans( assay(vsd.mirna) )), stringsAsFactors = FALSE)
df.mirna <- left_join(df.mirna, df.mirna.expr, by = "uniName")

all(colnames(vsd.mirna) == colnames(vsd.gene))

mat.mirna.expr <- assay(vsd.mirna)[df.mirna$uniName,]
mat.gene.expr <- assay(vsd.gene)[df.mirna$gene_id,]

df.mirna$corr <- sapply(1:dim(mat.mirna.expr)[1], function(i) cor(mat.mirna.expr[i,], mat.gene.expr[i,]))

df.mirna$eqtl <- df.mirna$uniName %in% df.eqtl$emiR

df.mirna %>%
    ggplot(mapping = aes(x = gene.expr, y = mirna.expr, color = corr)) +
    geom_point() +
    #geom_line(aes(group = uniName)) +
    labs(y = "Mean miRNA Expression (VST)",
         x = "Mean Host Gene Expression (VST)",
         title = "miRNA-Host Expression Correlation",
         color = "Pearson Corr.") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

ggsave(paste0(dir.pdfs, "mirna_host_expression_correlation.pdf"))

df.mirna %>%
    dplyr::filter(eqtl) %>%
    ggplot(mapping = aes(x = gene.expr, y = mirna.expr, color = corr)) +
    geom_point() +
    geom_line(aes(group = uniName)) +
    labs(y = "Mean emiR Expression (VST)",
         x = "Mean Host Gene Expression (VST)",
         title = "emiR-Host Expression Correlation",
         color = "Pearson Corr.") +
    scale_color_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0)

ggsave(paste0(dir.pdfs, "emir_host_expression_correlation.pdf"))

# 196 emiRs: 146 have host gene





# Scratch ###############


df.coloc <- readRDS(here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_colocalizations.rds"))
df.coloc$MIR_ENSG <- paste(df.coloc$emiR.mirQTL, df.coloc$ENSG.mQTL, sep = "_")

length(unique(df.coloc$MIR_ENSG))

df.coloc$COLOC <- paste(df.coloc$eQTL.mirQTL, df.coloc$SNP_ENSG.mQTL, sep = "_")
length(unique(df.coloc$COLOC))

df.host.mir.pairs <- df.mirna

df.host.mir.pairs$MIR_HOST <- paste(df.host.mir.pairs$uniName, df.host.mir.pairs$gene_id, sep = "_")
sum(duplicated(df.host.mir.pairs$MIR_HOST))

df.coloc.unique.mir.gene.pairs <- df.coloc[!duplicated(df.coloc$MIR_ENSG),]

sum(df.coloc.unique.mir.gene.pairs$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST)


df.coloc.unique.mir.gene.pairs$emir_hosted_by_egene <- df.coloc.unique.mir.gene.pairs$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST


df.coloc.unique.mir.gene.pairs %<>%
    left_join(select(df.host.mir.pairs, MIR_HOST, gene.expr, mirna.expr, corr, eqtl), by = c("MIR_ENSG" = "MIR_HOST"))


df.coloc.unique.mir.gene.pairs %>%
    group_by(emir_hosted_by_egene) %>%
    summarise(mean_corr = mean(emiR_eGene_Corr)) -> df.mean

df.coloc.unique.mir.gene.pairs %>%
    ggplot(aes(x = emiR_eGene_Corr,
               fill = emir_hosted_by_egene)) +
    geom_histogram() +
    facet_wrap(~emir_hosted_by_egene) +
    geom_vline(data = df.mean, aes(xintercept = mean_corr)) +
    geom_text(data = df.mean, aes(x = mean_corr, y = Inf, label = round(mean_corr, digits = 3)), vjust = 1.1, hjust = -0.1, size = 7) +
    scale_fill_manual(values = cbPalette) +
    labs(fill = "emiR\nHosted By\neGene",
         x = "emiR / eGene Expression Correlation") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))

ggsave("~/Desktop/hosted_emir_egene_coloc.pdf", width = 8, height = 5)


df.coloc.unique.esnp <- filter(df.coloc, !duplicated(eSNP.mirQTL))

df.coloc.unique.esnp$emir_hosted_by_egene <- df.coloc.unique.esnp$MIR_ENSG %in% df.host.mir.pairs$MIR_HOST

df.coloc.unique.esnp %>%
    ggplot(aes(x = BETA.mirQTL, y = BETA.mQTL, color = emir_hosted_by_egene)) +
    geom_point(size = 3) +
    scale_color_manual(values = cbPalette) +
    labs(color = "emiR\nHosted By\neGene",
         x = "miRNA-eQTL BETA",
         y = "mRNA-eQTL BETA") +
    theme(axis.text = element_text(size = 14),
          axis.title = element_text(size = 16),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 16))

ggsave("~/Desktop/hosted_emir_egene_beta.pdf", width = 8, height = 5)


