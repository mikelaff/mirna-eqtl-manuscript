# classification of expressed miRNAs

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(Gviz)
library(gridExtra)
library(biomaRt)
library(DESeq2)
library(GenomicRanges)
library(RColorBrewer)
library(readxl)

library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(GenomicFeatures)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(annotatr)

# intergenic

# intragenic:
# intronic
# exonic
# antisense
# junction
# readthrough < 4kb?
# divergent < 2kb?

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/eqtl_annotation/pdfs/")

# INPUT ################################################################################################################
# eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold
nom.p.val.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# summarized experiment with gene expression data
rse.gene.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# Huan et. al. 2015, eQTLs (not corrected for cell type?)
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


# Import #########################################################################################################

vsd <- readRDS(vsd.rds)

# subset for only miRBase miRNAs
vsd.mirbase <- vsd[rowData(vsd)$source == "miRBase_v22",]
# granges of mirbase expressed mature mirnas
gr.mirbase.expressed <- rowRanges(vsd.mirbase)

# import mirna annotations
gr.mirna <- readRDS(gr.mirna.rds)
# remove non-standard chroms
seqlevels(gr.mirna, pruning.mode = "coarse") <- CHROMS

gr.mirbase <- gr.mirna[gr.mirna$source == "miRBase_v22"]

# mirbase expressed primary mirnas
gr.mirbase.primir.expressed <- unique(gr.mirna[match(gr.mirbase.expressed$Derives_from, gr.mirna$ID)])


# gene annotations
rse.gene <- readRDS(rse.gene.rds)

# subset samples to those in vsd.mirna
rse.gene <- rse.gene[,colnames(vsd)]

# set seqnames for gene expression to UCSC
seqlevelsStyle(rse.gene) <- "UCSC"
genome(rse.gene) <- "hg38"
# remove non-standard chroms
seqlevels(rse.gene, pruning.mode = "coarse") <- CHROMS

# remove biotype = mirna
rse.gene <- rse.gene[rowRanges(rse.gene)$gene_biotype != "miRNA"]

# expressed or not
rowRanges(rse.gene)$expressed <- rowSums(assay(rse.gene) > 20) > 10

# # normalize expression (VST)
# dds.gene <- DESeqDataSet(rse.gene, design = ~1)
# vsd.gene <- vst(dds.gene)

# all gene annotations
gr.gene <- rowRanges(rse.gene)

gr.gene.expressed <- gr.gene[gr.gene$expressed]

# overlap

hits.primir.to.gene <- findOverlaps(gr.mirbase.primir.expressed, gr.gene.expressed)

# 525 expressed primary miRNAs
# 322 are overlapped with expressed genes (61.3%)




# mirQTLs ##############3
df.eqtl <- readRDS(eqtls.dataframe.rds)

df.eqtl.mirbase <- df.eqtl[df.eqtl$SOURCE == "miRBase_v22",]

emirs <- unique(df.eqtl.mirbase$NAME)

emirs.nums <- unique(sapply(strsplit(emirs, "-"), `[`, 3))


mirbase.nums <- unique(sapply(strsplit(gr.mirbase$Name, "-"), `[`, 3))

# Import Huan eQTLs ############

df.eqtl.huan <- read_xlsx(huan.eqtls.xlsx, skip = 1, na = "NA")

emirs.huan <- unique(df.eqtl.huan$hsa_miR_name)
emirs.huan <- emirs.huan[!is.na(emirs.huan)]


emirs.nums.huan <- unique(sapply(strsplit(emirs.huan, "-"), `[`, 3))

sum(emirs.huan %in% emirs)

sum(emirs %in% emirs.huan)


# overlap?

sum(emirs.nums %in% emirs.nums.huan)

sum()

# 109 miRBase emiRs
# 67 Huan emiRs

# 23 miRs overlap

not.emirs.nums <- mirbase.nums[!mirbase.nums %in% emirs.nums]



num.emirs <- length(emirs.nums)
num.emirs.notHuan <- sum(!emirs.nums %in% emirs.nums.huan)
num.emirs.Huan <- sum(emirs.nums %in% emirs.nums.huan)

num.not.emirs.Huan <- sum(not.emirs.nums %in% emirs.nums.huan)
num.not.emirs.notHuan <- sum(!not.emirs.nums %in% emirs.nums.huan)

ft <- fisher.test(matrix(c(num.emirs.Huan, num.emirs.notHuan,
                           num.not.emirs.Huan, num.not.emirs.notHuan),
                         byrow = TRUE,
                         nrow = 2))

ft

# Overlap with Gene annotations ################

# supportedFilters(EnsDb.Hsapiens.v86)
#
tx <- transcripts(EnsDb.Hsapiens.v86, filter = GeneIdFilter(gr.gene$gene_id[gr.gene$expressed]))



annotations <- build_annotations(genome = "hg38",
                                 annotations = c("hg38_genes_1to5kb", "hg38_genes_promoters", "hg38_genes_cds", "hg38_genes_5UTRs",
                                                 "hg38_genes_exons", "hg38_genes_introns", "hg38_genes_intronexonboundaries",
                                                 "hg38_genes_exonintronboundaries", "hg38_genes_3UTRs", "hg38_genes_intergenic"))

annotations$tx_id <- sapply(strsplit(annotations$tx_id, "\\."), `[`, 1)

# filter for expressed transcripts
annotations.expressed <- annotations[annotations$tx_id %in% names(tx)]

gr.annotated <- annotate_regions(regions = gr.mirbase.primir.expressed,
                                 annotations = annotations.expressed)

df.annotated <- data.frame(gr.annotated)

plot_annotation(gr.annotated)

ggsave("~/Desktop/annot_type.pdf")



