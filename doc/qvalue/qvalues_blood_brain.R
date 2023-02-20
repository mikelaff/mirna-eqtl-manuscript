# qvalues or pi one values?

library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(qvalue)
library(boot)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/qvalue/pdfs/")
dir.create(dir.pdfs, showWarnings = FALSE, recursive = TRUE)


# INPUT ################################################################################################################
# brain miRNA-eQTL primary results
rank1.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# blood miRNA-eQTLs
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# brain mRNA-eQTL raw data dir
mrna.eqtl.raw.dir <- here("results/external_data/fetal_brain_mRNA-eQTL/raw/")

# GTEx whole blood sig variants
gtex.sig.variants.txt <- here("data/GTEx/Whole_Blood.v8.signif_variant_gene_pairs.txt.gz")
# GTEx whole blood egenes
gtex.egenes.txt <- here("data/GTEx/Whole_Blood.v8.egenes.txt.gz")

# hg19 to hg38 chain file
chain.file <- here("data/hg19ToHg38.over.chain")


# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import and Lift Blood eQTLs ##########################################################################################

df.blood <- read_xlsx(huan.eqtls.xlsx, skip = 1)

# chain for hg19 to hg38 conversion
ch <- import.chain(chain.file)

# GRanges of huan eqtls
gr <- makeGRangesFromDataFrame(df.blood,
                               seqnames.field = "chr.SNP",
                               start.field = "SNP.pos",
                               end.field = "SNP.pos",
                               strand.field = "SNP.strand",
                               ignore.strand = TRUE,
                               keep.extra.columns = TRUE)
# GRangesList of GRanges conversion
lo <- liftOver(gr, ch)

gr.huan.hg38 <- unlist(lo)

rm(lo, gr, ch)

df.blood <- as_tibble(as.data.frame(gr.huan.hg38))

df.blood %<>%
    mutate(variant.pos.name = paste(seqnames, start, hsa_miR_name, sep = ":"),
           SNP.position = paste(seqnames, start, sep = ":"))

df.blood %>%
    group_by(hsa_miR_name) %>%
    filter(Pval == min(Pval)) -> df.blood.minP


# Import Brain eQTLs ###################################################################################################

df.results <- as_tibble(read_rds(rank1.results.dataframe.rds))

df.results %<>%
    mutate(Name = sapply(strsplit(UniName, "_"), `[`, 2))

df.results %<>%
    mutate(variant.pos.name = paste(CHR, BP.hg38, Name, sep = ":"),
           SNP.position = paste(CHR, BP.hg38, sep = ":"))

# Shared eQTLs #########################################################################################################
df.results %>%
    filter(variant.pos.name %in% df.blood.minP$variant.pos.name) -> df.shared

# mirs that didnt have a corrisponding variant
unique(df.results$Name[df.results$Name %in% df.blood.minP$hsa_miR_name])[!unique(df.results$Name[df.results$Name %in% df.blood.minP$hsa_miR_name]) %in% df.shared$Name]

# df.results %>%
#     filter(Name %in% df.blood.minP$hsa_miR_name) -> df.shared
#
# df.shared %<>%
#     filter(SNP.position %in% df.blood.minP$SNP.position & Name %in% df.blood.minP$hsa_miR_name)
#
# df.shared %<>%
#     group_by(Name) %>%
#     filter(P == min(P))

# df.blood %>%
#     filter(SNP.position %in% df.results$SNP.position) -> df.shared.blood


# df.shared %>%
#     group_by(SNP.position) %>%
#     filter(P == min(P)) -> df.shared.minP

# qvalue ###############################################################################################################

#q.results <- qvalue(p = df.shared$P, lambda = seq(0.05, 0.95, 0.05))
q.results <- qvalue(p = df.shared$P, pi0.method = "bootstrap", lambda = seq(0.1, 0.7, 0.05))

pi1 <- 1 - q.results$pi0

# q.results.minP <- qvalue(p = df.shared.minP$P)
#
# pi1.minP <- 1 - q.results.minP$pi0

# q.results.blood <- qvalue(p = df.shared.blood$Pval)

# bootstrap ############################################################################################################

set.seed(34)

bstrap <- boot(data = df.shared$P, function(x, i) 1-qvalue(x[i], pi0.method = "bootstrap", lambda = seq(0.1, 0.7, 0.05))$pi0, R = 100)

# std. error
std.error <- sd(bstrap$t)

# 95% conf interval
pi1.low <- pi1 - 1.96 * std.error

pi1.high <- pi1 + 1.96 * std.error

c(pi1.low, pi1, pi1.high)

print(paste("Pi_1 value in miRNA-eQTL analysis (+/- 95% conf):", signif(pi1, digits = 3), "+/-", signif(1.96 * std.error, digits = 3)))
print(paste("Number of p-values tested in mRNA-eQTL analysis:", length(df.shared$P)))

df.plot <- tibble(dataset = "miRNA", pi_1 = pi1, low = pi1.low, high = pi1.high)

# mRNA Data ############################################################################################################

# # GTEx sig variants
# df.gtex.variants <- read_tsv(gtex.sig.variants.txt)
#
# df.gtex.variants %>%
#     group_by(gene_id)

# GTEx egenes
df.gtex.egenes <- read_tsv(gtex.egenes.txt)
# significant genes
df.gtex.egenes %<>%
    filter(qval <= 0.05)
# remove ENSG version number
df.gtex.egenes %<>%
    mutate(gene_id = sapply(strsplit(gene_id, "\\."), `[`, 1))
# unique variant/gene identifier (eSNP-eGene or eQTL)
df.gtex.egenes %<>%
    mutate(variant.pos.ensg = paste(chr, variant_pos, gene_id, sep = ":"))

df.mrna.results <- tibble()

# loop over chroms and load p-vals for mRNA data at the gtex index positions
for (chr in CHROMS) {
    printMessage(chr)

    df.mrna.tmp <- NULL

    # brain mrna data
    df.mrna.tmp <- read_rds(paste0(mrna.eqtl.raw.dir, chr, "_fetal_brain_mRNA-eQTL_raw.rds"))

    # unique variant/gene identifier
    df.mrna.tmp %<>%
        mutate(variant.pos.ensg = paste(CHR, BP, ENSG, sep = ":"))

    # filter for only gtex eqtls
    df.mrna.tmp %<>%
        filter(variant.pos.ensg %in% df.gtex.egenes$variant.pos.ensg)

    df.mrna.results %<>%
        bind_rows(df.mrna.tmp)

}

rm(df.mrna.tmp)

q.results.mrna <- qvalue(p = df.mrna.results$P, pi0.method = "bootstrap")

pi1.mrna <- 1 - q.results.mrna$pi0

# q.results.minP <- qvalue(p = df.shared.minP$P)
#
# pi1.minP <- 1 - q.results.minP$pi0

# q.results.blood <- qvalue(p = df.shared.blood$Pval)

set.seed(34)

bstrap.mrna <- boot(data = df.mrna.results$P, function(x, i) 1-qvalue(x[i])$pi0, R = 100)

# std. error
std.error.mrna <- sd(bstrap.mrna$t)

# 95% conf interval
pi1.low.mrna <- pi1.mrna - 1.96 * std.error.mrna

pi1.high.mrna <- pi1.mrna + 1.96 * std.error.mrna

c(pi1.low.mrna, pi1.mrna, pi1.high.mrna)

print(paste("Pi_1 value in mRNA-eQTL analysis (+/- 95% conf):", signif(pi1.mrna, digits = 3), "+/-", signif(1.96 * std.error.mrna, digits = 3)))
print(paste("Number of p-values tested in mRNA-eQTL analysis:", length(df.mrna.results$P)))



df.plot %<>%
    bind_rows(tibble(dataset = "mRNA", pi_1 = pi1.mrna, low = pi1.low.mrna, high = pi1.high.mrna))

df.plot %>%
    ggplot(aes(x = dataset, y = pi_1)) +
    geom_point() +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.1) +
    plotTheme("figure")

ggsave("~/Desktop/pi_1_plot.pdf", width = 2, height = 3)





