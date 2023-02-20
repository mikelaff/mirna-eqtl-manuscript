# sex calls by genotype and xist expression

library(dplyr)
library(magrittr)
library(readr)
library(here)
library(ggplot2)
library(DESeq2)
library(mikelaffr)

# OUTPUT FILES ########################################################################################################
dir.pdf <- here("doc/eda/pdfs")

# INPUT FILES #########################################################################################################
# this should be a database query, however, here is the table of all small rna-seq samples
samples.small.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")
samples.total.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_totalRNAseq_metadata.tsv")

# plink .fam file to get sexes by genotype
samples.fam <- here("results/genotypes/AllSamplesQC.fam")

# plink sexcheck file
sexcheck.file <- here("results/genotypes/AllSamplesQC.sexcheck")

# total-rna-seq expression data
total.rna.exprs.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# Import Data ###################

rse <- readRDS(total.rna.exprs.rds)
dds <- DESeqDataSet(rse, design = ~1)

# get normalized XIST expression across all samples
# XIST = ENSG00000229807
df.xist <- data.frame(xist = assay(normTransform(dds))["ENSG00000229807",],
                      rnaid = names(assay(normTransform(dds))["ENSG00000229807",]),
                      donorid = paste0("D",dds$donor_id),
                      dnaid = dds$dnaid,
                      stringsAsFactors = FALSE)

# call sex
df.xist$sex_by_xist <- ifelse(df.xist$xist > 13, "Female", "Male")

# plot
df.xist %>%
    ggplot(aes(x=reorder(rnaid, xist), y=xist, color=sex_by_xist)) +
    geom_point() +
    geom_hline(yintercept = 13) +
    labs(x = "Rank Expression, Fetal Tissue Samples",
         y = "log2(Scaled XIST Expression)",
         #title = "Fetal Tissue Total RNA-Seq",
         color = "Sex By XIST") +
    plotTheme(theme = "presentation",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
    scale_color_manual(values = c("darkorange", "navy"))

ggsave(paste(dir.pdf, "presentation_sex_by_xist.pdf", sep="/"), height = 5, width = 8)

df.xist %>%
    ggplot(aes(x=reorder(rnaid, xist), y=xist, color=sex_by_xist)) +
    geom_point() +
    geom_hline(yintercept = 13) +
    labs(x = "Rank Expression, Fetal Tissue Samples",
         y = "log2(Scaled XIST Expression)",
         #title = "Fetal Tissue Total RNA-Seq",
         color = "Sex By XIST") +
    plotTheme(theme = "presentation",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
    scale_color_manual(values = c("darkorange", "navy")) +
    geom_point(aes(x=reorder(rnaid, xist), y=xist), data=subset(df.xist, rnaid == "RNAID1561"), color="black", shape=21, size=4, stroke=2) +
    geom_point(aes(x=reorder(rnaid, xist), y=xist), data=subset(df.xist, rnaid == "RNAID1560"), color="black", shape=21, size=4, stroke=2)

ggsave(paste(dir.pdf, "presentation_sex_by_xist_labs.pdf", sep="/"), height = 6, width = 6)




sex.check <- read_table(here("results/genotypes/AllSamplesQC.sexcheck"))

samps <- read_tsv(here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv"))

tmp <- left_join(samps, sex.check, by = c("VerifiedDNAID" = "IID"))

tmp <- left_join(tmp, df.xist, by = c("RNAID" = "rnaid"))

tmp %>%
    ggplot(aes(x=reorder(RNAID, `F`), y=`F`)) +
    geom_point() +
    plotTheme(theme = "presentation",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank()) +
    labs(x = "240 Small RNA-Seq Samples")

ggsave(paste(dir.pdf, "f_by_samples.pdf", sep="/"), height=6, width = 6)

tmp %>%
    ggplot(aes(y=`F`, x=xist, color=Sex.by.XIST)) +
    geom_point() +
    plotTheme(theme = "presentation",
              axis.text.x = element_blank(),
              axis.ticks.x = element_blank(),
              legend.position = "bottom") +
    scale_color_manual(values = c("darkorange", "navy")) +
    labs(x = "log2(Scaled XIST Expression)",
         color = "Sex By XIST")

ggsave(paste(dir.pdf, "f_by_xist.pdf", sep="/"), height=6, width = 6)
