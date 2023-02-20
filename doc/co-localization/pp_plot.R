
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(mikelaffr)


# OUTPUT FILES ########################################################################################################
# output dir for pdf files
dir.pdf <- here("doc/co-localization/pdfs/")
dir.create(dir.pdf, showWarnings = FALSE, recursive = TRUE)

# INPUT FILES ##########################################################################################################
# Conditional eQTLs
conditional.eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# Summarized PRIMARY association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# GWAS Summary Stats: Educational Attain.
edu.sumstats.granges.rds <- here("data/gwas_datasets/educational_attainment/educational_attainment.hg38.GRanges.rds")

# Association results mRNA-eQTLs for chr14
mrna.results.chr14.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/raw/chr14_fetal_brain_mRNA-eQTL_raw.rds")

# Edu Att Index SNP, LD using eQTL population
eqtl.snps.rs10432209.ld <- here("results/coloc/plink_filesets/eqtl.snps.rs1043209.ld")

# GLOBALS ##############################################################################################################
# miR-4707-3p eQTL
chr <- "chr14"
eSNP <- "chr14:22953244:A:G"
emiR <- "hsa-mir-4707_hsa-miR-4707-3p"
esnp.bp <- 22953244

rsid.index.snp.gwas <- "rs1043209"

# Import Data ##################################################################################################

# Primary association data
df.eqtl.results <- as_tibble(readRDS(summarized.results.dataframe.rds))

# Filter for miR-4707-3p associations
df.eqtl.results %<>%
    filter(UniName == emiR)

max.bp <- max(df.eqtl.results$BP.hg38)
min.bp <- min(df.eqtl.results$BP.hg38)

# Educational Attain sum stats
gr.sumstats.edu <- readRDS(edu.sumstats.granges.rds)

range.overlap <- GRanges(seqnames = chr,
                         ranges = IRanges(start = min.bp, end = max.bp),
                         strand = "*")

gr.sumstats.edu <- subsetByOverlaps(gr.sumstats.edu, range.overlap)

df.sumstats.edu <- as_tibble(gr.sumstats.edu)

rm(gr.sumstats.edu, range.overlap)

df.eqtl.ld.to.index <- read_table(eqtl.snps.rs10432209.ld)
df.eqtl.ld.to.index %<>%
    select(SNP = SNP_B,
           R2_to_rs1043209_eqtl = R2)

# Combine Data ###############

df.eqtl.results %<>%
    select(SNP,
           CHR,
           BP = BP.hg38,
           REF,
           ALT,
           EFFECT.ALLELE,
           BETA,
           P)

df.sumstats.edu %<>%
    select(RSID,
           CHR = seqnames,
           BP = start,
           A1.GWAS = A1.effect,
           A2.GWAS = A2,
           EFFECT.ALLELE = A1.effect,
           BETA = Beta,
           P = Pval)

df.combined <- left_join(df.eqtl.results, df.sumstats.edu, by = c("CHR", "BP"), suffix = c(".eQTL", ".GWAS"))
df.combined %<>%
    filter(!is.na(RSID))

df.combined %<>%
    left_join(df.eqtl.ld.to.index, by = "SNP")


df.combined$color_to_rs1043209_eqtl <- "navy"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.8)] <- "red"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.6 & df.combined$R2_to_rs1043209_eqtl < 0.8)] <- "orange"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.4 & df.combined$R2_to_rs1043209_eqtl < 0.6)] <- "green"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.2 & df.combined$R2_to_rs1043209_eqtl < 0.4)] <- "lightblue"
df.combined$color_to_rs1043209_eqtl[which(df.combined$RSID == rsid.index.snp.gwas)] <- "purple"


# PP Plot ##############

df.combined %>%
    filter(BP > esnp.bp - 100000 & BP < esnp.bp + 100000) %>%
    mutate(color_to_rs1043209_eqtl = factor(color_to_rs1043209_eqtl, levels = c("purple", "red", "orange", "green", "lightblue", "navy"), ordered = TRUE)) %>%
    ggplot(aes(x = -log10(P.eQTL), y = -log10(P.GWAS))) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "navy"), color = "navy") +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "lightblue"), color = "lightblue") +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "green"), color = "green") +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "orange"), color = "orange") +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "red"), color = "red") +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "purple"), color = "purple") +
    theme(legend.position = "none") +
    plotTheme("figure")


df.combined %>%
    filter(BP > esnp.bp - 100000 & BP < esnp.bp + 100000) -> df.tmp

df.tmp %<>%
    mutate(l10p.eqtl = -log10(P.eQTL),
           l10p.gwas = -log10(P.GWAS))

p.cor <- cor.test(x = -log10(df.tmp$P.eQTL), y = -log10(df.tmp$P.GWAS), alternative = "two.sided")

pp.lm <- lm(l10p.gwas ~ l10p.eqtl, df.tmp)

summary(pp.lm)

int <- coef(pp.lm)[1]
slope <- coef(pp.lm)[2]
r2 <- summary(pp.lm)$r.squared

df.tmp %>%
    mutate(color_to_rs1043209_eqtl = factor(color_to_rs1043209_eqtl, levels = c("purple", "red", "orange", "green", "lightblue", "navy"), ordered = TRUE)) %>%
    ggplot(aes(x = -log10(P.eQTL), y = -log10(P.GWAS))) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "navy"), color = "navy", size = 1) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "lightblue"), color = "lightblue", size = 1) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "green"), color = "green", size = 1) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "orange"), color = "orange", size = 1) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "red"), color = "red", size = 1) +
    geom_point(data = function(x) subset(x, color_to_rs1043209_eqtl == "purple"), color = "purple", size = 1) +
    geom_point(data = function(x) subset(x, RSID == "rs4981455"), color = "purple", size = 1) +
    theme(legend.position = "none") +
    geom_text(x = 1, y = 12, label = paste0("pearson corr: ", format(p.cor$estimate, digits = 3), "\np-value: ", format(p.cor$p.value, digits = 2))) +
    theme(axis.title = ggplot2::element_text(size = 6,
                                             face = "bold"),
          axis.text = ggplot2::element_text(size = 5),
          title = ggplot2::element_text(size = 6,
                                        face = "bold"),
          legend.title = ggplot2::element_text(size = 6),
          legend.text = ggplot2::element_text(size = 6),
          panel.background = ggplot2::element_rect(fill = "white",
                                                   linetype = "solid",
                                                   color = "black",
                                                   linewidth = 1),
          panel.grid = ggplot2::element_blank(),) +
    scale_x_continuous(expand = expansion(.01), limits = c(0,15)) +
    scale_y_continuous(expand = expansion(.01), limits = c(0,15))

ggsave(paste0(dir.pdf, "pp-plot2.pdf"), height = 3, width = 3)





