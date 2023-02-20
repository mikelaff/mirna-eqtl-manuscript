
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(mikelaffr)
library(here)

library(rtracklayer)
library(liftOver)

# OUTPUT ##############################################################################################################
dir.pdfs <- here("doc/imputation/pdfs/")

# INPUT ###############################################################################################################
# topmed imputation (subsampled to 1% of imputed snps on each chrom)
df.info.topmed.tsv <- here("results/genotypes/imputation_topmed_freeze5/results/df_info_subsample_topmed.tsv.gz")
df.info.typed.topmed.tsv <- here("results/genotypes/imputation_topmed_freeze5/results/df_info_typed_topmed.tsv.gz")
df.info.low.topmed.tsv <- here("results/genotypes/imputation_topmed_freeze5/results/df_info_low_topmed.tsv.gz")

# 1kg imputation (subsampled to 10% of imputed snps on each chrom)
df.info.1kg.tsv <- here("results/genotypes/imputation_1000G_phase3v5/results/df_info_subsample_1kg.tsv.gz")
df.info.typed.1kg.tsv <- here("results/genotypes/imputation_1000G_phase3v5/results/df_info_typed_1kg.tsv.gz")
df.info.low.1kg.tsv <- here("results/genotypes/imputation_1000G_phase3v5/results/df_info_low_1kg.tsv.gz")

# GLOBALS #############################################################################################################


# Load Data ###########################################################################################################
# load subsampled imputated snp data
df.info.1kg <- read_tsv(df.info.1kg.tsv)
df.info.topmed <- read_tsv(df.info.topmed.tsv)

# load all typed imputed snp data
df.info.typed.1kg <- read_tsv(df.info.typed.1kg.tsv)
df.info.typed.topmed <- read_tsv(df.info.typed.topmed.tsv)

# load low MAF imputed snp data
df.info.low.1kg <- read_tsv(df.info.low.1kg.tsv)
df.info.low.topmed <- read_tsv(df.info.low.topmed.tsv)

# label sets
df.info.1kg$set <- "1kg"
df.info.typed.1kg$set <- "1kg"
df.info.low.1kg$set <- "1kg"

df.info.topmed$set <- "topmed"
df.info.typed.topmed$set <- "topmed"
df.info.low.topmed$set <- "topmed"

# combine tables
df.info <- bind_rows(df.info.1kg, df.info.topmed)
df.info.typed <- bind_rows(df.info.typed.1kg, df.info.typed.topmed)
df.info.low <- bind_rows(df.info.low.1kg, df.info.low.topmed)
rm(df.info.1kg, df.info.topmed, df.info.typed.1kg, df.info.typed.topmed, df.info.low.1kg, df.info.low.topmed)

# Plot ################################################################################################################

# imputed SNPs

df.info %>%
    ggplot(mapping = aes(x = MAF, y = Rsq, color = set)) +
    #geom_point(data = sample_frac(df.info, 0.01), alpha = 1/5) +
    stat_summary_bin(data = sample_frac(df.info, 1), fun.y = "mean", geom = "line", size = 1.5, bins = 50) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2)) +
    scale_color_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    labs(caption = "50bins",
         title = "Imputed SNPs (binned)",
         color = "Reference Panel",
         x = "Minor Allele Frequency",
         y = "Mean Rsq")

ggsave(paste0(dir.pdfs, "pres_imputed_snps_rsq_by_maf.pdf"), height = 5.5, width = 5.5)

df.info %>%
    ggplot(mapping = aes(x = MAF, y = AvgCall, color = set)) +
    #geom_point(data = sample_frac(df.info, 0.01), alpha = 1/5) +
    stat_summary_bin(data = sample_frac(df.info, 1), fun.y = "mean", geom = "line", size = 1.5, bins = 50) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.05, .05),
          legend.justification = c("left", "bottom"),
          legend.box.just = "left",
          legend.margin = margin(2, 2, 2, 2)) +
    scale_color_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    labs(caption = "50bins",
         title = "Imputed SNPs (binned)",
         color = "Reference Panel",
         x = "Minor Allele Frequency",
         y = "Mean AvgCall")

ggsave(paste0(dir.pdfs, "pres_imputed_snps_avgcall_by_maf.pdf"), height = 5.5, width = 5.5)

df.info %>%
    filter(Rsq > 0.3) %>%
    ggplot(mapping = aes(x = MAF, y = Rsq, color = set)) +
    stat_summary_bin(fun.y = "mean", geom = "line", size = 1.5, bins = 50) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2)) +
    scale_color_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    labs(caption = "50bins",
         title = "Imputed SNPs: Rsq > 0.3",
         color = "Reference Panel",
         x = "Minor Allele Frequency",
         y = "Mean Rsq")

ggsave(paste0(dir.pdfs, "pres_imputed_snps_filtered_rsq_by_maf.pdf"), height = 5.5, width = 5.5)

# typed SNPs

df.info.typed %>%
    ggplot(mapping = aes(x = MAF, y = EmpRsq, color = set)) +
    #geom_point(data = sample_frac(df.info.typed, 0.01), alpha = 1/5) +
    stat_summary_bin(data = sample_frac(df.info.typed, 1), fun.y = "mean", geom = "line", size = 1.5, bins = 50) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2)) +
    scale_color_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    labs(caption = "50bins",
         title = "Typed SNPs (binned)",
         color = "Reference Panel",
         x = "Minor Allele Frequency",
         y = "Mean EmpRsq")

ggsave(paste0(dir.pdfs, "pres_typed_snps_rsq_by_maf.pdf"), height = 5.5, width = 5.5)


# imputed snps quantity

df.info %>%
    mutate(rsq.bin = cut(Rsq, breaks = seq(0,1,.1), include.lowest = TRUE)) %>%
    group_by(rsq.bin, set) %>%
    summarise(n = n()) %>%
    ggplot(aes(x = rsq.bin, y = n, fill = set)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5"))

#
#     ggplot(aes(Rsq, fill = set)) +
#     #geom_histogram(bins = 50, position = "dodge") +
#     geom_density() +
#     plotTheme() +
#     scale_fill_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5"))

# zoom in

df.info.low %>%
    filter(MAF >= 0.001, MAF <= 0.05) %>%
    ggplot(mapping = aes(x = MAF, y = Rsq, color = set)) +
    stat_summary_bin(data = sample_frac(df.info, 1), fun.y = "mean", geom = "line", size = 1.5, bins = 10) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .05),
          legend.justification = c("right", "bottom"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2)) +
    scale_color_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    labs(caption = "10 bins",
         title = "Imputed SNPs: low MAF",
         color = "Reference Panel",
         x = "Minor Allele Frequency",
         y = "Mean Rsq") +
    scale_x_continuous(limits = c(0.001, 0.05))

ggsave(paste0(dir.pdfs, "pres_imputed_snps_rsq_by_maf_lowMAF.pdf"), height = 5.5, width = 5.5)




# Scratch #############################################################################################################
stop()

p.hat <- seq(from=0, to=1, by=.01)
one.minus.p.hat <- 1-p.hat
denom <- p.hat * one.minus.p.hat

plot(p.hat, denom)

n.haplotypes <- 328 * 2
d.prob <- runif(n = n.haplotypes, min = 0.8, max = 0.9)

r.hat.squared <- sapply(p.hat, function(x) ( (1 / n.haplotypes) * sum( (d.prob - x)^2 ) ) / (x * (1 - x)))

plot(p.hat, r.hat.squared)

# Number of Variants ###################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")


df.1kg.hard <- read_table(here("results/genotypes/imputation_1000G_phase3v5/results/1kg.hardcall.num.vars"), col_names = c("num.vars", "file"), n_max = 23)
df.1kg.hard$chr <- sapply(strsplit(df.1kg.hard$file, "\\."), `[`, 1)
df.1kg.hard$chr <- factor(df.1kg.hard$chr, levels = CHROMS, ordered = TRUE)
df.1kg.hard$source <- "1kg"
df.1kg.hard$type <- "hardcall"

df.1kg.dose <- read_table(here("results/genotypes/imputation_1000G_phase3v5/results/1kg.dosage.num.vars"), col_names = c("num.vars", "file"), n_max = 23)
df.1kg.dose$chr <- sapply(strsplit(df.1kg.dose$file, "\\."), `[`, 1)
df.1kg.dose$chr <- factor(df.1kg.dose$chr, levels = CHROMS, ordered = TRUE)
df.1kg.dose$source <- "1kg"
df.1kg.dose$type <- "dosage"

df.topmed.hard <- read_table(here("results/genotypes/imputation_topmed_freeze5/results/topmed.hardcall.num.vars"), col_names = c("num.vars", "file"), n_max = 23)
df.topmed.hard$chr <- sapply(strsplit(df.topmed.hard$file, "\\."), `[`, 1)
df.topmed.hard$chr <- factor(df.topmed.hard$chr, levels = CHROMS, ordered = TRUE)
df.topmed.hard$source <- "topmed"
df.topmed.hard$type <- "hardcall"

df.topmed.dose <- read_table(here("results/genotypes/imputation_topmed_freeze5/results/topmed.dosage.num.vars"), col_names = c("num.vars", "file"), n_max = 23)
df.topmed.dose$chr <- sapply(strsplit(df.topmed.dose$file, "\\."), `[`, 1)
df.topmed.dose$chr <- factor(df.topmed.dose$chr, levels = CHROMS, ordered = TRUE)
df.topmed.dose$source <- "topmed"
df.topmed.dose$type <- "dosage"

df <- bind_rows(df.1kg.hard, df.1kg.dose, df.topmed.hard, df.topmed.dose)

df %>%
    filter(type == "hardcall") %>%
    ggplot(aes(x = chr, y = num.vars/1e6, fill = source)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    plotTheme(theme = "presentation",
              legend.position = "bottom",
              axis.text.x = element_text(size = 14, angle = 45, vjust = 0.6)) +
    scale_y_continuous(expand = expand_scale(c(0,.05))) +
    labs(y = "Number of Variants (million)",
         x = NULL,
         fill = "Reference Panel",
         title = "Imputed Variants after QC Filtering",
         caption = "QC Filters: R2 > 0.3, MAF > 0.01, HWE > 1e-6")

ggsave(paste0(dir.pdfs, "imputed_vars_post_filt.pdf"), width = 11, height = 6)

df %>%
    filter(type == "hardcall") %>%
    group_by(source) %>%
    summarise(total = sum(num.vars))


df.num <- read_csv(here("results/genotypes/imputation_num_vars.csv"))

df.num %>%
    ggplot(aes(x=factor(filter, labels = c("Raw", "R2 > 0.3", "R2 > 0.3\n+\nMAF > 0.01", "R2 > 0.3\n+\nMAF > 0.01\n+\nHWE > 1e-6")), y=num.vars/1e6, fill = source)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2),
          axis.text.x = element_text(size = 14)) +
    scale_y_continuous(expand = expand_scale(c(0,.07))) +
    labs(y = "Number of Variants (million)",
         x = NULL,
         fill = "Reference Panel",
         title = "Imputed Variants") +
    geom_text(aes(label=round(num.vars/1e6, digits = 1)), vjust=-0.5, size=5, fontface=2, position = position_dodge(width = 1))

ggsave(paste0(dir.pdfs, "total_imputed_vars_by_filt.pdf"), width = 6, height = 7)

# Compare MAF 1KG to TOPMed #################

# load subsampled imputated snp data
df.info.1kg <- read_tsv(df.info.1kg.tsv)
df.info.topmed <- read_tsv(df.info.topmed.tsv)

# label sets
df.info.1kg$set <- "1kg"
df.info.topmed$set <- "topmed"

# get base pair
df.info.1kg$BP <- as.numeric(sapply(strsplit(df.info.1kg$SNP, ":"), `[`, 2))
df.info.topmed$BP <- as.numeric(sapply(strsplit(df.info.topmed$SNP, ":"), `[`, 2))

# converto to granges
gr.1kg.hg19 <- GRanges(seqnames = df.info.1kg$CHR,
                       ranges = IRanges(start = df.info.1kg$BP,
                                        end = df.info.1kg$BP),
                       seqinfo = Seqinfo(genome = "hg19"),
                       mcols = df.info.1kg)

gr.topmed.hg38 <- GRanges(seqnames = df.info.topmed$CHR,
                          ranges = IRanges(start = df.info.topmed$BP,
                                           end = df.info.topmed$BP),
                          seqinfo = Seqinfo(genome = "hg38"),
                          mcols = df.info.topmed)

# chain for hg19 to hg38 conversion
path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch <- import.chain(path)
# liftOver hg19 to hg38
grl.1kg.hg38 <- liftOver(gr.1kg.hg19, ch)
gr.1kg.hg38 <- unlist(grl.1kg.hg38)


# find overlaps
overlaps <- findOverlaps(gr.1kg.hg38, gr.topmed.hg38)

gr.1kg.overlaps <- gr.1kg.hg38[queryHits(overlaps)]
gr.topmed.overlaps <- gr.topmed.hg38[subjectHits(overlaps)]

df.1kg.overlaps <- as.data.frame(gr.1kg.overlaps)
df.topmed.overlaps <- as.data.frame(gr.topmed.overlaps)

df.1kg.overlaps$pos <- paste(df.1kg.overlaps$seqnames, df.1kg.overlaps$start, sep = ":")
df.topmed.overlaps$pos <- paste(df.topmed.overlaps$seqnames, df.topmed.overlaps$start, sep = ":")

df.combo <- full_join(df.1kg.overlaps, df.topmed.overlaps, by = "pos", suffix = c(".1kg", ".topmed"))
df.combo$allele.match <- df.combo$mcols.REF.0..1kg == df.combo$mcols.REF.0..topmed & df.combo$mcols.ALT.1..1kg == df.combo$mcols.ALT.1..topmed

sum(df.combo$allele.match)

pdf(paste0(dir.pdfs, "overlaping_variants.pdf"), width = 7, height = 7)

df.combo %>%
    ggplot(aes(x = mcols.MAF.1kg, y = mcols.MAF.topmed)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "MAF : 1KG",
         y = "MAF : TOPMed",
         title = "Variants in Both Datasets: Unfiltered")

df.combo %>%
    filter(mcols.Rsq.1kg > 0.3 & mcols.Rsq.topmed > 0.3) %>%
    ggplot(aes(x = mcols.MAF.1kg, y = mcols.MAF.topmed)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "MAF : 1KG",
         y = "MAF : TOPMed",
         title = "Variants in Both Datasets: R2 > 0.3")

df.combo %>%
    filter(mcols.Rsq.1kg > 0.3 & mcols.Rsq.topmed > 0.3) %>%
    filter(mcols.MAF.1kg > 0.01 & mcols.MAF.topmed > 0.01) %>%
    ggplot(aes(x = mcols.MAF.1kg, y = mcols.MAF.topmed)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "MAF : 1KG",
         y = "MAF : TOPMed",
         title = "Variants in Both Datasets: R2 > 0.3 + MAF > 0.01")

df.combo %>%
    filter(mcols.Rsq.1kg > 0.3 & mcols.Rsq.topmed > 0.3) %>%
    filter(mcols.MAF.1kg < 0.05 & mcols.MAF.topmed < 0.05) %>%
    ggplot(aes(x = mcols.MAF.1kg, y = mcols.MAF.topmed)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "MAF : 1KG",
         y = "MAF : TOPMed",
         title = "Variants in Both Datasets: R2 > 0.3, Low Frequency")

df.combo %>%
    filter(mcols.Rsq.1kg > 0.3 & mcols.Rsq.topmed > 0.3) %>%
    filter(mcols.MAF.1kg < 0.02 & mcols.MAF.topmed < 0.02) %>%
    ggplot(aes(x = mcols.MAF.1kg, y = mcols.MAF.topmed)) +
    geom_point(size = 0.5) +
    geom_smooth(method = "lm") +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "MAF : 1KG",
         y = "MAF : TOPMed",
         title = "Variants in Both Datasets: R2 > 0.3, Low Frequency")

df.combo %>%
    ggplot(aes(x = mcols.Rsq.1kg, y = mcols.Rsq.topmed)) +
    geom_point(size = 0.5) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    plotTheme() +
    scale_y_continuous(expand = expand_scale(c(0,0))) +
    scale_x_continuous(expand = expand_scale(c(0,0))) +
    labs(x = "Rsq : 1KG",
         y = "Rsq : TOPMed",
         title = "Variants in Both Datasets: Unfiltered")

dev.off()


# Lego Plot ############
library(brickr)

df.num %>%
    ggplot(aes(x=factor(filter, labels = c("Raw", "R2 > 0.3", "R2 > 0.3\n+\nMAF > 0.01", "R2 > 0.3\n+\nMAF > 0.01\n+\nHWE > 1e-6")), y=num.vars/1e6, fill = source)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = cbPalette[c(1,2)], labels = c("1000G Phase3v5", "TOPMed Freeze5")) +
    plotTheme(theme = "presentation") +
    theme(legend.position = c(.95, .95),
          legend.justification = c("right", "top"),
          legend.box.just = "right",
          legend.margin = margin(2, 2, 2, 2),
          axis.text.x = element_text(size = 14)) +
    scale_y_continuous(expand = expand_scale(c(0,.07))) +
    labs(y = "Number of Variants (million)",
         x = NULL,
         fill = "Reference Panel",
         title = "Imputed Variants") +
    geom_text(aes(label=round(num.vars/1e6, digits = 1)), vjust=-0.5, size=5, fontface=2, position = position_dodge(width = 1))

df.num %>%
    ggplot(aes(x = filter, y = num.vars/1e6)) +
    geom_brick_col(aes(fill = source)) +
    scale_fill_brick() +
    coord_brick() +
    theme_brick()

ggsave("~/Downloads/brick_plot.png", height = 6, width = 6)




