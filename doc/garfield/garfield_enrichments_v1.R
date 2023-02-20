# look at garfield results

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
dir.pdf <- here("doc/garfield/pdfs/")

# INPUT FILES ##########################################################################################################
# mirQTL enrichment on chromHMM using 1000G-EUR maf/ld
eur.results.txt <- here("results/garfield/hg38-1000G-EUR_mirQTL_chromHMM/archive/garfield.test.hg38-1000G-EUR_mirQTL_chromHMM.out")
# p-value file
eur.pvalue.txt <- here("results/garfield/hg38-1000G-EUR_mirQTL_chromHMM/archive/garfield.Meff.hg38-1000G-EUR_mirQTL_chromHMM.out")

# mirQTL enrichment on chromHMM using mirQTL-MIXED maf/ld
mix.results.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL_chromHMM/archive/garfield.test.hg38-mirQTL-MIXED_mirQTL_chromHMM.out")
# p-value file
mix.pvalue.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL_chromHMM/archive/garfield.Meff.hg38-mirQTL-MIXED_mirQTL_chromHMM.out")

# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################
df.eur.results <- read_table2(eur.results.txt)
eur.pvalue.thresh <- as.numeric(strsplit(read_lines(eur.pvalue.txt), "\t")[[2]][2])

df.mix.results <- read_table2(mix.results.txt)
mix.pvalue.thresh <- as.numeric(strsplit(read_lines(mix.pvalue.txt), "\t")[[2]][2])

# Plot #################################################################################################################

df.eur.results %<>%
    mutate(Annotation_Name = sapply(strsplit(Annotation, "_"), `[`, 2),
           Annotation_Sex = sapply(strsplit(Annotation, "_"), `[`, 3),
           Annotation_Number = as.integer(sapply(strsplit(Annotation, "_"), `[`, 1)),
           SIG = Pvalue <= eur.pvalue.thresh)

df.mix.results %<>%
    mutate(Annotation_Name = sapply(strsplit(Annotation, "_"), `[`, 2),
           Annotation_Sex = sapply(strsplit(Annotation, "_"), `[`, 3),
           Annotation_Number = as.integer(sapply(strsplit(Annotation, "_"), `[`, 1)),
           SIG = Pvalue <= mix.pvalue.thresh)

pdf(file = paste0(dir.pdf, "mirQTL_pval_chromHMM_annotation_mirQTL-MAFLD_and_1000G-MAFLD.pdf"), width = 10, height = 4)

# mirQTL MIXED
df.mix.results %>%
    filter(Annotation_Name != "TxFlnk") %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(0.7), shape = 15, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = cbPalette) +
    geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
    geom_hline(yintercept = 0) +
    labs(y = "Log Odds Ratio (95% Conf. Int.)",
         x = "ChromHMM 15-state Mnemonic",
         color = "Sex",
         alpha = "Significant Enrichment P-value",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         caption = "TxFlnk annotation removed (high error). Using mirQTL (MIXED ancestry) MAF and LD. P-value threshold for mirQTL variants: 2.804e-5")

df.mix.results %>%
    filter(Annotation_Name != "TxFlnk") %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = cbPalette) +
    geom_hline(aes(yintercept = -log10(mix.pvalue.thresh)), linetype = "dashed") +
    labs(y = "-Log10(P-value)",
         x = "ChromHMM 15-state Mnemonic",
         fill = "Sex",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         caption = "TxFlnk annotation removed (high error). Using mirQTL (MIXED ancestry) MAF and LD. P-value threshold for mirQTL variants: 2.804e-5")

# 1000Genomes EUR
df.eur.results %>%
    filter(Annotation_Name != "TssBiv") %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(0.7), shape = 15, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = cbPalette) +
    geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
    geom_hline(yintercept = 0) +
    labs(y = "Log Odds Ratio (95% Conf. Int.)",
         x = "ChromHMM 15-state Mnemonic",
         color = "Sex",
         alpha = "Significant Enrichment P-value",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         caption = "TssBiv annotation removed (high error). Using 1KGenomes (EUR ancestry) MAF and LD. P-value threshold for mirQTL variants: 2.804e-5")

df.eur.results %>%
    filter(Annotation_Name != "TssBiv") %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = cbPalette) +
    geom_hline(aes(yintercept = -log10(eur.pvalue.thresh)), linetype = "dashed") +
    labs(y = "-Log10(P-value)",
         x = "ChromHMM 15-state Mnemonic",
         fill = "Sex",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         caption = "TssBiv annotation removed (high error). Using 1KGenomes (EUR ancestry) MAF and LD. P-value threshold for mirQTL variants: 2.804e-5")

dev.off()




#ggsave(paste0(dir.pdf, "chrom_states_or.pdf"), width = 10, height = 6)

