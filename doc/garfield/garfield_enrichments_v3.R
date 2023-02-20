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
eur.results.txt <- here("results/garfield/hg38-1000G-EUR_mirQTL_chromHMM/garfield.test.hg38-1000G-EUR_mirQTL_chromHMM.out")
# p-value file
eur.pvalue.txt <- here("results/garfield/hg38-1000G-EUR_mirQTL_chromHMM/garfield.Meff.hg38-1000G-EUR_mirQTL_chromHMM.out")

# mirQTL enrichment on chromHMM using mirQTL-MIXED maf/ld
mix.results.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL_chromHMM/garfield.test.hg38-mirQTL-MIXED_mirQTL_chromHMM.out")
# p-value file
mix.pvalue.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL_chromHMM/garfield.Meff.hg38-mirQTL-MIXED_mirQTL_chromHMM.out")

# mirQTL enrichment on chromHMM using mirQTL-MIXED maf/ld and minimum P-value
minP.results.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL-minP_chromHMM/garfield.test.hg38-mirQTL-MIXED_mirQTL-minP_chromHMM.out")
# p-value file
minP.pvalue.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL-minP_chromHMM/garfield.Meff.hg38-mirQTL-MIXED_mirQTL-minP_chromHMM.out")

# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################
df.eur.results <- read_table2(eur.results.txt)
eur.pvalue.thresh <- as.numeric(strsplit(read_lines(eur.pvalue.txt), "\t")[[2]][2])

df.mix.results <- read_table2(mix.results.txt)
mix.pvalue.thresh <- as.numeric(strsplit(read_lines(mix.pvalue.txt), "\t")[[2]][2])

df.minP.results <- read_table2(minP.results.txt)
minP.pvalue.thresh <- as.numeric(strsplit(read_lines(minP.pvalue.txt), "\t")[[2]][2])

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

df.minP.results %<>%
    mutate(Annotation_Name = sapply(strsplit(Annotation, "_"), `[`, 2),
           Annotation_Sex = sapply(strsplit(Annotation, "_"), `[`, 3),
           Annotation_Number = as.integer(sapply(strsplit(Annotation, "_"), `[`, 1)),
           SIG = Pvalue <= minP.pvalue.thresh)

# pdf(file = paste0(dir.pdf, "mirQTL_pval_chromHMM_annotation_mirQTL-MAFLD_and_1000G-MAFLD_3thresholds.pdf"), width = 10, height = 4)
#
# # mirQTL MIXED
# df.mix.results %>%
#     filter(SE <= 10) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
#     geom_point(position = position_dodge(0.7), shape = 15, size = 1.5) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_alpha_manual(values = c(0.3,1)) +
#     scale_color_manual(values = cbPalette) +
#     geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(~PThresh) +
#     labs(y = "Log Odds Ratio (95% Conf. Int.)",
#          x = "ChromHMM 15-state Mnemonic",
#          color = "Sex",
#          alpha = "Significant Enrichment P-value",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "at 3 mirQTL p-value thresholds",
#          caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")
#
# df.mix.results %>%
#     filter(SE <= 10) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
#     geom_col(position = "dodge") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_fill_manual(values = cbPalette) +
#     geom_hline(aes(yintercept = -log10(mix.pvalue.thresh)), linetype = "dashed") +
#     facet_wrap(~PThresh) +
#     labs(y = "-Log10(P-value)",
#          x = "ChromHMM 15-state Mnemonic",
#          fill = "Sex",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "at 3 mirQTL p-value thresholds",
#          caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")
#
# # mirQTL MIXED minP
# df.minP.results %>%
#     filter(SE <= 10) %>%
#     filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
#     geom_point(position = position_dodge(0.7), shape = 15, size = 3) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_alpha_manual(values = c(0.3,1)) +
#     scale_color_manual(values = cbPalette) +
#     geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(~PThresh, nrow = 2) +
#     labs(y = "Log Odds Ratio (95% Conf. Int.)",
#          x = "ChromHMM 15-state Mnemonic",
#          color = "Sex",
#          alpha = "Significant Enrichment P-value",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "(minimum eQTL P-value at each variant)",
#          caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")
#
# ggsave("~/Desktop/gar_log_odds.pdf", width = 10, height = 7)
#
# df.minP.results %>%
#     filter(SE <= 10) %>%
#     filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
#     geom_col(position = "dodge") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_fill_manual(values = cbPalette) +
#     geom_hline(aes(yintercept = -log10(minP.pvalue.thresh)), linetype = "dashed") +
#     facet_wrap(~PThresh, nrow = 2) +
#     labs(y = "-Log10(P-value)",
#          x = "ChromHMM 15-state Mnemonic",
#          fill = "Sex",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "(minimum eQTL P-value at each variant)",
#          caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")
#
# ggsave("~/Desktop/gar_pval.pdf", width = 10, height = 7)
#
# # 1000Genomes EUR
# df.eur.results %>%
#     filter(SE <= 10) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
#     geom_point(position = position_dodge(0.7), shape = 15, size = 1.5) +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_alpha_manual(values = c(0.3,1)) +
#     scale_color_manual(values = cbPalette) +
#     geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
#     geom_hline(yintercept = 0) +
#     facet_wrap(~PThresh) +
#     labs(y = "Log Odds Ratio (95% Conf. Int.)",
#          x = "ChromHMM 15-state Mnemonic",
#          color = "Sex",
#          alpha = "Significant Enrichment P-value",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "at 3 mirQTL p-value thresholds",
#          caption = "Annotations with high error (SE>10) removed. Using 1000Genomes (EUR ancestry) MAF and LD.")
#
# df.eur.results %>%
#     filter(SE <= 10) %>%
#     ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
#     geom_col(position = "dodge") +
#     theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
#           legend.position = "bottom") +
#     scale_fill_manual(values = cbPalette) +
#     geom_hline(aes(yintercept = -log10(eur.pvalue.thresh)), linetype = "dashed") +
#     facet_wrap(~PThresh) +
#     labs(y = "-Log10(P-value)",
#          x = "ChromHMM 15-state Mnemonic",
#          fill = "Sex",
#          title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
#          subtitle = "at 3 mirQTL p-value thresholds",
#          caption = "Annotations with high error (SE>10) removed. Using 1000Genomes (EUR ancestry) MAF and LD.")
#
# dev.off()




#ggsave(paste0(dir.pdf, "chrom_states_or.pdf"), width = 10, height = 6)

#pdf(file = paste0(dir.pdf, "garfield_mirQTLcond_pval_chromHMM_annotation_mirQTL-MAFLD.pdf"), width = 10, height = 4)

# mirQTL MIXED minP
df.minP.results %>%
    #filter(SE <= 10) %>%
    #filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
    filter(PThresh == 1.434e-6) %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(0.7), shape = 15, size = 2) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = cbPalette) +
    geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2) +
    geom_hline(yintercept = 0) +
    #facet_wrap(~PThresh, ncol = 1) +
    labs(y = "Log Odds Ratio (95% Conf. Int.)",
         x = "ChromHMM 15-state Mnemonic",
         color = "Sex",
         alpha = "Significant Enrichment P-value",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         subtitle = "(minimum eQTL P-value at each variant)",
         caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")

#ggsave("~/Desktop/enrichment2.pdf", height = 5, width = 7)

df.minP.results %>%
    #filter(SE <= 10) %>%
    filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = cbPalette) +
    geom_hline(aes(yintercept = -log10(minP.pvalue.thresh)), linetype = "dashed") +
    facet_wrap(~PThresh, ncol = 1) +
    labs(y = "-Log10(P-value)",
         x = "ChromHMM 15-state Mnemonic",
         fill = "Sex",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         subtitle = "(minimum eQTL P-value at each variant)",
         caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")

#dev.off()

df.minP.results %>%
    #filter(SE <= 10) %>%
    #filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
    filter(PThresh == 1.434e-6) %>%
    mutate(Beta = if_else(SE > 10, 0, Beta)) %>%
    ggplot(aes(y = reorder(Annotation_Name, -Annotation_Number), x = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(width = 0.7), shape = 15, size = 1) +
    geom_errorbar(aes(xmin=CI95_lower, xmax=CI95_upper), position = position_dodge(0.7), width = 0.5) +
    geom_vline(xintercept = 0) +
    #facet_wrap(~PThresh) +
    coord_cartesian(xlim = c(-2.1, 2.1), clip = "off") +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = paperCatEight) +
    plotTheme("figure")

ggsave("~/Desktop/enrichment.pdf", height = 2.5, width = 3.5)

df.minP.results %>%
    #filter(SE <= 10) %>%
    #filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
    filter(PThresh == 1.434e-6) %>%
    #mutate(Beta = if_else(SE > 10, 0, Beta)) %>%
    ggplot(aes(y = reorder(Annotation_Name, -Annotation_Number), x = -log10(Pvalue), fill = Annotation_Sex)) +
    geom_col(position = position_dodge()) +
    geom_vline(xintercept = -log10(minP.pvalue.thresh), linetype = "dashed") +
    scale_fill_manual(values = paperCatEight) +
    scale_x_continuous(expand = expansion(mult = c(0,.1))) +
    plotTheme("figure")

ggsave("~/Desktop/pvals.pdf", height = 2.5, width = 3.5)

df.minP.results %>%
    #filter(SE <= 10) %>%
    #filter(PThresh == 1.434e-6 | PThresh == 2.804e-5) %>%
    #filter(PThresh == 1.434e-6) %>%
    mutate(Beta = if_else(SE > 10, 0, Beta)) %>%
    ggplot(aes(y = reorder(Annotation_Name, -Annotation_Number), x = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(width = 0.7), shape = 15, size = 1) +
    geom_errorbar(aes(xmin=CI95_lower, xmax=CI95_upper), position = position_dodge(0.7), width = 0.5) +
    geom_vline(xintercept = 0) +
    facet_wrap(~PThresh) +
    coord_cartesian(xlim = c(-2.1, 2.1), clip = "off") +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = paperCatEight) +
    plotTheme("figure")

ggsave("~/Desktop/enrichment_thresholds.pdf", height = 7, width = 7)
