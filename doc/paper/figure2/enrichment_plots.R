
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
#library(DESeq2)
#library(RColorBrewer)
#library(qqman)
library(mikelaffr)

# FIGURE ###############################################################################################################
# figure 2B: miRNA-eQTL enrichment in chrom annotations
output.odds.pdf <- paste0(here("doc/paper/figure2/pdfs/"), "figure2B_enrichment_odds.pdf")
output.pval.pdf <- paste0(here("doc/paper/figure2/pdfs/"), "figure2B_enrichment_pval.pdf")


# OUTPUT FILES #########################################################################################################
# output directory for pdf files
#dir.pdf <- here("doc/paper/figure2/pdfs/")


# INPUT FILES ##########################################################################################################
# mirQTL enrichment on chromHMM using mirQTL-MIXED maf/ld and minimum P-value
minP.results.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL-minP_chromHMM/garfield.test.hg38-mirQTL-MIXED_mirQTL-minP_chromHMM.out")
# p-value file
minP.pvalue.txt <- here("results/garfield/hg38-mirQTL-MIXED_mirQTL-minP_chromHMM/garfield.Meff.hg38-mirQTL-MIXED_mirQTL-minP_chromHMM.out")

# GLOBALS ##############################################################################################################


# Import Results #######################################################################################################

df.minP.results <- read_table2(minP.results.txt)
minP.pvalue.thresh <- as.numeric(strsplit(read_lines(minP.pvalue.txt), "\t")[[2]][2])

df.minP.results %<>%
    mutate(Annotation_Name = sapply(strsplit(Annotation, "_"), `[`, 2),
           Annotation_Sex = sapply(strsplit(Annotation, "_"), `[`, 3),
           Annotation_Number = as.integer(sapply(strsplit(Annotation, "_"), `[`, 1)),
           SIG = Pvalue <= minP.pvalue.thresh)

df.minP.results$Beta[df.minP.results$SE > 10] <- NA
df.minP.results$CI95_lower[df.minP.results$SE > 10] <- NA
df.minP.results$CI95_upper[df.minP.results$SE > 10] <- NA
df.minP.results$Pvalue[df.minP.results$SE > 10] <- NA

df.minP.results %>%
    filter(PThresh == 1.434e-6) %>%
    ggplot(mapping = aes(x = reorder(Annotation_Name, Annotation_Number), y = Beta, color = Annotation_Sex, alpha = SIG)) +
    geom_point(position = position_dodge(0.7), shape = 15, size = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom",
          legend.title = element_text(size = 4),
          legend.text = element_text(size = 2),
          axis.text = element_text(size = 4),
          axis.title = element_text(size = 4),
          plot.caption = element_text(size = 4)) +
    scale_alpha_manual(values = c(0.3,1)) +
    scale_color_manual(values = c("green", "navy")) +
    geom_errorbar(aes(ymin=CI95_lower, ymax=CI95_upper), position = position_dodge(0.7), width = 0.2, size = 0.5) +
    geom_hline(yintercept = 0) +
    labs(y = "Log Odds Ratio (95% Conf. Int.)",
         x = "ChromHMM 15-state Mnemonic",
         color = "Sex",
         alpha = "Significant Enrichment P-value",
         #title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         #subtitle = "(minimum eQTL P-value at each variant)",
         caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")

ggsave(filename = output.odds.pdf, width = 6.5, height = 2, units = "in")

df.minP.results %>%
    filter(PThresh == 1.434e-6) %>%
    ggplot(aes(x = reorder(Annotation_Name, Annotation_Number), y = -log10(Pvalue), fill = Annotation_Sex)) +
    geom_col(position = "dodge") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
          legend.position = "bottom") +
    scale_fill_manual(values = c("green", "navy")) +
    geom_hline(aes(yintercept = -log10(minP.pvalue.thresh)), linetype = "dashed") +
    labs(y = "-Log10(P-value)",
         x = "ChromHMM 15-state Mnemonic",
         fill = "Sex",
         title = "GARFIELD: mirQTL Enrichment within Fetal Brain Chromatin States",
         subtitle = "(minimum eQTL P-value at each variant)",
         caption = "Annotations with high error (SE>10) removed. Using mirQTL (MIXED ancestry) MAF and LD.")

ggsave(filename = output.pval.pdf, width = 6.5, height = 2, units = "in")

theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.spacing.x = unit(1, "mm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 10)) +


