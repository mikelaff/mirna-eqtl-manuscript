# Dot plots and Geno x Pheno plots for mirQTLor associations

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_plots/pdfs/")
dir.pngs <- here("doc/mirQTL_plots/pngs/")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")


# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################

# compiled eQTLs
df.eqtls <- as_tibble(readRDS(eqtls.dataframe.rds))

# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# # expression data
# vsd <- readRDS(vsd.rds)

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# Manhattan Plot #######################################################################################################



# chr number for plotting
df.results$chrLabel <- factor(sapply(strsplit(df.results$CHR, "chr"), `[`, 2), levels = c(1:22,"X"), ordered = TRUE)

# dummy points to get proper chrom lengths
tmp <- data.frame(chr = names(CHROM.LENGTHS), BP.hg38 = CHROM.LENGTHS, stringsAsFactors = FALSE)
tmp$chrLabel <- factor(sapply(strsplit(tmp$chr, "chr"), `[`, 2), levels = c(1:22,"X"), ordered = TRUE)
tmp$P <- 0.09


df.results %>%
    bind_rows(tmp) %>%
    filter(P < 0.1) %>%
    ggplot() +
    geom_point(aes(x = BP.hg38, y = -log10(P), color = chrLabel), alpha = 0.7, size = 0.2) +
    facet_grid(~chrLabel, scales = "free_x", space = "free_x", switch = "x") +
    geom_hline(yintercept = -log10(nom.p.val), color = "black", linetype = "solid") +
    labs(x = "Chromosome",
         y = expression(paste("-lo",g[10],"(",italic("p-value"),")"))) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0, "mm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 10)) +
    scale_color_manual(values = c(rep(c("navy","cornflowerblue"),12))) +
    scale_y_continuous(expand=expansion(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(mapping = aes(x = BP.hg38, y = -log10(P)), data = filter(df.results, UniName_SNP %in% df.eqtls$eQTL), size=1, color="red")

ggsave(paste0(dir.pngs,"eqtl_manhattan_OLD-bh_rast.png"), height = 5, width = 10, units = "in", dpi = "retina")

