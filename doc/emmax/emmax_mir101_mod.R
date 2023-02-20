# Summarized EMMAX association results from the final pipeline

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.rds <- here("results/emmax/association_results/20200102_mirQTLor_mir101_mod/compiled/20200102_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200102_mirQTLor_mir101_mod/compiled/20200102_mirQTLor_dataFrame.rds")

# nominal p-value threshold
nom.p.val.txt <- here("results/emmax/association_results/20200102_mirQTLor_mir101_mod/compiled/20200102_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200102_mirQTLor_mir101_mod/20200102_mirQTLor_mir101_mod.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################
# compiled eqtls
df.eqtls <- as_tibble(readRDS(eqtls.dataframe.rds))
df.eqtls %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# nominal p-value
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# expression data
vsd <- readRDS(vsd.rds)

# mirna granges
gr.mirs <- readRDS(gr.mirna.rds)

# Find Overlaps ########################################################################################################
df.esnps <- df.eqtls

df.esnps <- left_join(df.esnps, df.results, by = c("esnp" = "SNP"))
df.esnps %<>%
    filter(eqtl == UniName_SNP)

df.esnps %<>%
    mutate(seqnames = CHR, start = BP.hg38, end = BP.hg38)

gr.esnps <- makeGRangesFromDataFrame(df.esnps, keep.extra.columns = TRUE, seqnames.field = "seqnames", ignore.strand = TRUE)

gr.emirs <- gr.mirs[gr.mirs$uniqueName %in% df.eqtls$emir]

findOverlaps(gr.esnps, gr.emirs)

# Manhattan Plot #######################################################################################################


#check the male dosage genotypes on the x
#large LD blocks on chr6 and chr17
# plot all eQTLs dots and geno X pheno
# look into kinship for JC

df.results <- as_tibble(df.results)


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
          strip.text = element_text(size = 12)) +
    scale_color_manual(values = c(rep(c("navy","cornflowerblue"),12))) +
    scale_y_continuous(expand=expand_scale(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(mapping = aes(x = BP.hg38, y = -log10(P)), data = filter(df.results, UniName_SNP %in% df.eqtls$eqtl), size=1, color="red")

#ggsave("~/Desktop/eqtl_man_rast.png", height = 4, width = 10, units = "in", dpi = "retina")




















