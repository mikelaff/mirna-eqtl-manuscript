# Summarized EMMAX association results comparison of pregression to final pipeline

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(DESeq2)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor.eQTLs.dataFrame.rds")
eqtls.dataframe.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_dataFrame.rds")
summarized.results.dataframe.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_dataFrame.rds")

# nominal p-value threshold
nom.p.val.fp.txt <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_nomPvalue.txt")
nom.p.val.pr.txt <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.fp.rds <- here("results/emmax/phenotype_files/20191204_mirQTLor_VST_miRNA_expression/20191204_mirQTLor_VST_miRNA_expression_rse.rds")
vsd.pr.rds <- here("results/emmax/phenotype_files/20200103_mirQTLor_VST_miRNA_expression_preregress/20200103_mirQTLor_VST_miRNA_expression_preregress_rse.rds")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################
# compiled eqtls
df.eqtls.fp <- as_tibble(readRDS(eqtls.dataframe.fp.rds))
df.eqtls.fp %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

df.eqtls.pr <- as_tibble(readRDS(eqtls.dataframe.pr.rds))
df.eqtls.pr %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

# all variants
df.results.fp <- as_tibble(readRDS(summarized.results.dataframe.fp.rds))
df.results.fp %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.pr <- as_tibble(readRDS(summarized.results.dataframe.pr.rds))
df.results.pr %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# nominal p-value
nom.p.val.fp <- as.numeric(read_lines(nom.p.val.fp.txt))
nom.p.val.pr <- as.numeric(read_lines(nom.p.val.pr.txt))

# expression data
vsd.fp <- readRDS(vsd.fp.rds)
vsd.pr <- readRDS(vsd.pr.rds)

# mirna granges
gr.mirs <- readRDS(gr.mirna.rds)

# Find Overlaps ########################################################################################################
# df.esnps <- df.eqtls
#
# df.esnps <- left_join(df.esnps, df.results, by = c("esnp" = "SNP"))
# df.esnps %<>%
#     filter(eqtl == UniName_SNP)
#
# df.esnps %<>%
#     mutate(seqnames = CHR, start = BP.hg38, end = BP.hg38)
#
# gr.esnps <- makeGRangesFromDataFrame(df.esnps, keep.extra.columns = TRUE, seqnames.field = "seqnames", ignore.strand = TRUE)
#
# gr.emirs <- gr.mirs[gr.mirs$uniqueName %in% df.eqtls$emir]
#
# findOverlaps(gr.esnps, gr.emirs)

# Compare Datasets #####################################################################################################
sum(df.eqtls.fp$eqtl %in% df.eqtls.pr$eqtl)

sum(df.eqtls.pr$eqtl %in% df.eqtls.fp$eqtl)

venn.diagram(x = list(df.eqtls.fp$eqtl, df.eqtls.pr$eqtl),
             category.names = c("EMMAX Cov.", "Pre-regress Cov."),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emmax_v_preregress.png"),
             output = TRUE,

             # Output features
             imagetype = "png",
             height = 480,
             width = 480,
             resolution = 300,
             compression = "lzw",

             # Circles
             lwd = 3,
             lty = 'blank',
             fill = brewer.pal(3, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.4,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-105, 70),
             cat.dist = c(0.1, 0.1),
             cat.fontfamily = "sans",
             cat.just = list(c(-.5,1), c(1,1)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, .95)
             )

venn.diagram(x = list(df.eqtls.fp$emir, df.eqtls.pr$emir),
             category.names = c("EMMAX Cov.", "Pre-regress Cov."),
             main = "emiRs",
             filename = paste0(dir.pdfs, "venn_emirs_emmax_v_preregress.png"),
             output = TRUE,

             # Output features
             imagetype = "png",
             height = 480,
             width = 480,
             resolution = 300,
             compression = "lzw",

             # Circles
             lwd = 3,
             lty = 'blank',
             fill = brewer.pal(3, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.4,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-105, 70),
             cat.dist = c(0.1, 0.1),
             cat.fontfamily = "sans",
             cat.just = list(c(-.5,1), c(1,1)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, .95)
)

venn.diagram(x = list(df.eqtls.fp$esnp, df.eqtls.pr$esnp),
             category.names = c("EMMAX Cov.", "Pre-regress Cov."),
             main = "eSNPs",
             filename = paste0(dir.pdfs, "venn_esnps_emmax_v_preregress.png"),
             output = TRUE,

             # Output features
             imagetype = "png",
             height = 480,
             width = 480,
             resolution = 300,
             compression = "lzw",

             # Circles
             lwd = 3,
             lty = 'blank',
             fill = brewer.pal(3, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.4,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-105, 70),
             cat.dist = c(0.1, 0.1),
             cat.fontfamily = "sans",
             cat.just = list(c(-.5,1), c(1,1)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, .95)
)

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
          strip.text = element_text(size = 12)) +
    scale_color_manual(values = c(rep(c("navy","cornflowerblue"),12))) +
    scale_y_continuous(expand=expand_scale(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(mapping = aes(x = BP.hg38, y = -log10(P)), data = filter(df.results, SNP_UniName %in% df.eqtls$eqtl), size=1, color="red")

ggsave("~/Desktop/eqtl_man_rast.png", height = 4, width = 10, units = "in", dpi = "retina")




















