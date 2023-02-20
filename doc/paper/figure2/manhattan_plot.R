
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(RColorBrewer)
#library(qqman)
library(mikelaffr)

# FIGURE ###############################################################################################################
# figure 2A: manhattan plot of all miRNA-eQTLs
output.pdf <- paste0(here("doc/paper/figure2/pdfs/"), "figure2A_manhattan.pdf")
output.png <- paste0(here("doc/paper/figure2/pdfs/"), "figure2A_manhattan.png")
output_labeled.pdf <- paste0(here("doc/paper/figure2/pdfs/"), "figure2A_manhattan_labeled_max.pdf")


# OUTPUT FILES #########################################################################################################
# output directory for pdf files
#dir.pdf <- here("doc/paper/figure2/pdfs/")


# INPUT FILES ##########################################################################################################
# primary association results
primary.results.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")
# secondary association results
secondary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/secondary/association_results/compiled/20200120_mirQTLor_secondary_variants_dataFrame.rds")
# tertiary association results
tertiary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/tertiary/association_results/compiled/20200120_mirQTLor_tertiary_variants_dataFrame.rds")
# quarternary association results
# quarternary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/quarternary/association_results/compiled/20200120_mirQTLor_quarternary_variants_dataFrame.rds")

# primary eqtls
primary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/primary/20200120_mirQTLor_primary_eQTLs_dataFrame.rds")
# secondary eqtls
secondary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/secondary/20200120_mirQTLor_secondary_eQTLs_dataFrame.rds")
# tertiary eqtls
tertiary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/tertiary/20200120_mirQTLor_tertiary_eQTLs_dataFrame.rds")

# nominal p-value from eigenMT-BH procedure
nom.p.value.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Import Summarized Results ############################################################################################
# summarized eQTLs
df.primary.eqtls <- read_rds(primary.eqtls.rds)
df.secondary.eqtls <- read_rds(secondary.eqtls.rds)
df.tertiary.eqtls <- read_rds(tertiary.eqtls.rds)

# label eQTLs
df.primary.eqtls$degree <- "primary"
df.secondary.eqtls$degree <- "secondary"
df.tertiary.eqtls$degree <- "tertiary"

# combine into one df
df.eqtls <- bind_rows(df.primary.eqtls,
                      df.secondary.eqtls,
                      df.tertiary.eqtls)
rm(df.primary.eqtls, df.secondary.eqtls, df.tertiary.eqtls)

# eqtl name for plotting
df.eqtls %<>%
    mutate(UniName_SNP_degree = paste(UniName, SNP, degree, sep = "+"))

df.eqtls$color <- NA
df.eqtls$color[df.eqtls$degree == "primary"] <- paperCatEight[2]
df.eqtls$color[df.eqtls$degree == "secondary"] <- paperCatEight[3]
df.eqtls$color[df.eqtls$degree == "tertiary"] <- paperCatEight[4]


# all variants
df.results.primary <- as_tibble(readRDS(primary.results.rds))
df.results.secondary <- as_tibble(readRDS(secondary.results.rds))
df.results.tertiary <- as_tibble(readRDS(tertiary.results.rds))

# label variants
df.results.primary$degree <- "primary"
df.results.secondary$degree <- "secondary"
df.results.tertiary$degree <- "tertiary"

# combine into one df
df.results <- bind_rows(df.results.primary,
                        df.results.secondary,
                        df.results.tertiary)
rm(df.results.primary, df.results.secondary, df.results.tertiary)

# eqtl name for plotting
df.results %<>%
    mutate(UniName_SNP_degree = paste(UniName, SNP, degree, sep = "+"))

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.value.txt))

# Manhattan Plot #######################################################################################################

# chr number for plotting
df.results$chrLabel <- factor(sapply(strsplit(df.results$CHR, "chr"), `[`, 2), levels = c(1:22,"X"), ordered = TRUE)

# # qqman for Manhattan
# df.results$chrNum <- sapply(strsplit(df.results$CHR, "chr"), `[`, 2)
# df.results$chrNum[df.results$chrNum == "X"] <- 23
# df.results$chrNum <- as.numeric(df.results$chrNum)
#
# manhattan(x = df.results,
#           chr = "chrNum",
#           bp = "BP.hg38",
#           p = "P",
#           chrlabs = "chrLabel")


# ggplot Manhattan

# dummy points to get proper chrom lengths
tmp <- tibble(CHR = names(CHROM.LENGTHS), BP.hg38 = CHROM.LENGTHS)
tmp$chrLabel <- factor(sapply(strsplit(tmp$CHR, "chr"), `[`, 2), levels = c(1:22,"X"), ordered = TRUE)
tmp$P <- 0.009

df.results %<>%
    bind_rows(tmp)
rm(tmp)

df.results %>%
    #filter(P < 0.00001) %>%
    ggplot() +
    geom_point(aes(x = BP.hg38, y = -log10(P), color = chrLabel), alpha = 1, size = 0.15) +
    facet_grid(~chrLabel, scales = "free_x", space = "free_x", switch = "x") +
    geom_hline(yintercept = -log10(nom.p.val), color = "black", linetype = "dashed") +
    labs(x = "Chromosome",
         y = expression(paste("-lo",g[10],"(",italic("p-value"),")"))) +
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
    scale_color_manual(values = c(rep(c(paperLightBlue,paperDarkBlue),12))) +
    scale_y_continuous(expand=expansion(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    coord_cartesian(clip = "off") +
    geom_point(mapping = aes(x = BP.hg38, y = -log10(P)),
               data = filter(df.results, UniName_SNP_degree %in% df.eqtls$UniName_SNP_degree),
               size = 1,
               color = df.eqtls$color,
               shape = 17) -> plot1

ggsave(filename = output.png, plot = plot1, height = 2.5, width = 6.5, units = "in", dpi = 600)

# PDF with labeled miRs

# eqtls for labeling
df.eqtls$label <- TRUE
df.eqtls$label[df.eqtls$UniName == "hsa-mir-4707_hsa-miR-4707-3p"] <- TRUE
df.eqtls$label[df.eqtls$UniName == "hsa-mir-5683_hsa-miR-5683" & df.eqtls$degree == "primary"] <- TRUE

df.results %>%
    filter(P < 0.0001) %>%
    ggplot() +
    geom_point(aes(x = BP.hg38, y = -log10(P), color = chrLabel), alpha = 1, size = 0.15) +
    facet_grid(~chrLabel, scales = "free_x", space = "free_x", switch = "x") +
    geom_hline(yintercept = -log10(nom.p.val), color = "black", linetype = "dashed") +
    labs(x = "Chromosome",
         y = expression(paste("-lo",g[10],"(",italic("p-value"),")"))) +
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
    scale_color_manual(values = c(rep(c(paperLightBlue,paperDarkBlue),12))) +
    scale_y_continuous(expand=expansion(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    coord_cartesian(clip = "off") +
    geom_point(mapping = aes(x = BP.hg38, y = -log10(P)),
               data = filter(df.results, UniName_SNP_degree %in% df.eqtls$UniName_SNP_degree),
               size = 1,
               color = df.eqtls$color,
               shape = 17) +
    geom_label(data = filter(df.results, UniName_SNP_degree %in% filter(df.eqtls, label)$UniName_SNP_degree),
               mapping = aes(x = BP.hg38, y = -log10(P), label = UniName),
               hjust = "right",
               vjust = "bottom",
               size = 2)

ggsave(filename = output_labeled.pdf, height = 10, width = 20, units = "in")
