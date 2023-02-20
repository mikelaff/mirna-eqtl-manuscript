# Summarized mlma association results and compare to emmax results

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
#library(DESeq2)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mlma/pdfs/")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis from mlma
eqtls.dataframe.mlma.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor.eQTLs.dataFrame.rds")
eqtls.dataframe.mlma.npc.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor.eQTLs.dataFrame.rds")
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis from emmax
eqtls.dataframe.emmax.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor.eQTLs.dataFrame.rds")
eqtls.dataframe.emmax.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA from mlma
summarized.results.dataframe.mlma.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor_dataFrame.rds")
summarized.results.dataframe.mlma.npc.rds <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor_dataFrame.rds")
# Summarized association results for every variant within 1MB of each expressed miRNA from emmax
summarized.results.dataframe.emmax.fp.rds <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_dataFrame.rds")
summarized.results.dataframe.emmax.pr.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_dataFrame.rds")

# nominal p-value threshold from mlma
nom.p.val.mlma.txt <- here("results/mlma/association_results/20200114_mirQTLor_mlma/compiled/20200114_mirQTLor_nomPvalue.txt")
nom.p.val.mlma.npc.txt <- here("results/mlma/association_results/20200114_mirQTLor_mlma_npc/compiled/20200114_mirQTLor_nomPvalue.txt")
# nominal p-value threshold from emmax
nom.p.val.emmax.fp.txt <- here("results/emmax/association_results/20191204_mirQTLor/compiled/20191204_mirQTLor_nomPvalue.txt")
nom.p.val.emmax.pr.txt <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_nomPvalue.txt")

# GLOBALS ##############################################################################################################


# Import Summarized Results ############################################################################################
# compiled eqtls
df.eqtls.emmax.fp <- as_tibble(readRDS(eqtls.dataframe.emmax.fp.rds))
df.eqtls.emmax.fp %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+")) %>%
    filter(!grepl("chrX", esnp))

df.eqtls.emmax.pr <- as_tibble(readRDS(eqtls.dataframe.emmax.pr.rds))
df.eqtls.emmax.pr %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+")) %>%
    filter(!grepl("chrX", esnp))

df.eqtls.mlma <- as_tibble(readRDS(eqtls.dataframe.mlma.rds))
df.eqtls.mlma %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+")) %>%
    filter(!grepl("chrX", esnp))

df.eqtls.mlma.npc <- as_tibble(readRDS(eqtls.dataframe.mlma.npc.rds))
df.eqtls.mlma.npc %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+")) %>%
    filter(!grepl("chrX", esnp))

# summarized results
df.results.emmax.fp <- as_tibble(readRDS(summarized.results.dataframe.emmax.fp.rds))
df.results.emmax.fp %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.emmax.pr <- as_tibble(readRDS(summarized.results.dataframe.emmax.pr.rds))
df.results.emmax.pr %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.mlma <- as_tibble(readRDS(summarized.results.dataframe.mlma.rds))
df.results.mlma %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

df.results.mlma.npc <- as_tibble(readRDS(summarized.results.dataframe.mlma.npc.rds))
df.results.mlma.npc %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# nom p-values
nom.p.val.emmax.fp <- as.numeric(read_lines(nom.p.val.emmax.fp.txt))
nom.p.val.emmax.pr <- as.numeric(read_lines(nom.p.val.emmax.pr.txt))
nom.p.val.mlma <- as.numeric(read_lines(nom.p.val.mlma.txt))
nom.p.val.mlma.npc <- as.numeric(read_lines(nom.p.val.mlma.npc.txt))
# P-value Comparison ###################################################################################################

tmp <- left_join(df.results.emmax.pr, df.results.mlma, by = "UniName_SNP", suffix = c(".EMMAX.PR", ".MLMA"))

tmp %>%
    filter(UniName_SNP %in% df.eqtls.emmax.pr$eqtl | UniName_SNP %in% df.eqtls.mlma$eqtl) -> tmp2

tmp2 %<>%
    filter(!is.na(P.MLMA), CHR.EMMAX.PR != "chrX")

tmp2$sig.emmax.pr <- tmp2$P.EMMAX.PR <= nom.p.val.emmax.pr
tmp2$sig.mlma <- tmp2$P.MLMA <= nom.p.val.mlma

tmp2$sig <- "neither"
tmp2$sig <- ifelse(tmp2$sig.emmax.pr & tmp2$sig.mlma, "both", ifelse(tmp2$sig.emmax.pr, "emmax.pr", ifelse(tmp2$sig.mlma, "mlma", "none")))

tmp2$sig <- factor(tmp2$sig)

tmp2 %>%
    ggplot(aes(x = -log10(P.EMMAX.PR), y = -log10(P.MLMA), color = sig)) +
    geom_point(size = 1) +
    geom_hline(yintercept = -log10(nom.p.val.mlma), linetype = "dashed") +
    geom_vline(xintercept = -log10(nom.p.val.emmax.pr), linetype = "dashed") +
    scale_color_manual(values = cbPalette[c(1,2,4)])

ggsave(paste0(dir.pdfs, "p-vales_mlma_v_emmax_pr.pdf"))

tmp2 %>%
    filter(sig == "mlma") -> tmp3

# Venn Diagrams ########################################################################################################
sum(df.eqtls.emmax.fp$eqtl %in% df.eqtls.emmax.pr$eqtl)

sum(df.eqtls.mlma$eqtl %in% df.eqtls.mlma.npc$eqtl)

venn.diagram(x = list(df.eqtls.emmax.fp$eqtl, df.eqtls.emmax.pr$eqtl),
             category.names = c("EMMAX", "EMMAX.PR"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emmax.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             #cat.pos = c(-105, 70),
             cat.dist = c(0.05, 0.07),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.mlma$eqtl, df.eqtls.mlma.npc$eqtl),
             category.names = c("MLMA", "MLMA.NPC"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_mlma.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-20, 160),
             cat.dist = c(0.02, 0.01),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.emmax.pr$eqtl, df.eqtls.mlma.npc$eqtl),
             category.names = c("EMMAX.PR", "MLMA.NPC"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emmapr_mlmanpc.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-10, 160),
             cat.dist = c(0.02, 0.01),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.emmax.pr$eqtl, df.eqtls.mlma$eqtl),
             category.names = c("EMMAX.PR", "MLMA"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emmapr_mlma.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-30, 150),
             cat.dist = c(0.04, 0.05),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.emmax.fp$eqtl, df.eqtls.mlma.npc$eqtl),
             category.names = c("EMMAX", "MLMA.NPC"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emma_mlmanpc.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-30, 160),
             cat.dist = c(0.02, 0.01),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.emmax.fp$eqtl, df.eqtls.mlma$eqtl),
             category.names = c("EMMAX", "MLMA"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emma_mlma.png"),
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
             fill = brewer.pal(4, "Dark2")[1:2],

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             cat.pos = c(-30, 160),
             cat.dist = c(0.04, 0.04),
             cat.fontfamily = "sans",
             cat.just = list(c(0,0), c(1,0)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
)

venn.diagram(x = list(df.eqtls.emmax.fp$eqtl, df.eqtls.emmax.pr$eqtl, df.eqtls.mlma$eqtl, df.eqtls.mlma.npc$eqtl),
             category.names = c("EMMAX", "EMMAX.PR", "MLMA", "MLMA.NPC"),
             main = "eQTLs",
             filename = paste0(dir.pdfs, "venn_eqtls_emmax_v_mlma.png"),
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
             fill = brewer.pal(4, "Dark2"),

             # Numbers
             cex = .4,
             fontface = "bold",
             fontfamily = "sans",

             # Set names
             cat.cex = 0.35,
             cat.fontface = "bold",
             cat.default.pos = "outer",
             #cat.pos = c(-105, 70),
             cat.dist = c(0.2, 0.2, 0.1, 0.1),
             cat.fontfamily = "sans",
             cat.just = list(c(0.4,0), c(0.65,0), c(0.5,0.5), c(0.5,0.5)),

             # Main
             main.fontfamily = "sans",
             main.fontface = "bold",
             main.pos = c(.5, 1)
             )


















