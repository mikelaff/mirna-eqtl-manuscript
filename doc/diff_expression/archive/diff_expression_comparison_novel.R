# diff expression gz/cp

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(readxl)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/diff_expression/pdfs/")
# write plots
write_plots <- FALSE
# combined deseq2 shrunken results as a dataframe with labels for significance and category (neuro vs maturation)
combined_df_rds <- here("doc/diff_expression/rdata/combined_deseq2_shrunken_results_novel_df.rds")


# INPUT FILES #####################################################################################
# rds with a DESeqDataSet with completed diff. expression analysis on gz/cp mirna expression
dds_gz_cp_rds <- here("doc/diff_expression/rdata/dds_gz_cp_novel.rds")
# rds with a DESeqDataSet with completed diff. expression analysis on gw mirna expression
dds_gw_rds <- here("doc/diff_expression/rdata/dds_gw_novel.rds")
# mirtarbase v7 hsa targets
mirna_targets_xlsx <- here("data/mirtarbase_v7/hsa_MTI.xlsx")

# Load DDS #######################################################################################
dds.gzcp <- readRDS(dds_gz_cp_rds)
dds.gw <- readRDS(dds_gw_rds)

# Results #########################################################################################
res.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")
res.gw <- lfcShrink(dds.gw, coef="gestation_week", type="apeglm")

# data frame from results
res.gzcp.df <- as.data.frame(res.gzcp)
res.gw.df <- as.data.frame(res.gw)

# create mirna column
res.gzcp.df$mirna <- rownames(res.gzcp)
res.gw.df$mirna <- rownames(res.gw)

# create column to specify significance
res.gzcp.df$sig <- FALSE
res.gw.df$sig <- FALSE

# label mirnas that have sig. diff expression
padj_threshold <- 0.1
res.gzcp.df$sig[which(res.gzcp.df$padj < padj_threshold)] <- TRUE
res.gw.df$sig[which(res.gw.df$padj < padj_threshold)] <- TRUE

# order data frame for plotting
res.gzcp.df <- res.gzcp.df[order(res.gzcp.df$padj, decreasing = TRUE),]
res.gw.df <- res.gw.df[order(res.gw.df$padj, decreasing = TRUE),]

# MA Plot #########################################################################################
circColor <- "black"
labelColor <- "black"

res.gzcp.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression GZ/CP") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey50", "red")) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(res.gzcp.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gzcp.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

res.gw.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point() +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across Gestation Week") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey50", "red")) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-124-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-124-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-9-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-9-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-92b-3p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(res.gw.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange, label=mirna), hjust=1, vjust=1.5, color=labelColor) +
  geom_point(data=subset(res.gw.df, mirna == "hsa-miR-92b-5p"), aes(baseMean, log2FoldChange), shape=21, size=4, stroke=2, color=circColor)

#if(write_plots){ggsave(file.path(pdf_dir,"diff_expression_gw_shrunklfc_labels.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Compare to GZ/CP Diff. Expression ###############################################################
combined_shrunk_results_df <- full_join(res.gw.df, res.gzcp.df, suffix = c(".gw", ".gzcp"), by = "mirna")

# after full join, some are NA if they arent in both tables, label these as not significant
combined_shrunk_results_df$sig.gw[which(is.na(combined_shrunk_results_df$sig.gw))] <- FALSE
combined_shrunk_results_df$sig.gzcp[which(is.na(combined_shrunk_results_df$sig.gzcp))] <- FALSE

# create column for significant in both analyses, none should be "pend" after proper labeling
combined_shrunk_results_df$sig.both <- "pend"

# label significant in both analyses column
combined_shrunk_results_df$sig.both[which(combined_shrunk_results_df$sig.gw & combined_shrunk_results_df$sig.gzcp)] <- "both"
combined_shrunk_results_df$sig.both[which(combined_shrunk_results_df$sig.gw & !combined_shrunk_results_df$sig.gzcp)] <- "gw only"
combined_shrunk_results_df$sig.both[which(!combined_shrunk_results_df$sig.gw & combined_shrunk_results_df$sig.gzcp)] <- "gzcp only"
combined_shrunk_results_df$sig.both[which(!combined_shrunk_results_df$sig.gw & !combined_shrunk_results_df$sig.gzcp)] <- "neither"
# factor sig.both variable
combined_shrunk_results_df$sig.both <- factor(combined_shrunk_results_df$sig.both)
# set log2fold change of NA values to 0
combined_shrunk_results_df$log2FoldChange.gw[which(is.na(combined_shrunk_results_df$log2FoldChange.gw))] <- 0
combined_shrunk_results_df$log2FoldChange.gzcp[which(is.na(combined_shrunk_results_df$log2FoldChange.gzcp))] <- 0

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gw)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gzcp)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

# order rows for plotting
combined_shrunk_results_df <- combined_shrunk_results_df[order(combined_shrunk_results_df$sig.both, decreasing = TRUE),]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=1) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNA Differential Expression",
       color="Significant") +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
  #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=1, color=circColor) +
  #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=1, color=circColor)
  #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
  #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

ggsave("~/Desktop/combined_de.pdf", width = 7, height = 6, useDingbats = FALSE)

if(write_plots){ggsave(file.path(pdf_dir,"diff_expression_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Subset miRNAs into Groups #######################################################################
# mirnas that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
miRNAs_maturation_late <- combined_shrunk_results_df$mirna[which(combined_shrunk_results_df$sig.both == "gw only" & combined_shrunk_results_df$log2FoldChange.gw > 0)]
# mirnas that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
miRNAs_maturation_early <- combined_shrunk_results_df$mirna[which(combined_shrunk_results_df$sig.both == "gw only" & combined_shrunk_results_df$log2FoldChange.gw < 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
miRNAs_neurogenesis_gz <- combined_shrunk_results_df$mirna[which((combined_shrunk_results_df$sig.both == "both" | combined_shrunk_results_df$sig.both == "gzcp only") & combined_shrunk_results_df$log2FoldChange.gzcp > 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
miRNAs_neurogenesis_cp <- combined_shrunk_results_df$mirna[which((combined_shrunk_results_df$sig.both == "both" | combined_shrunk_results_df$sig.both == "gzcp only") & combined_shrunk_results_df$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
combined_shrunk_results_df$category <- "none"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$mirna %in% miRNAs_maturation_late)] <- "maturation_late"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$mirna %in% miRNAs_maturation_early)] <- "maturation_early"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$mirna %in% miRNAs_neurogenesis_gz)] <- "neurogenesis_gz"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$mirna %in% miRNAs_neurogenesis_cp)] <- "neurogenesis_cp"
# factor category variable
combined_shrunk_results_df$category <- factor(combined_shrunk_results_df$category)

# plot the categories
combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=1) +
  plotTheme +
    labs(y="Log2 Fold Change GZ/CP",
         x="Log2 Fold Change Gestation Week",
         title="miRNA Differential Expression",
         color="Significant")+
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
    geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
    #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
    #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
    #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=1.5, color=labelColor) +
    geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=1, color=circColor) +
    #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-124-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
    #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
    #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-9-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
    geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=1, vjust=-0.5, color=labelColor) +
    geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=1, color=circColor)
    #geom_label(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label=mirna), hjust=0, vjust=-0.5, color=labelColor) +
    #geom_point(data=subset(combined_shrunk_results_df, mirna == "hsa-miR-92b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

ggsave("~/Desktop/combined_de_cat.pdf", width = 7, height = 6, useDingbats = FALSE)

combined_shrunk_results_df2 <- left_join(combined_shrunk_results_df, rows, by = "mirna")
# plot the categories
combined_shrunk_results_df2 %>%
    ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
    geom_point(size=1) +
    plotTheme +
    labs(y="Log2 Fold Change GZ/CP",
         x="Log2 Fold Change Gestation Week",
         title="miRNA Differential Expression",
         color="Significant")+
    scale_color_manual(values=cbPalette[c(1,2,4,5,3)])  +
    geom_point(data=subset(combined_shrunk_results_df2, novel == "novel" & category != "none"), aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp), shape=21, size=2, stroke=1.5, color=circColor)

ggsave("~/Desktop/combined_de_cat_novel.pdf", width = 9, height = 6, useDingbats = FALSE)


# Save Combined Results DF ##############################################
saveRDS(combined_shrunk_results_df, combined_df_rds)

# Plot miRNAs of Specific Targets ########################
# load mirna target database
targets <- read_xlsx(mirna_targets_xlsx)

circColor <- "black"

#pdf(file.path(pdf_dir,"layer_specific_mirna_targets.pdf"), width=11, height=6)

presTheme <- presTheme + theme(legend.position = "bottom")

# CTIP2 (BCL11B)
tmp <- filter(targets, `Target Gene` == "BCL11B")

ctip2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting CTIP2",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% ctip2_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"ctip2_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# TBR1
tmp <- filter(targets, `Target Gene` == "TBR1")

tbr1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting TBR1",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% tbr1_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"tbr1_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# FOXP2
tmp <- filter(targets, `Target Gene` == "FOXP2")

foxp2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting FOXP2",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% foxp2_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"foxp2_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# CUX1
tmp <- filter(targets, `Target Gene` == "CUX1")

cux1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting CUX1",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% cux1_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"cux1_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# CUX2
tmp <- filter(targets, `Target Gene` == "CUX2")

cux2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting CUX2",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% cux2_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"cux2_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# SATB2
tmp <- filter(targets, `Target Gene` == "SATB2")

satb2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting SATB2",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% satb2_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"satb2_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# BRN1
tmp <- filter(targets, `Target Gene` == "POU3F3")

brn1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting BRN1",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% brn1_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"brn1_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# BRN2
tmp <- filter(targets, `Target Gene` == "POU3F2")

brn2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting BRN2",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% brn2_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"brn2_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# RELN
tmp <- filter(targets, `Target Gene` == "RELN")

reln_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=2) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="miRNAs Targeting RELN",
       color="Significant") +
  presTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)]) +
  geom_point(data=subset(combined_shrunk_results_df, mirna %in% reln_mirnas),
             aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

if(write_plots){ggsave(file.path(pdf_dir,"reln_targets_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

#dev.off()
