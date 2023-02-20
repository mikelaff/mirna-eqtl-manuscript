# diff expression gene expression gz/cp and gestation week

library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(AnnotationHub)

source(here("src/utils/lafferty_utils.R"))

# OUTPUT FILES ####################################################################################
# output directory for pdf files
pdf_dir <- here("doc/diff_expression/pdfs/")
# write plots
write_plots <- FALSE
# intermediate rdata files
dds_gw_rds <- here("doc/diff_expression/rdata/dds_gene_expression_gw.rds")
dds_gzcp_rds <- here("doc/diff_expression/rdata/dds_gene_expression_gzcp.rds")
# combined deseq2 shrunken results as a dataframe with labels for significance and category (neuro vs maturation)
combined_df_rds <- here("doc/diff_expression/rdata/combined_gene_expression_deseq2_shrunken_results_df.rds")

# INPUT FILES #####################################################################################
# ranged summarized experiment with gene expression counts
expression_rse_rds <- here("results/rdata_files/20181018_fetalTissue_ranged_summarized_experiment_gene_counts.rds")
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")

# miRNA data #######################################################################################
# read in summarized experiment
se <- readRDS(mirbase_se_rds)

# Build DDS gene expression data #######################################################################################
# read in ranged summarized experiment
rse <- readRDS(expression_rse_rds)

# build DESeq data set
dds <- DESeqDataSet(rse, design = ~1)
# remove rows with low counts
dds <- dds[rowSums(counts(dds)) > 10, ]
# filter for samples within mirna dataset
dds <- dds[,which(colnames(dds) %in% se$rnaid)]
# estimate size factors for normalization
#dds <- estimateSizeFactors(dds)
# norm. transformed counts (DESeq size factors, log2 transformed)
#ntc <- assay(normTransform(dds))
# variance stabilizing tranformed data
#vsd <- vst(dds)
rm(rse, se)

# Subset for GW or GZ/CP ########################################################################
# remove GZ/CP samples
dds.gw <- dds[, dds$tissue_section == "CW"]
dds.gw$tissue_section <- factor(dds.gw$tissue_section)
# remove CW samples
dds.gzcp <- dds[, dds$tissue_section == "GZ" | dds$tissue_section == "CP"]
dds.gzcp$tissue_section <- factor(dds.gzcp$tissue_section)

# Design Equations ################################################################################

# variable of interest: gestation_week
# round
# rin
# sex

dds.gw$round <- factor(dds.gw$round)
dds.gw$sex <- factor(dds.gw$sex)

# relabel round factor labels
levels(dds.gw$round)[levels(dds.gw$round) == "round2-redo"] <- "round2_redo"

design(dds.gw) <- formula(~ round + rin + sex + gestation_week)

# variable of interest: tissue_section
# round (cannot use because all samples from round1 come from the same donor, which is not present in the other rounds)
# rin
# donor_id

dds.gzcp$round <- factor(dds.gzcp$round)
dds.gzcp$sex <- factor(dds.gzcp$sex)
dds.gzcp$donor_id <- factor(dds.gzcp$donor_id)

# relabel round factor labels
levels(dds.gzcp$round)[levels(dds.gzcp$round) == "round2-redo"] <- "round2_redo"

design(dds.gzcp) <- formula(~ rin + donor_id + tissue_section)

# Run DESeq2 ######################################################################################
#dds.gw <- DESeq(dds.gw)
#dds.gzcp <- DESeq(dds.gzcp)

#saveRDS(dds.gw, dds_gw_rds)
#saveRDS(dds.gzcp, dds_gzcp_rds)

# load dds
dds.gw <- readRDS(dds_gw_rds)
dds.gzcp <- readRDS(dds_gzcp_rds)

# Results #########################################################################################
res.gw <- results(dds.gw)
res.gzcp <- results(dds.gzcp)
# shrunken results
shrunkres.gw <- lfcShrink(dds.gw, coef="gestation_week", type="apeglm")
shrunkres.gzcp <- lfcShrink(dds.gzcp, coef="tissue_section_GZ_vs_CP", type="apeglm")

#plotMA(res.gw)
#summary(res.gw)

#plotMA(res.gzcp)
#summary(res.gzcp)

plotMA(shrunkres.gw)
summary(shrunkres.gw)

plotMA(shrunkres.gzcp)
summary(shrunkres.gzcp)

# data frame from shrunken results
shrunkres.gw.df <- as.data.frame(res.gw)
shrunkres.gw.df$ensg <- rownames(shrunkres.gw.df)
shrunkres.gw.df$sig <- FALSE
shrunkres.gw.df$sig[which(shrunkres.gw.df$padj < 0.05)] <- TRUE
# order data frame for plotting
shrunkres.gw.df <- shrunkres.gw.df[order(shrunkres.gw.df$padj, decreasing = TRUE),]

# data frame from shrunken results
shrunkres.gzcp.df <- as.data.frame(res.gzcp)
shrunkres.gzcp.df$ensg <- rownames(shrunkres.gzcp.df)
shrunkres.gzcp.df$sig <- FALSE
shrunkres.gzcp.df$sig[which(shrunkres.gzcp.df$padj < 0.05)] <- TRUE
# order data frame for plotting
shrunkres.gzcp.df <- shrunkres.gzcp.df[order(shrunkres.gzcp.df$padj, decreasing = TRUE),]

rm(dds, dds.gw, dds.gzcp, res.gw, res.gzcp, shrunkres.gw, shrunkres.gzcp)

shrunkres.gzcp.df %>%
  filter(sig == TRUE, baseMean > 500) %>%
  arrange(log2FoldChange) -> tmp
write_lines(tmp$ensg[1:200], "~/Desktop/genes_cp_over_gz.txt")
  
shrunkres.gzcp.df %>%
  filter(sig == TRUE, baseMean > 500) %>%
  arrange(-log2FoldChange) -> tmp
write_lines(tmp$ensg[1:200], "~/Desktop/genes_gz_over_cp.txt")

# MA Plot #########################################################################################
circColor <- "black"
labelColor <- "black"

shrunkres.gw.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point(size=.5) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across Gestation Week") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"diff_gene_expression_gw_shrunklfc.pdf"), width = 11, height = 6, useDingbats = FALSE)}

shrunkres.gzcp.df %>%
  ggplot(aes(x=baseMean, y=log2FoldChange, col=sig)) +
  geom_point(size=.5) +
  labs(x="Mean of Normalized Counts",
       y="Shrunken Log2 Fold Change",
       title="Diff. Expression Across GZ/CP") +
  scale_x_log10() +
  geom_hline(yintercept = 0, size=.5, color="black") +
  theme(legend.position = "none") +
  plotTheme +
  scale_color_manual(values=c("grey70", "red"))

if(write_plots){ggsave(file.path(pdf_dir,"diff_gene_expression_gzcp_shrunklfc.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Compare GW to GZ/CP Differences #####################################################################

combined_shrunk_results_df <- left_join(shrunkres.gw.df, shrunkres.gzcp.df, suffix = c(".gw", ".gzcp"), by = "ensg")

# filter out rows not in both datasets (result of removing low expressors before diff. expression analysis)
combined_shrunk_results_df %<>% dplyr::filter(!is.na(log2FoldChange.gw) & !is.na(log2FoldChange.gzcp))

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gw)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.gzcp)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=c("grey60", "darkorange"))

# create column for significant in both analyses, none should be "pend" after proper labeling
combined_shrunk_results_df$sig.both <- "pend"

# label significant in both analyses column
combined_shrunk_results_df$sig.both[which(combined_shrunk_results_df$sig.gw & combined_shrunk_results_df$sig.gzcp)] <- "both"
combined_shrunk_results_df$sig.both[which(combined_shrunk_results_df$sig.gw & !combined_shrunk_results_df$sig.gzcp)] <- "gw only"
combined_shrunk_results_df$sig.both[which(!combined_shrunk_results_df$sig.gw & combined_shrunk_results_df$sig.gzcp)] <- "gzcp only"
combined_shrunk_results_df$sig.both[which(!combined_shrunk_results_df$sig.gw & !combined_shrunk_results_df$sig.gzcp)] <- "neither"
# factor sig.both variable
combined_shrunk_results_df$sig.both <- factor(combined_shrunk_results_df$sig.both)

# order rows for plotting
combined_shrunk_results_df <- combined_shrunk_results_df[order(combined_shrunk_results_df$sig.both, decreasing = TRUE),]

combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=sig.both)) +
  geom_point(size=.5) +
  labs(y="Log2 Fold Change GZ/CP",
       x="Log2 Fold Change Gestation Week",
       title="Differential Expression",
       color="Significant") +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,3)])

if(write_plots){ggsave(file.path(pdf_dir,"diff_gene_expression_gw_v_gzcp.pdf"), width = 11, height = 6, useDingbats = FALSE)}

# Subset Genes into Groups #######################################################################
# genes that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
genes_maturation_late <- combined_shrunk_results_df$ensg[which(combined_shrunk_results_df$sig.both == "gw only" & combined_shrunk_results_df$log2FoldChange.gw > 0)]
# genes that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
genes_maturation_early <- combined_shrunk_results_df$ensg[which(combined_shrunk_results_df$sig.both == "gw only" & combined_shrunk_results_df$log2FoldChange.gw < 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
genes_neurogenesis_gz <- combined_shrunk_results_df$ensg[which((combined_shrunk_results_df$sig.both == "both" | combined_shrunk_results_df$sig.both == "gzcp only") & combined_shrunk_results_df$log2FoldChange.gzcp > 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
genes_neurogenesis_cp <- combined_shrunk_results_df$ensg[which((combined_shrunk_results_df$sig.both == "both" | combined_shrunk_results_df$sig.both == "gzcp only") & combined_shrunk_results_df$log2FoldChange.gzcp < 0)]

# create category variable to each of these gene groups
combined_shrunk_results_df$category <- "none"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$ensg %in% genes_maturation_late)] <- "maturation_late"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$ensg %in% genes_maturation_early)] <- "maturation_early"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$ensg %in% genes_neurogenesis_gz)] <- "neurogenesis_gz"
combined_shrunk_results_df$category[which(combined_shrunk_results_df$ensg %in% genes_neurogenesis_cp)] <- "neurogenesis_cp"
# factor category variable
combined_shrunk_results_df$category <- factor(combined_shrunk_results_df$category)

# plot the categories
combined_shrunk_results_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=.5) +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

# Save Combined Results DF ##############################################
saveRDS(combined_shrunk_results_df, combined_df_rds)


# Gene Ontology ##################
df <- readRDS(combined_df_rds)

ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]

df %>%
  filter(category == "maturation_early") %>%
  arrange(log2FoldChange.gw) -> tmp

rowdata <- AnnotationDbi::select(orgdb, 
                                           keys = tmp$ensg, 
                                           columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                           keytype = "ENSEMBL")

rowdata <- rowdata[!is.na(rowdata$ENTREZID),]

write_lines(rowdata$ENTREZID[1:500], "~/Desktop/genes_mat_early.txt")

df %>%
  filter(category == "maturation_late") %>%
  arrange(log2FoldChange.gw) -> tmp

rowdata <- AnnotationDbi::select(orgdb, 
                                 keys = tmp$ensg, 
                                 columns = c("ENSEMBL", "SYMBOL", "ENTREZID"), 
                                 keytype = "ENSEMBL")

rowdata <- rowdata[!is.na(rowdata$ENTREZID),]

write_lines(rowdata$ENTREZID[1:500], "~/Desktop/genes_mat_late.txt")

