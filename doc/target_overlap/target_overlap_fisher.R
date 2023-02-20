# look at miRNA targets and their overlap with differentially expressed genes, and genes defined
# as being enriched in a cell type by scRNA-seq (from Luis at UCLA)

library(here)
library(dplyr)
library(magrittr)
library(readr)
library(ggplot2)
library(AnnotationHub)
library(miRNAtap)
library(EnsDb.Hsapiens.v86)
library(reshape2)

source(here("src/utils/lafferty_utils.R"))
prefix <- paste(format(Sys.time(), "%Y%m%d"), "target_overlap", sep="_")

# OUTPUT FILES ########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/target_overlap/pdfs/")
# write plots
write.plots <- FALSE
# data frame with fisher test results
df.results.rds <- paste(paste(here("doc/target_overlap/rdata/"), prefix, "_fisher_df.rds", sep=""))

# INPUT FILES #########################################################################################################
# total RNA-seq (gene expression) DE data
df.gene.diff.expression.rds <- here("doc/diff_expression/rdata/20190506_diff_gene_expression_df.rds")
# small RNA-seq (miRNA expression) DE data
df.mirna.diff.expression.rds <- here("doc/diff_expression/rdata/20190506_diff_expression_df.rds")
# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.csv <- here("data/ucla_single_cell/TableS6 Cluster enriched genes.csv")

# Import Data #########################################################################################################
df.gene <- readRDS(df.gene.diff.expression.rds)
df.mirna <- readRDS(df.mirna.diff.expression.rds)

df.gene$gene_biotype <- factor(df.gene$gene_biotype)

df.genes.by.cluster <- read_csv(genes.by.cluster.csv)
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Filter Genes ########################################################################################################
# edb <- EnsDb.Hsapiens.v86
# prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

df.gene %<>%
  dplyr::filter(gene_biotype == "protein_coding")

df.mirna %<>%
  dplyr::filter(sig.gzcp | sig.gw)

edb <- EnsDb.Hsapiens.v86
prot_genes_ensg <- names(genes(edb, filter = ~ gene_biotype == "protein_coding"))

df.genes.by.cluster %<>%
  dplyr::filter(Ensembl %in% prot_genes_ensg)

# Label miRNAs and Genes by Category ##################################################################################
# label significant in both analyses column
df.mirna$sig.both[which(df.mirna$sig.gw & df.mirna$sig.gzcp)] <- "both"
df.mirna$sig.both[which(df.mirna$sig.gw & !df.mirna$sig.gzcp)] <- "gw only"
df.mirna$sig.both[which(!df.mirna$sig.gw & df.mirna$sig.gzcp)] <- "gzcp only"
df.mirna$sig.both[which(!df.mirna$sig.gw & !df.mirna$sig.gzcp)] <- "neither"
df.mirna$sig.both <- factor(df.mirna$sig.both)

# mirnas that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
miRNAs.maturation.late <- df.mirna$Name[which(df.mirna$sig.both == "gw only" & df.mirna$log2FoldChange.gw > 0)]
# mirnas that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
miRNAs.maturation.early <- df.mirna$Name[which(df.mirna$sig.both == "gw only" & df.mirna$log2FoldChange.gw < 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
miRNAs.neurogenesis.gz <- df.mirna$Name[which((df.mirna$sig.both == "both" | df.mirna$sig.both == "gzcp only") & df.mirna$log2FoldChange.gzcp > 0)]
# mirnas that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
miRNAs.neurogenesis.cp <- df.mirna$Name[which((df.mirna$sig.both == "both" | df.mirna$sig.both == "gzcp only") & df.mirna$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df.mirna$category <- "none"
df.mirna$category[which(df.mirna$Name %in% miRNAs.maturation.late)] <- "maturation_late"
df.mirna$category[which(df.mirna$Name %in% miRNAs.maturation.early)] <- "maturation_early"
df.mirna$category[which(df.mirna$Name %in% miRNAs.neurogenesis.gz)] <- "neurogenesis_gz"
df.mirna$category[which(df.mirna$Name %in% miRNAs.neurogenesis.cp)] <- "neurogenesis_cp"
# factor category variable
df.mirna$category <- factor(df.mirna$category)

# label significant in both analyses column
df.gene$sig.both <- NA
df.gene$sig.both[which(df.gene$sig.gw & df.gene$sig.gzcp)] <- "both"
df.gene$sig.both[which(df.gene$sig.gw & !df.gene$sig.gzcp)] <- "gw only"
df.gene$sig.both[which(!df.gene$sig.gw & df.gene$sig.gzcp)] <- "gzcp only"
df.gene$sig.both[which(!df.gene$sig.gw & !df.gene$sig.gzcp)] <- "neither"
df.gene$sig.both <- factor(df.gene$sig.both)

# genes that are significant in only the gest. week analysis and have a positive log2 fold change (possibly genes associated with maturation in late neurogenesis)
genes.maturation.late <- df.gene$gene_id[which(df.gene$sig.both == "gw only" & df.gene$log2FoldChange.gw > 0)]
# genes that are significant in only the gest. week analysis and have a negative log2 fold change (possibly genes associated with maturation in early neurogenesis)
genes.maturation.early <- df.gene$gene_id[which(df.gene$sig.both == "gw only" & df.gene$log2FoldChange.gw < 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a positive log2 fold change (possibly genes associated with neurogensis and upregulated in gz tissue)
genes.neurogenesis.gz <- df.gene$gene_id[which((df.gene$sig.both == "both" | df.gene$sig.both == "gzcp only") & df.gene$log2FoldChange.gzcp > 0)]
# genes that are significant in only gz/cp analysis or significant in both analyses and have a negative log2 fold change (possibly genes associated with neurogensis and upregulated in cp tissue)
genes.neurogenesis.cp <- df.gene$gene_id[which((df.gene$sig.both == "both" | df.gene$sig.both == "gzcp only") & df.gene$log2FoldChange.gzcp < 0)]

# create category variable to each of these mirna groups
df.gene$category <- "none"
df.gene$category[which(df.gene$gene_id %in% genes.maturation.late)] <- "maturation_late"
df.gene$category[which(df.gene$gene_id %in% genes.maturation.early)] <- "maturation_early"
df.gene$category[which(df.gene$gene_id %in% genes.neurogenesis.gz)] <- "neurogenesis_gz"
df.gene$category[which(df.gene$gene_id %in% genes.neurogenesis.cp)] <- "neurogenesis_cp"
# factor category variable
df.gene$category <- factor(df.gene$category)

# Find Target Overlaps ######################################
# minimum number of sources required for a target to be considered
min_sources <- 2

# short name used for DB lookup
df.mirna$mir <- NA
# number of targets
df.mirna$num_targets <- NA
# number of targets overlapping gz>cp genes
df.mirna$num_targets_gz_over_cp <- NA
# number of targets overlapping cp>gz genes
df.mirna$num_targets_cp_over_gz <- NA
# number of targets overlapping maturation early genes
df.mirna$num_targets_mat_early <- NA
# number of targets overlapping maturation late genes
df.mirna$num_targets_mat_late <- NA
# number of targets not overlapping gz>cp genes
df.mirna$num_targets_not_gz_over_cp <- NA
# number of targets not overlapping cp>gz genes
df.mirna$num_targets_not_cp_over_gz <- NA
# number of targets not overlapping maturation early genes
df.mirna$num_targets_not_mat_early <- NA
# number of targets not overlapping maturation late genes
df.mirna$num_targets_not_mat_late <- NA

# number of targets overlapping cluster genes: ExN
df.mirna$num_targets_ExN <- NA
df.mirna$num_targets_not_ExN <- NA
# number of targets overlapping cluster genes: End
df.mirna$num_targets_End <- NA
df.mirna$num_targets_not_End <- NA
# number of targets overlapping cluster genes: ExCal
df.mirna$num_targets_ExCal <- NA
df.mirna$num_targets_not_ExCal <- NA
# number of targets overlapping cluster genes: ExDp1
df.mirna$num_targets_ExDp1 <- NA
df.mirna$num_targets_not_ExDp1 <- NA
# number of targets overlapping cluster genes: ExDp2
df.mirna$num_targets_ExDp2 <- NA
df.mirna$num_targets_not_ExDp2 <- NA
# number of targets overlapping cluster genes: ExM
df.mirna$num_targets_ExM <- NA
df.mirna$num_targets_not_ExM <- NA
# number of targets overlapping cluster genes: InCALB2
df.mirna$num_targets_InCALB2 <- NA
df.mirna$num_targets_not_InCALB2 <- NA
# number of targets overlapping cluster genes: InSST
df.mirna$num_targets_InSST <- NA
df.mirna$num_targets_not_InSST <- NA
# number of targets overlapping cluster genes: IP
df.mirna$num_targets_IP <- NA
df.mirna$num_targets_not_IP <- NA
# number of targets overlapping cluster genes: Mic
df.mirna$num_targets_Mic <- NA
df.mirna$num_targets_not_Mic <- NA
# number of targets overlapping cluster genes: OPC
df.mirna$num_targets_OPC <- NA
df.mirna$num_targets_not_OPC <- NA
# number of targets overlapping cluster genes: oRG
df.mirna$num_targets_oRG <- NA
df.mirna$num_targets_not_oRG <- NA
# number of targets overlapping cluster genes: Per
df.mirna$num_targets_Per <- NA
df.mirna$num_targets_not_Per <- NA
# number of targets overlapping cluster genes: PgG2M
df.mirna$num_targets_PgG2M <- NA
df.mirna$num_targets_not_PgG2M <- NA
# number of targets overlapping cluster genes: PgS
df.mirna$num_targets_PgS <- NA
df.mirna$num_targets_not_PgS <- NA
# number of targets overlapping cluster genes: vRG
df.mirna$num_targets_vRG <- NA
df.mirna$num_targets_not_vRG <- NA

# results of fisher's exact for gz>cp enrichment
df.mirna$FT_odds_gz_cp <- NA
df.mirna$FT_pval_gz_cp <- NA
df.mirna$FT_ci_low_gz_cp <- NA
df.mirna$FT_ci_high_gz_cp <- NA
# results of fisher's exact for cp>gz enrichment
df.mirna$FT_odds_cp_gz <- NA
df.mirna$FT_pval_cp_gz <- NA
df.mirna$FT_ci_low_cp_gz <- NA
df.mirna$FT_ci_high_cp_gz <- NA
# results of fisher's exact for mat_early
df.mirna$FT_odds_mat_early <- NA
df.mirna$FT_pval_mat_early <- NA
df.mirna$FT_ci_low_mat_early <- NA
df.mirna$FT_ci_high_mat_early <- NA
# results of fisher's exact for mat_late
df.mirna$FT_odds_mat_late <- NA
df.mirna$FT_pval_mat_late <- NA
df.mirna$FT_ci_low_mat_late <- NA
df.mirna$FT_ci_high_mat_late <- NA

# results of fisher's exact for cluster: ExN
df.mirna$FT_odds_ExN <- NA
df.mirna$FT_pval_ExN <- NA
df.mirna$FT_ci_low_ExN <- NA
df.mirna$FT_ci_high_ExN <- NA
# results of fisher's exact for cluster: End
df.mirna$FT_odds_End <- NA
df.mirna$FT_pval_End <- NA
df.mirna$FT_ci_low_End <- NA
df.mirna$FT_ci_high_End <- NA
# results of fisher's exact for cluster: ExCal
df.mirna$FT_odds_ExCal <- NA
df.mirna$FT_pval_ExCal <- NA
df.mirna$FT_ci_low_ExCal <- NA
df.mirna$FT_ci_high_ExCal <- NA
# results of fisher's exact for cluster: ExDp1
df.mirna$FT_odds_ExDp1 <- NA
df.mirna$FT_pval_ExDp1 <- NA
df.mirna$FT_ci_low_ExDp1 <- NA
df.mirna$FT_ci_high_ExDp1 <- NA
# results of fisher's exact for cluster: ExDp2
df.mirna$FT_odds_ExDp2 <- NA
df.mirna$FT_pval_ExDp2 <- NA
df.mirna$FT_ci_low_ExDp2 <- NA
df.mirna$FT_ci_high_ExDp2 <- NA
# results of fisher's exact for cluster: ExM
df.mirna$FT_odds_ExM <- NA
df.mirna$FT_pval_ExM <- NA
df.mirna$FT_ci_low_ExM <- NA
df.mirna$FT_ci_high_ExM <- NA
# results of fisher's exact for cluster: InCALB2
df.mirna$FT_odds_InCALB2 <- NA
df.mirna$FT_pval_InCALB2 <- NA
df.mirna$FT_ci_low_InCALB2 <- NA
df.mirna$FT_ci_high_InCALB2 <- NA
# results of fisher's exact for cluster: InSST
df.mirna$FT_odds_InSST <- NA
df.mirna$FT_pval_InSST <- NA
df.mirna$FT_ci_low_InSST <- NA
df.mirna$FT_ci_high_InSST <- NA
# results of fisher's exact for cluster: IP
df.mirna$FT_odds_IP <- NA
df.mirna$FT_pval_IP <- NA
df.mirna$FT_ci_low_IP <- NA
df.mirna$FT_ci_high_IP <- NA
# results of fisher's exact for cluster: Mic
df.mirna$FT_odds_Mic <- NA
df.mirna$FT_pval_Mic <- NA
df.mirna$FT_ci_low_Mic <- NA
df.mirna$FT_ci_high_Mic <- NA
# results of fisher's exact for cluster: OPC
df.mirna$FT_odds_OPC <- NA
df.mirna$FT_pval_OPC <- NA
df.mirna$FT_ci_low_OPC <- NA
df.mirna$FT_ci_high_OPC <- NA
# results of fisher's exact for cluster: oRG
df.mirna$FT_odds_oRG <- NA
df.mirna$FT_pval_oRG <- NA
df.mirna$FT_ci_low_oRG <- NA
df.mirna$FT_ci_high_oRG <- NA
# results of fisher's exact for cluster: Per
df.mirna$FT_odds_Per <- NA
df.mirna$FT_pval_Per <- NA
df.mirna$FT_ci_low_Per <- NA
df.mirna$FT_ci_high_Per <- NA
# results of fisher's exact for cluster: PgG2M
df.mirna$FT_odds_PgG2M <- NA
df.mirna$FT_pval_PgG2M <- NA
df.mirna$FT_ci_low_PgG2M <- NA
df.mirna$FT_ci_high_PgG2M <- NA
# results of fisher's exact for cluster: PgS
df.mirna$FT_odds_PgS <- NA
df.mirna$FT_pval_PgS <- NA
df.mirna$FT_ci_low_PgS <- NA
df.mirna$FT_ci_high_PgS <- NA
# results of fisher's exact for cluster: vRG
df.mirna$FT_odds_vRG <- NA
df.mirna$FT_pval_vRG <- NA
df.mirna$FT_ci_low_vRG <- NA
df.mirna$FT_ci_high_vRG <- NA

# all genes in each group tested
genes_gz_greater_than_cp <- df.gene$gene_id[which(df.gene$category == "neurogenesis_gz")]
genes_cp_greater_than_gz <- df.gene$gene_id[which(df.gene$category == "neurogenesis_cp")]
genes_maturation_early <- df.gene$gene_id[which(df.gene$category == "maturation_early")]
genes_maturation_late <- df.gene$gene_id[which(df.gene$category == "maturation_late")]

genes_not_gz_greater_than_cp <- df.gene$gene_id[which(df.gene$category != "neurogenesis_gz")]
genes_not_cp_greater_than_gz <- df.gene$gene_id[which(df.gene$category != "neurogenesis_cp")]
genes_not_maturation_early <- df.gene$gene_id[which(df.gene$category != "maturation_early")]
genes_not_maturation_late <- df.gene$gene_id[which(df.gene$category != "maturation_late")]

genes_ExN <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "ExN")]
genes_End <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "End")]
genes_ExCal <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "ExCal")]
genes_ExDp1 <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "ExDp1")]
genes_ExDp2 <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "ExDp2")]
genes_ExM <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "ExM")]
genes_InCALB2 <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "InCALB2")]
genes_InSST <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "InSST")]
genes_IP <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "IP")]
genes_Mic <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "Mic")]
genes_OPC <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "OPC")]
genes_oRG <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "oRG")]
genes_Per <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "Per")]
genes_PgG2M <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "PgG2M")]
genes_PgS <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "PgS")]
genes_vRG <- df.genes.by.cluster$Ensembl[which(df.genes.by.cluster$Cluster == "vRG")]

genes_not_ExN <- prot_genes_ensg[!prot_genes_ensg %in% genes_ExN]
genes_not_End <- prot_genes_ensg[!prot_genes_ensg %in% genes_End]
genes_not_ExCal <- prot_genes_ensg[!prot_genes_ensg %in% genes_ExCal]
genes_not_ExDp1 <- prot_genes_ensg[!prot_genes_ensg %in% genes_ExDp1]
genes_not_ExDp2 <- prot_genes_ensg[!prot_genes_ensg %in% genes_ExDp2]
genes_not_ExM <- prot_genes_ensg[!prot_genes_ensg %in% genes_ExM]
genes_not_InCALB2 <- prot_genes_ensg[!prot_genes_ensg %in% genes_InCALB2]
genes_not_InSST <- prot_genes_ensg[!prot_genes_ensg %in% genes_InSST]
genes_not_IP <- prot_genes_ensg[!prot_genes_ensg %in% genes_IP]
genes_not_Mic <- prot_genes_ensg[!prot_genes_ensg %in% genes_Mic]
genes_not_OPC <- prot_genes_ensg[!prot_genes_ensg %in% genes_OPC]
genes_not_oRG <- prot_genes_ensg[!prot_genes_ensg %in% genes_oRG]
genes_not_Per <- prot_genes_ensg[!prot_genes_ensg %in% genes_Per]
genes_not_PgG2M <- prot_genes_ensg[!prot_genes_ensg %in% genes_PgG2M]
genes_not_PgS <- prot_genes_ensg[!prot_genes_ensg %in% genes_PgS]
genes_not_vRG <- prot_genes_ensg[!prot_genes_ensg %in% genes_vRG]

# totals
num_genes_gz_greater_than_cp <- length(genes_gz_greater_than_cp)
num_genes_cp_greater_than_gz <- length(genes_cp_greater_than_gz)
num_genes_maturation_early <- length(genes_maturation_early)
num_genes_maturation_late <- length(genes_maturation_late)

num_genes_not_gz_greater_than_cp <- length(genes_not_gz_greater_than_cp)
num_genes_not_cp_greater_than_gz <- length(genes_not_cp_greater_than_gz)
num_genes_not_maturation_early <- length(genes_not_maturation_early)
num_genes_not_maturation_late <- length(genes_not_maturation_late)

num_genes_ExN <- length(genes_ExN)
num_genes_End <- length(genes_End)
num_genes_ExCal <- length(genes_ExCal)
num_genes_ExDp1 <- length(genes_ExDp1)
num_genes_ExDp2 <- length(genes_ExDp2)
num_genes_ExM <- length(genes_ExM)
num_genes_InCALB2 <- length(genes_InCALB2)
num_genes_InSST <- length(genes_InSST)
num_genes_IP <- length(genes_IP)
num_genes_Mic <- length(genes_Mic)
num_genes_OPC <- length(genes_OPC)
num_genes_oRG <- length(genes_oRG)
num_genes_Per <- length(genes_Per)
num_genes_PgG2M <- length(genes_PgG2M)
num_genes_PgS <- length(genes_PgS)
num_genes_vRG <- length(genes_vRG)

num_genes_not_ExN <- length(genes_not_ExN)
num_genes_not_End <- length(genes_not_End)
num_genes_not_ExCal <- length(genes_not_ExCal)
num_genes_not_ExDp1 <- length(genes_not_ExDp1)
num_genes_not_ExDp2 <- length(genes_not_ExDp2)
num_genes_not_ExM <- length(genes_not_ExM)
num_genes_not_InCALB2 <- length(genes_not_InCALB2)
num_genes_not_InSST <- length(genes_not_InSST)
num_genes_not_IP <- length(genes_not_IP)
num_genes_not_Mic <- length(genes_not_Mic)
num_genes_not_OPC <- length(genes_not_OPC)
num_genes_not_oRG <- length(genes_not_oRG)
num_genes_not_Per <- length(genes_not_Per)
num_genes_not_PgG2M <- length(genes_not_PgG2M)
num_genes_not_PgS <- length(genes_not_PgS)
num_genes_not_vRG <- length(genes_not_vRG)

# looking at only sig. expressed mirnas (all categories)
# loop and find targets and overlaps
indexes <- 1:length(df.mirna$Name)
for (i in indexes) {

  # if novel, skip
  if (df.mirna$source[i] != "miRBase_v22") {
    next
  }
  ft_gzcp <- NULL
  ft_cpgz <- NULL
  ft_early <- NULL
  ft_late <- NULL

  ft_End <- NULL
  ft_ExCal <- NULL
  ft_ExDp1 <- NULL
  ft_ExDp2 <- NULL
  ft_ExM <- NULL
  ft_ExN <- NULL
  ft_InCALB2 <- NULL
  ft_InSST <- NULL
  ft_IP <- NULL
  ft_Mic <- NULL
  ft_OPC <- NULL
  ft_oRG <- NULL
  ft_Per <- NULL
  ft_PgG2M <- NULL
  ft_PgS <- NULL
  ft_vRG <- NULL

  predictions <- NULL
  rowdata <- NULL

  # short name used for database lookup
  df.mirna$mir[i] <- strsplit(strsplit(df.mirna$Name[i], '/')[[1]][1], "hsa-")[[1]][2]
  print(paste(i, df.mirna$mir[i], sep=": "))
  # get predicted targets
  predictions <- data.frame(getPredictedTargets(df.mirna$mir[i], species = "hsa", method = "geom", min_src = min_sources))
  if (length(predictions) == 0) {
    print(paste("No Targets:", df.mirna$mirna[i]))
    next
  }
  predictions$ENTREZID <- rownames(predictions)
  # get ensembl ids
  rowdata <- AnnotationDbi::select(orgdb,
                                   keys = predictions$ENTREZID,
                                   keytype = "ENTREZID",
                                   columns = c("ENSEMBL"))
  # join with predictions
  predictions <- left_join(predictions, rowdata, by = "ENTREZID")

  # populate table
  # number of targets
  df.mirna$num_targets[i] <- length(predictions$ENTREZID)
  # number of targets overlapping gz>cp genes
  df.mirna$num_targets_gz_over_cp[i] <- sum(predictions$ENSEMBL %in% genes_gz_greater_than_cp)
  # number of targets overlapping cp>gz genes
  df.mirna$num_targets_cp_over_gz[i] <- sum(predictions$ENSEMBL %in% genes_cp_greater_than_gz)
  # number of targets overlapping maturation early genes
  df.mirna$num_targets_mat_early[i] <- sum(predictions$ENSEMBL %in% genes_maturation_early)
  # number of targets overlapping maturation late genes
  df.mirna$num_targets_mat_late[i] <- sum(predictions$ENSEMBL %in% genes_maturation_late)
  # number of targets not overlapping gz>cp genes
  df.mirna$num_targets_not_gz_over_cp[i] <- sum(predictions$ENSEMBL %in% genes_not_gz_greater_than_cp)
  # number of targets not overlapping cp>gz genes
  df.mirna$num_targets_not_cp_over_gz[i] <- sum(predictions$ENSEMBL %in% genes_not_cp_greater_than_gz)
  # number of targets not overlapping maturation early genes
  df.mirna$num_targets_not_mat_early[i] <- sum(predictions$ENSEMBL %in% genes_not_maturation_early)
  # number of targets not overlapping maturation late genes
  df.mirna$num_targets_not_mat_late[i] <- sum(predictions$ENSEMBL %in% genes_not_maturation_late)

  df.mirna$num_targets_End[i] <- sum(predictions$ENSEMBL %in% genes_End)
  df.mirna$num_targets_not_End[i] <- sum(predictions$ENSEMBL %in% genes_not_End)
  df.mirna$num_targets_ExCal[i] <- sum(predictions$ENSEMBL %in% genes_ExCal)
  df.mirna$num_targets_not_ExCal[i] <- sum(predictions$ENSEMBL %in% genes_not_ExCal)
  df.mirna$num_targets_ExDp1[i] <- sum(predictions$ENSEMBL %in% genes_ExDp1)
  df.mirna$num_targets_not_ExDp1[i] <- sum(predictions$ENSEMBL %in% genes_not_ExDp1)
  df.mirna$num_targets_ExDp2[i] <- sum(predictions$ENSEMBL %in% genes_ExDp2)
  df.mirna$num_targets_not_ExDp2[i] <- sum(predictions$ENSEMBL %in% genes_not_ExDp2)
  df.mirna$num_targets_ExM[i] <- sum(predictions$ENSEMBL %in% genes_ExM)
  df.mirna$num_targets_not_ExM[i] <- sum(predictions$ENSEMBL %in% genes_not_ExM)
  df.mirna$num_targets_ExN[i] <- sum(predictions$ENSEMBL %in% genes_ExN)
  df.mirna$num_targets_not_ExN[i] <- sum(predictions$ENSEMBL %in% genes_not_ExN)
  df.mirna$num_targets_InCALB2[i] <- sum(predictions$ENSEMBL %in% genes_InCALB2)
  df.mirna$num_targets_not_InCALB2[i] <- sum(predictions$ENSEMBL %in% genes_not_InCALB2)
  df.mirna$num_targets_InSST[i] <- sum(predictions$ENSEMBL %in% genes_InSST)
  df.mirna$num_targets_not_InSST[i] <- sum(predictions$ENSEMBL %in% genes_not_InSST)
  df.mirna$num_targets_IP[i] <- sum(predictions$ENSEMBL %in% genes_IP)
  df.mirna$num_targets_not_IP[i] <- sum(predictions$ENSEMBL %in% genes_not_IP)
  df.mirna$num_targets_Mic[i] <- sum(predictions$ENSEMBL %in% genes_Mic)
  df.mirna$num_targets_not_Mic[i] <- sum(predictions$ENSEMBL %in% genes_not_Mic)
  df.mirna$num_targets_OPC[i] <- sum(predictions$ENSEMBL %in% genes_OPC)
  df.mirna$num_targets_not_OPC[i] <- sum(predictions$ENSEMBL %in% genes_not_OPC)
  df.mirna$num_targets_oRG[i] <- sum(predictions$ENSEMBL %in% genes_oRG)
  df.mirna$num_targets_not_oRG[i] <- sum(predictions$ENSEMBL %in% genes_not_oRG)
  df.mirna$num_targets_Per[i] <- sum(predictions$ENSEMBL %in% genes_Per)
  df.mirna$num_targets_not_Per[i] <- sum(predictions$ENSEMBL %in% genes_not_Per)
  df.mirna$num_targets_PgG2M[i] <- sum(predictions$ENSEMBL %in% genes_PgG2M)
  df.mirna$num_targets_not_PgG2M[i] <- sum(predictions$ENSEMBL %in% genes_not_PgG2M)
  df.mirna$num_targets_PgS[i] <- sum(predictions$ENSEMBL %in% genes_PgS)
  df.mirna$num_targets_not_PgS[i] <- sum(predictions$ENSEMBL %in% genes_not_PgS)
  df.mirna$num_targets_vRG[i] <- sum(predictions$ENSEMBL %in% genes_vRG)
  df.mirna$num_targets_not_vRG[i] <- sum(predictions$ENSEMBL %in% genes_not_vRG)

  # # cont. table: genes are either targets or not targets, and genes are either gz>cp or other
  # tab <- matrix(c(num_genes_not_gz_greater_than_cp-df.mirna$num_targets_not_gz_over_cp[i],
  #                 num_genes_gz_greater_than_cp-df.mirna$num_targets_gz_over_cp[i],
  #                 df.mirna$num_targets_not_gz_over_cp[i],
  #                 df.mirna$num_targets_gz_over_cp[i]), byrow = TRUE, nrow = 2)
  #
  # # labels rows and columns for this table
  # colnames(tab) <- c("other", "gz>cp")
  # rownames(tab) <- c("not_targets", "targets")
  # tab

  # gz>cp fisher's exact test
  ft_gzcp <- fisher.test(matrix(c(num_genes_not_gz_greater_than_cp-df.mirna$num_targets_not_gz_over_cp[i],
                                  num_genes_gz_greater_than_cp-df.mirna$num_targets_gz_over_cp[i],
                                  df.mirna$num_targets_not_gz_over_cp[i],
                                  df.mirna$num_targets_gz_over_cp[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_gz_cp[i] <- ft_gzcp$estimate
  df.mirna$FT_pval_gz_cp[i] <- ft_gzcp$p.value
  df.mirna$FT_ci_low_gz_cp[i] <- ft_gzcp$conf.int[1]
  df.mirna$FT_ci_high_gz_cp[i] <- ft_gzcp$conf.int[2]

  # cp>gz fisher's exact test
  ft_cpgz <- fisher.test(matrix(c(num_genes_not_cp_greater_than_gz-df.mirna$num_targets_not_cp_over_gz[i],
                                  num_genes_cp_greater_than_gz-df.mirna$num_targets_cp_over_gz[i],
                                  df.mirna$num_targets_not_cp_over_gz[i],
                                  df.mirna$num_targets_cp_over_gz[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_cp_gz[i] <- ft_cpgz$estimate
  df.mirna$FT_pval_cp_gz[i] <- ft_cpgz$p.value
  df.mirna$FT_ci_low_cp_gz[i] <- ft_cpgz$conf.int[1]
  df.mirna$FT_ci_high_cp_gz[i] <- ft_cpgz$conf.int[2]

  # maturation early fisher's exact test
  ft_early <- fisher.test(matrix(c(num_genes_not_maturation_early-df.mirna$num_targets_not_mat_early[i],
                                   num_genes_maturation_early-df.mirna$num_targets_mat_early[i],
                                   df.mirna$num_targets_not_mat_early[i],
                                   df.mirna$num_targets_mat_early[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_mat_early[i] <- ft_early$estimate
  df.mirna$FT_pval_mat_early[i] <- ft_early$p.value
  df.mirna$FT_ci_low_mat_early[i] <- ft_early$conf.int[1]
  df.mirna$FT_ci_high_mat_early[i] <- ft_early$conf.int[2]

  # maturation late fisher's exact test
  ft_late <- fisher.test(matrix(c(num_genes_not_maturation_late-df.mirna$num_targets_not_mat_late[i],
                                  num_genes_maturation_late-df.mirna$num_targets_mat_late[i],
                                  df.mirna$num_targets_not_mat_late[i],
                                  df.mirna$num_targets_mat_late[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_mat_late[i] <- ft_late$estimate
  df.mirna$FT_pval_mat_late[i] <- ft_late$p.value
  df.mirna$FT_ci_low_mat_late[i] <- ft_late$conf.int[1]
  df.mirna$FT_ci_high_mat_late[i] <- ft_late$conf.int[2]

  # End
  ft_End <- fisher.test(matrix(c(num_genes_not_End - df.mirna$num_targets_not_End[i],
                                  num_genes_End - df.mirna$num_targets_End[i],
                                  df.mirna$num_targets_not_End[i],
                                  df.mirna$num_targets_End[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_End[i] <- ft_End$estimate
  df.mirna$FT_pval_End[i] <- ft_End$p.value
  df.mirna$FT_ci_low_End[i] <- ft_End$conf.int[1]
  df.mirna$FT_ci_high_End[i] <- ft_End$conf.int[2]

  # ExCal
  ft_ExCal <- fisher.test(matrix(c(num_genes_not_ExCal - df.mirna$num_targets_not_ExCal[i],
                                 num_genes_ExCal - df.mirna$num_targets_ExCal[i],
                                 df.mirna$num_targets_not_ExCal[i],
                                 df.mirna$num_targets_ExCal[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_ExCal[i] <- ft_ExCal$estimate
  df.mirna$FT_pval_ExCal[i] <- ft_ExCal$p.value
  df.mirna$FT_ci_low_ExCal[i] <- ft_ExCal$conf.int[1]
  df.mirna$FT_ci_high_ExCal[i] <- ft_ExCal$conf.int[2]

  # ExDp1
  ft_ExDp1 <- fisher.test(matrix(c(num_genes_not_ExDp1 - df.mirna$num_targets_not_ExDp1[i],
                              num_genes_ExDp1 - df.mirna$num_targets_ExDp1[i],
                              df.mirna$num_targets_not_ExDp1[i],
                              df.mirna$num_targets_ExDp1[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_ExDp1[i] <- ft_ExDp1$estimate
  df.mirna$FT_pval_ExDp1[i] <- ft_ExDp1$p.value
  df.mirna$FT_ci_low_ExDp1[i] <- ft_ExDp1$conf.int[1]
  df.mirna$FT_ci_high_ExDp1[i] <- ft_ExDp1$conf.int[2]

  # ExDp2
  ft_ExDp2 <- fisher.test(matrix(c(num_genes_not_ExDp2 - df.mirna$num_targets_not_ExDp2[i],
                              num_genes_ExDp2 - df.mirna$num_targets_ExDp2[i],
                              df.mirna$num_targets_not_ExDp2[i],
                              df.mirna$num_targets_ExDp2[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_ExDp2[i] <- ft_ExDp2$estimate
  df.mirna$FT_pval_ExDp2[i] <- ft_ExDp2$p.value
  df.mirna$FT_ci_low_ExDp2[i] <- ft_ExDp2$conf.int[1]
  df.mirna$FT_ci_high_ExDp2[i] <- ft_ExDp2$conf.int[2]

  # ExM
  ft_ExM <- fisher.test(matrix(c(num_genes_not_ExM - df.mirna$num_targets_not_ExM[i],
                              num_genes_ExM - df.mirna$num_targets_ExM[i],
                              df.mirna$num_targets_not_ExM[i],
                              df.mirna$num_targets_ExM[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_ExM[i] <- ft_ExM$estimate
  df.mirna$FT_pval_ExM[i] <- ft_ExM$p.value
  df.mirna$FT_ci_low_ExM[i] <- ft_ExM$conf.int[1]
  df.mirna$FT_ci_high_ExM[i] <- ft_ExM$conf.int[2]

  # ExN
  ft_ExN <- fisher.test(matrix(c(num_genes_not_ExN - df.mirna$num_targets_not_ExN[i],
                              num_genes_ExN - df.mirna$num_targets_ExN[i],
                              df.mirna$num_targets_not_ExN[i],
                              df.mirna$num_targets_ExN[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_ExN[i] <- ft_ExN$estimate
  df.mirna$FT_pval_ExN[i] <- ft_ExN$p.value
  df.mirna$FT_ci_low_ExN[i] <- ft_ExN$conf.int[1]
  df.mirna$FT_ci_high_ExN[i] <- ft_ExN$conf.int[2]

  # InCALB2
  ft_InCALB2 <- fisher.test(matrix(c(num_genes_not_InCALB2 - df.mirna$num_targets_not_InCALB2[i],
                              num_genes_InCALB2 - df.mirna$num_targets_InCALB2[i],
                              df.mirna$num_targets_not_InCALB2[i],
                              df.mirna$num_targets_InCALB2[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_InCALB2[i] <- ft_InCALB2$estimate
  df.mirna$FT_pval_InCALB2[i] <- ft_InCALB2$p.value
  df.mirna$FT_ci_low_InCALB2[i] <- ft_InCALB2$conf.int[1]
  df.mirna$FT_ci_high_InCALB2[i] <- ft_InCALB2$conf.int[2]

  # InSST
  ft_InSST <- fisher.test(matrix(c(num_genes_not_InSST - df.mirna$num_targets_not_InSST[i],
                              num_genes_InSST - df.mirna$num_targets_InSST[i],
                              df.mirna$num_targets_not_InSST[i],
                              df.mirna$num_targets_InSST[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_InSST[i] <- ft_InSST$estimate
  df.mirna$FT_pval_InSST[i] <- ft_InSST$p.value
  df.mirna$FT_ci_low_InSST[i] <- ft_InSST$conf.int[1]
  df.mirna$FT_ci_high_InSST[i] <- ft_InSST$conf.int[2]

  # IP
  ft_IP <- fisher.test(matrix(c(num_genes_not_IP - df.mirna$num_targets_not_IP[i],
                              num_genes_IP - df.mirna$num_targets_IP[i],
                              df.mirna$num_targets_not_IP[i],
                              df.mirna$num_targets_IP[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_IP[i] <- ft_IP$estimate
  df.mirna$FT_pval_IP[i] <- ft_IP$p.value
  df.mirna$FT_ci_low_IP[i] <- ft_IP$conf.int[1]
  df.mirna$FT_ci_high_IP[i] <- ft_IP$conf.int[2]

  # Mic
  ft_Mic <- fisher.test(matrix(c(num_genes_not_Mic - df.mirna$num_targets_not_Mic[i],
                              num_genes_Mic - df.mirna$num_targets_Mic[i],
                              df.mirna$num_targets_not_Mic[i],
                              df.mirna$num_targets_Mic[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_Mic[i] <- ft_Mic$estimate
  df.mirna$FT_pval_Mic[i] <- ft_Mic$p.value
  df.mirna$FT_ci_low_Mic[i] <- ft_Mic$conf.int[1]
  df.mirna$FT_ci_high_Mic[i] <- ft_Mic$conf.int[2]

  # OPC
  ft_OPC <- fisher.test(matrix(c(num_genes_not_OPC - df.mirna$num_targets_not_OPC[i],
                              num_genes_OPC - df.mirna$num_targets_OPC[i],
                              df.mirna$num_targets_not_OPC[i],
                              df.mirna$num_targets_OPC[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_OPC[i] <- ft_OPC$estimate
  df.mirna$FT_pval_OPC[i] <- ft_OPC$p.value
  df.mirna$FT_ci_low_OPC[i] <- ft_OPC$conf.int[1]
  df.mirna$FT_ci_high_OPC[i] <- ft_OPC$conf.int[2]

  # oRG
  ft_oRG <- fisher.test(matrix(c(num_genes_not_oRG - df.mirna$num_targets_not_oRG[i],
                              num_genes_oRG - df.mirna$num_targets_oRG[i],
                              df.mirna$num_targets_not_oRG[i],
                              df.mirna$num_targets_oRG[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_oRG[i] <- ft_oRG$estimate
  df.mirna$FT_pval_oRG[i] <- ft_oRG$p.value
  df.mirna$FT_ci_low_oRG[i] <- ft_oRG$conf.int[1]
  df.mirna$FT_ci_high_oRG[i] <- ft_oRG$conf.int[2]

  # Per
  ft_Per <- fisher.test(matrix(c(num_genes_not_Per - df.mirna$num_targets_not_Per[i],
                              num_genes_Per - df.mirna$num_targets_Per[i],
                              df.mirna$num_targets_not_Per[i],
                              df.mirna$num_targets_Per[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_Per[i] <- ft_Per$estimate
  df.mirna$FT_pval_Per[i] <- ft_Per$p.value
  df.mirna$FT_ci_low_Per[i] <- ft_Per$conf.int[1]
  df.mirna$FT_ci_high_Per[i] <- ft_Per$conf.int[2]

  # PgG2M
  ft_PgG2M <- fisher.test(matrix(c(num_genes_not_PgG2M - df.mirna$num_targets_not_PgG2M[i],
                              num_genes_PgG2M - df.mirna$num_targets_PgG2M[i],
                              df.mirna$num_targets_not_PgG2M[i],
                              df.mirna$num_targets_PgG2M[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_PgG2M[i] <- ft_PgG2M$estimate
  df.mirna$FT_pval_PgG2M[i] <- ft_PgG2M$p.value
  df.mirna$FT_ci_low_PgG2M[i] <- ft_PgG2M$conf.int[1]
  df.mirna$FT_ci_high_PgG2M[i] <- ft_PgG2M$conf.int[2]

  # PgS
  ft_PgS <- fisher.test(matrix(c(num_genes_not_PgS - df.mirna$num_targets_not_PgS[i],
                              num_genes_PgS - df.mirna$num_targets_PgS[i],
                              df.mirna$num_targets_not_PgS[i],
                              df.mirna$num_targets_PgS[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_PgS[i] <- ft_PgS$estimate
  df.mirna$FT_pval_PgS[i] <- ft_PgS$p.value
  df.mirna$FT_ci_low_PgS[i] <- ft_PgS$conf.int[1]
  df.mirna$FT_ci_high_PgS[i] <- ft_PgS$conf.int[2]

  # vRG
  ft_vRG <- fisher.test(matrix(c(num_genes_not_vRG - df.mirna$num_targets_not_vRG[i],
                              num_genes_vRG - df.mirna$num_targets_vRG[i],
                              df.mirna$num_targets_not_vRG[i],
                              df.mirna$num_targets_vRG[i]), byrow = TRUE, nrow = 2))

  df.mirna$FT_odds_vRG[i] <- ft_vRG$estimate
  df.mirna$FT_pval_vRG[i] <- ft_vRG$p.value
  df.mirna$FT_ci_low_vRG[i] <- ft_vRG$conf.int[1]
  df.mirna$FT_ci_high_vRG[i] <- ft_vRG$conf.int[2]
}

# save data
saveRDS(df.mirna, df.results.rds)

stop()
# load data
df.mirna <- readRDS()



tmp <- basicProfile(predictions$ENTREZID, onto = "BP", level = 2, orgPackage = "org.Hs.eg.db")

stop()

# basic l2fc plot
mirna_df %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point() +
  plotTheme +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)])

pdf("~/Desktop/mirna_target_enrichment.pdf")

# num_gzcp_test <- length(dplyr::filter(mirna_df, category == "neurogenesis_gz", !is.na(FT_gz_cp_odds))$mirna)
# num_cpgz_test <- length(dplyr::filter(mirna_df, category == "neurogenesis_cp", !is.na(FT_cp_gz_odds))$mirna)
# num_early_test <- length(dplyr::filter(mirna_df, category == "maturation_early", !is.na(FT_mat_early_odds))$mirna)
# num_late_test <- length(dplyr::filter(mirna_df, category == "maturation_late", !is.na(FT_mat_late_odds))$mirna)

# mirna upregulated gz: GZ>CP
mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated cp: CP>GZ
mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated early: maturation_early
mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# mirna upregulated late: maturation_late
mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "BH")))) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

dev.off()

###########################
mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_gz_cp_odds,
                pval = FT_gz_cp_pval,
                ci_low = FT_gz_cp_ci_low,
                ci_high = FT_gz_cp_ci_high) %>%
  dplyr::mutate(target_cat = "gz>cp") -> tmp1

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_cp_gz_odds,
                pval = FT_cp_gz_pval,
                ci_low = FT_cp_gz_ci_low,
                ci_high = FT_cp_gz_ci_high) %>%
  dplyr::mutate(target_cat = "cp>gz") -> tmp2

tmp3 <- bind_rows(tmp1, tmp2)

tmp3 %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), color = target_cat)) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

#################
mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_gz_cp_odds,
                pval = FT_gz_cp_pval,
                ci_low = FT_gz_cp_ci_low,
                ci_high = FT_gz_cp_ci_high) %>%
  dplyr::mutate(target_cat = "gz>cp") -> tmp1

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(mirna = mir,
                odds = FT_cp_gz_odds,
                pval = FT_cp_gz_pval,
                ci_low = FT_cp_gz_ci_low,
                ci_high = FT_cp_gz_ci_high) %>%
  dplyr::mutate(target_cat = "cp>gz") -> tmp2

tmp3 <- bind_rows(tmp1, tmp2)

tmp3 %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), color = target_cat)) +
  geom_point() +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

# SfN Plotting ##########################

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="blue") +
  geom_text(aes(label=mir), hjust=-0.5)+
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","1.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_gz", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="red") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA GZ>CP: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","2.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_gz_cp_odds)) %>%
  dplyr::select(odds = FT_gz_cp_odds, pval = FT_gz_cp_pval, ci_low = FT_gz_cp_ci_low, ci_high = FT_gz_cp_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="green") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in GZ>CP Genes",
       caption="odds > 0: targets enriched among GZ>CP genes\nodds < 0: targets depleted among GZ>CP genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","3.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "neurogenesis_cp", !is.na(FT_cp_gz_odds)) %>%
  dplyr::select(odds = FT_cp_gz_odds, pval = FT_cp_gz_pval, ci_low = FT_cp_gz_ci_low, ci_high = FT_cp_gz_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="orange") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA CP>GZ: Target Enrichment in CP>GZ Genes",
       caption="odds > 0: targets enriched among CP>GZ genes\nodds < 0: targets depleted among CP>GZ genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1))

ggsave(file.path("~/Desktop/","4.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="purple") +
  geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","5.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_early", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="navy") +
    geom_text(aes(label=mir), hjust=-0.5) +
  plotTheme +
  labs(title="miRNA Up Early: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","6.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_early_odds)) %>%
  dplyr::select(odds = FT_mat_early_odds, pval = FT_mat_early_pval, ci_low = FT_mat_early_ci_low, ci_high = FT_mat_early_ci_high, mir=mir) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "fdr")), label=mir)) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="pink") +
  plotTheme +
  geom_text(aes(label=mir), hjust=-1, size=2) +
  labs(title="miRNA Up Late: Target Enrichment in Early Genes",
       caption="odds > 0: targets enriched among early genes\nodds < 0: targets depleted among early genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","7.pdf"), width = 5, height = 5, useDingbats = FALSE)

mirna_df %>%
  dplyr::filter(category == "maturation_late", !is.na(FT_mat_late_odds)) %>%
  dplyr::select(odds = FT_mat_late_odds, pval = FT_mat_late_pval, ci_low = FT_mat_late_ci_low, ci_high = FT_mat_late_ci_high) %>%
  ggplot(aes(x=log10(odds), y=-log10(p.adjust(pval, method = "BH")))) +
  geom_errorbarh(aes(xmin=log10(ci_low), xmax=log10(ci_high))) +
  geom_point(color="yellow") +
  plotTheme +
  labs(title="miRNA Up Late: Target Enrichment in Late Genes",
       caption="odds > 0: targets enriched among late genes\nodds < 0: targets depleted among late genes") +
  geom_hline(yintercept = -log10(0.05), linetype = 2) +
  geom_vline(xintercept = 0, linetype = 2) +
  scale_x_continuous(limits = c(-1,1)) +
  scale_y_continuous(limits = c(0, 5))

ggsave(file.path("~/Desktop/","8.pdf"), width = 5, height = 5, useDingbats = FALSE)

# basic l2fc plot
circColor <- "black"
labelColor <- "black"

mirna_df %>%
  dplyr::filter(category != "none") %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=3) +
  presTheme +
  labs(x="Log2 Fold Change Gest. Week",
       y="Log2 Fold Change GZ/CP") +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  theme(legend.position = "none") +
  geom_label(data=subset(mirna_df, mir == "miR-122-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-122-5p"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-122-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mir == "miR-106b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-106b-5p"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-106b-5p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(mirna_df, mir == "miR-19a-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="miR-19a-3p"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(mirna_df, mir == "miR-19a-3p"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

ggsave(file.path("~/Desktop/","mirna_category_mir122.pdf"), width = 11, height = 9, useDingbats = FALSE)


ggsave(file.path("~/Desktop/","mirna_category.pdf"), width = 11, height = 9, useDingbats = FALSE)

gene_df %>%
  dplyr::filter(category != "none") %>%
  ggplot(aes(x=log2FoldChange.gw, y=-log2FoldChange.gzcp, color=category)) +
  geom_point(size=1) +
  presTheme +
  labs(x="Log2 Fold Change Gest. Week",
       y="Log2 Fold Change GZ/CP") +
  scale_color_manual(values=cbPalette[c(1,2,4,5,3)]) +
  theme(legend.position = "none") +
  geom_label(data=subset(gene_df, ensg == "ENSG00000128573"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="FOXP2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000128573"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000111249"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="CUX2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000111249"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000184486"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="BRN2"), hjust=0, vjust=1.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000184486"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor) +
  geom_label(data=subset(gene_df, ensg == "ENSG00000127152"), aes(log2FoldChange.gw, -log2FoldChange.gzcp, label="CTIP2"), hjust=0, vjust=-0.5, color=labelColor) +
  geom_point(data=subset(gene_df, ensg == "ENSG00000127152"), aes(log2FoldChange.gw, -log2FoldChange.gzcp), shape=21, size=4, stroke=2, color=circColor)

ggsave(file.path("~/Desktop/","gene_category.pdf"), width = 9, height = 9, useDingbats = FALSE)

mir <- "miR-106b-5p"

# get predicted targets
predictions <- data.frame(getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2))
predictions$ENTREZID <- rownames(predictions)
# get ensembl ids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = predictions$ENTREZID,
                                 keytype = "ENTREZID",
                                 columns = c("ENSEMBL", "SYMBOL"))
# join with predictions
predictions <- left_join(predictions, rowdata, by = "ENTREZID")

# join with genes
predictions <- left_join(predictions, gene_df, by = c("ENSEMBL" = "ensg"))

predictions %<>% arrange(-log2FoldChange.gw)
write_lines(predictions$ENSEMBL[1:200], "~/Desktop/mir_106b_targets_up_late_200.txt")


mir <- "miR-19a-3p"

# get predicted targets
predictions <- data.frame(getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2))
predictions$ENTREZID <- rownames(predictions)
# get ensembl ids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = predictions$ENTREZID,
                                 keytype = "ENTREZID",
                                 columns = c("ENSEMBL", "SYMBOL"))
# join with predictions
predictions <- left_join(predictions, rowdata, by = "ENTREZID")

# join with genes
predictions <- left_join(predictions, gene_df, by = c("ENSEMBL" = "ensg"))

predictions %<>% arrange(-log2FoldChange.gw)
write_lines(predictions$ENSEMBL[1:200], "~/Desktop/mir_19a_targets_up_late_200.txt")


