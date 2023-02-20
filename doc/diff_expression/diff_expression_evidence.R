# compile evidence


library(here)
library(DESeq2)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(pheatmap)
library(RColorBrewer)

source(here("src/utils/lafferty_utils.R"))
prefix <- paste(format(Sys.time(), "%Y%m%d"), "diff_expression_evidence", sep="_")

# OUTPUT FILES ########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/diff_expression/pdfs/")
# write plots
write.plots <- FALSE

# data frame for compiling DE evidence
df.rds <- paste(paste(here("doc/diff_expression/rdata/"), prefix, "_df.rds", sep=""))

# INPUT FILES #########################################################################################################
# RangedSummarizedExperiment (rse) with known and novel counts. mirbaseV22 counts quantified by mirge2.0
# (mirbaseV22, friedlander, nowakowski, mirdeep2, mirge2.0)
rse.mirna.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# rse gene expression
rse.genes.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")
# evidence df
df.evidence.rds <- here("doc/target_overlap/rdata/20190507_target_overlap_fisher_df.rds")

# GRanges for mirbase and novel mirna
mirbase.granges.rds <- here("data/gtf_and_granges/20190430_mirbase_v22_mirna_granges.rds")
fried.granges.rds <- here("data/gtf_and_granges/20190430_friedlander_mirna_granges.rds")
nowa.granges.rds <- here("data/gtf_and_granges/20190430_nowakowski_mirna_granges.rds")
mirdeep.granges.rds <- here("data/gtf_and_granges/20190430_mirdeep_novel_mirna_granges.rds")
mirge.granges.rds <- here("data/gtf_and_granges/20190430_mirge_novel_mirna_granges.rds")

# single cell data from luis at UCLA showing genes upregulated in cell type clusters
genes.by.cluster.csv <- here("data/ucla_single_cell/TableS6 Cluster enriched genes.csv")

# Load GRanges Objects ################################################################################################
gr.mirbase <- readRDS(mirbase.granges.rds)
gr.fried <- readRDS(fried.granges.rds)
gr.nowa <- readRDS(nowa.granges.rds)
gr.mirdeep <- readRDS(mirdeep.granges.rds)
gr.mirge <- readRDS(mirge.granges.rds)

gr.nowa <- gr.nowa[!duplicated(gr.nowa$ID)]

# Import ##########
df <- readRDS(df.evidence.rds)

# import RSE for mirna counts
rse <- readRDS(rse.mirna.rds)
# granges for diff expressed mirnas
gr.de <- rowRanges(rse)[which(names(rowRanges(rse)) %in% df$Name)]

all(names(gr.de) == df$Name)

df.genes.by.cluster <- read_csv(genes.by.cluster.csv)
df.genes.by.cluster$Cluster <- factor(df.genes.by.cluster$Cluster)

clusters <- levels(df.genes.by.cluster$Cluster)

# Get Overlaps #########################3
# get overlaps into each of the types of mirnas
df$overlap_miRBase <- countOverlaps(gr.de, gr.mirbase[which(mcols(gr.mirbase)$type == "miRNA")])
df$overlap_fried <- countOverlaps(gr.de, gr.fried[which(mcols(gr.fried)$type != "miRNA_putative_precursor")])
df$overlap_nowa <- countOverlaps(gr.de, gr.nowa[which(mcols(gr.nowa)$type != "miRNA_putative_precursor")])
df$overlap_mirdeep <- countOverlaps(gr.de, gr.mirdeep[which(mcols(gr.mirdeep)$type != "miRNA_putative_precursor")])
df$overlap_mirge <- countOverlaps(gr.de, gr.mirge[which(mcols(gr.mirge)$type != "miRNA_putative_precursor")])


# Enriched or Not ###############
df$FT_adj_pval_gz_cp <- p.adjust(df$FT_pval_gz_cp, method = "fdr")
df$FT_adj_pval_cp_gz <- p.adjust(df$FT_pval_cp_gz, method = "fdr")
df$FT_adj_pval_mat_early <- p.adjust(df$FT_pval_mat_early, method = "fdr")
df$FT_adj_pval_mat_late <- p.adjust(df$FT_pval_mat_late, method = "fdr")
df$FT_adj_pval_End <- p.adjust(df$FT_pval_End, method = "fdr")
df$FT_adj_pval_ExCal <- p.adjust(df$FT_pval_ExCal, method = "fdr")
df$FT_adj_pval_ExDp1 <- p.adjust(df$FT_pval_ExDp1, method = "fdr")
df$FT_adj_pval_ExDp2 <- p.adjust(df$FT_pval_ExDp2, method = "fdr")
df$FT_adj_pval_ExM <- p.adjust(df$FT_pval_ExM, method = "fdr")
df$FT_adj_pval_ExN <- p.adjust(df$FT_pval_ExN, method = "fdr")
df$FT_adj_pval_InCALB2 <- p.adjust(df$FT_pval_InCALB2, method = "fdr")
df$FT_adj_pval_InSST <- p.adjust(df$FT_pval_InSST, method = "fdr")
df$FT_adj_pval_IP <- p.adjust(df$FT_pval_IP, method = "fdr")
df$FT_adj_pval_Mic <- p.adjust(df$FT_pval_Mic, method = "fdr")
df$FT_adj_pval_OPC <- p.adjust(df$FT_pval_OPC, method = "fdr")
df$FT_adj_pval_oRG <- p.adjust(df$FT_pval_oRG, method = "fdr")
df$FT_adj_pval_Per <- p.adjust(df$FT_pval_Per, method = "fdr")
df$FT_adj_pval_PgG2M <- p.adjust(df$FT_pval_PgG2M, method = "fdr")
df$FT_adj_pval_PgS <- p.adjust(df$FT_pval_PgS, method = "fdr")
df$FT_adj_pval_vRG <- p.adjust(df$FT_pval_vRG, method = "fdr")

pval.cutoff <- 0.05

df$enriched_gz_cp <- (df$FT_adj_pval_gz_cp <= pval.cutoff) & (df$FT_odds_gz_cp > 1) & (df$FT_ci_low_gz_cp > 1)
df$enriched_cp_gz <- (df$FT_adj_pval_cp_gz <= pval.cutoff) & (df$FT_odds_cp_gz > 1) & (df$FT_ci_low_cp_gz > 1)
df$enriched_mat_early <- (df$FT_adj_pval_mat_early <= pval.cutoff) & (df$FT_odds_mat_early > 1) & (df$FT_ci_low_mat_early > 1)
df$enriched_mat_late <- (df$FT_adj_pval_mat_late <= pval.cutoff) & (df$FT_odds_mat_late > 1) & (df$FT_ci_low_mat_late > 1)

df$enriched_End <- (df$FT_adj_pval_End <= pval.cutoff) & (df$FT_odds_End > 1) & (df$FT_ci_low_End > 1)
df$enriched_ExCal <- (df$FT_adj_pval_ExCal <= pval.cutoff) & (df$FT_odds_ExCal > 1) & (df$FT_ci_low_ExCal > 1)
df$enriched_ExDp1 <- (df$FT_adj_pval_ExDp1 <= pval.cutoff) & (df$FT_odds_ExDp1 > 1) & (df$FT_ci_low_ExDp1 > 1)
df$enriched_ExDp2 <- (df$FT_adj_pval_ExDp2 <= pval.cutoff) & (df$FT_odds_ExDp2 > 1) & (df$FT_ci_low_ExDp2 > 1)
df$enriched_ExM <- (df$FT_adj_pval_ExM <= pval.cutoff) & (df$FT_odds_ExM > 1) & (df$FT_ci_low_ExM > 1)
df$enriched_ExN <- (df$FT_adj_pval_ExN <= pval.cutoff) & (df$FT_odds_ExN > 1) & (df$FT_ci_low_ExN > 1)
df$enriched_InCALB2 <- (df$FT_adj_pval_InCALB2 <= pval.cutoff) & (df$FT_odds_InCALB2 > 1) & (df$FT_ci_low_InCALB2 > 1)
df$enriched_InSST <- (df$FT_adj_pval_InSST <= pval.cutoff) & (df$FT_odds_InSST > 1) & (df$FT_ci_low_InSST > 1)
df$enriched_IP <- (df$FT_adj_pval_IP <= pval.cutoff) & (df$FT_odds_IP > 1) & (df$FT_ci_low_IP > 1)
df$enriched_Mic <- (df$FT_adj_pval_Mic <= pval.cutoff) & (df$FT_odds_Mic > 1) & (df$FT_ci_low_Mic > 1)
df$enriched_OPC <- (df$FT_adj_pval_OPC <= pval.cutoff) & (df$FT_odds_OPC > 1) & (df$FT_ci_low_OPC > 1)
df$enriched_oRG <- (df$FT_adj_pval_oRG <= pval.cutoff) & (df$FT_odds_oRG > 1) & (df$FT_ci_low_oRG > 1)
df$enriched_Per <- (df$FT_adj_pval_Per <= pval.cutoff) & (df$FT_odds_Per > 1) & (df$FT_ci_low_Per > 1)
df$enriched_PgG2M <- (df$FT_adj_pval_PgG2M <= pval.cutoff) & (df$FT_odds_PgG2M > 1) & (df$FT_ci_low_PgG2M > 1)
df$enriched_PgS <- (df$FT_adj_pval_PgS <= pval.cutoff) & (df$FT_odds_PgS > 1) & (df$FT_ci_low_PgS > 1)
df$enriched_vRG <- (df$FT_adj_pval_vRG <= pval.cutoff) & (df$FT_odds_vRG > 1) & (df$FT_ci_low_vRG > 1)


df$depleted_gz_cp <- (df$FT_adj_pval_gz_cp <= pval.cutoff) & (df$FT_odds_gz_cp < 1) & (df$FT_ci_high_gz_cp < 1)
df$depleted_cp_gz <- (df$FT_adj_pval_cp_gz <= pval.cutoff) & (df$FT_odds_cp_gz < 1) & (df$FT_ci_high_cp_gz < 1)
df$depleted_mat_early <- (df$FT_adj_pval_mat_early <= pval.cutoff) & (df$FT_odds_mat_early < 1) & (df$FT_ci_high_mat_early < 1)
df$depleted_mat_late <- (df$FT_adj_pval_mat_late <= pval.cutoff) & (df$FT_odds_mat_late < 1) & (df$FT_ci_high_mat_late < 1)

df$depleted_End <- (df$FT_adj_pval_End <= pval.cutoff) & (df$FT_odds_End < 1) & (df$FT_ci_high_End < 1)
df$depleted_ExCal <- (df$FT_adj_pval_ExCal <= pval.cutoff) & (df$FT_odds_ExCal < 1) & (df$FT_ci_high_ExCal < 1)
df$depleted_ExDp1 <- (df$FT_adj_pval_ExDp1 <= pval.cutoff) & (df$FT_odds_ExDp1 < 1) & (df$FT_ci_high_ExDp1 < 1)
df$depleted_ExDp2 <- (df$FT_adj_pval_ExDp2 <= pval.cutoff) & (df$FT_odds_ExDp2 < 1) & (df$FT_ci_high_ExDp2 < 1)
df$depleted_ExM <- (df$FT_adj_pval_ExM <= pval.cutoff) & (df$FT_odds_ExM < 1) & (df$FT_ci_high_ExM < 1)
df$depleted_ExN <- (df$FT_adj_pval_ExN <= pval.cutoff) & (df$FT_odds_ExN < 1) & (df$FT_ci_high_ExN < 1)
df$depleted_InCALB2 <- (df$FT_adj_pval_InCALB2 <= pval.cutoff) & (df$FT_odds_InCALB2 < 1) & (df$FT_ci_high_InCALB2 < 1)
df$depleted_InSST <- (df$FT_adj_pval_InSST <= pval.cutoff) & (df$FT_odds_InSST < 1) & (df$FT_ci_high_InSST < 1)
df$depleted_IP <- (df$FT_adj_pval_IP <= pval.cutoff) & (df$FT_odds_IP < 1) & (df$FT_ci_high_IP < 1)
df$depleted_Mic <- (df$FT_adj_pval_Mic <= pval.cutoff) & (df$FT_odds_Mic < 1) & (df$FT_ci_high_Mic < 1)
df$depleted_OPC <- (df$FT_adj_pval_OPC <= pval.cutoff) & (df$FT_odds_OPC < 1) & (df$FT_ci_high_OPC < 1)
df$depleted_oRG <- (df$FT_adj_pval_oRG <= pval.cutoff) & (df$FT_odds_oRG < 1) & (df$FT_ci_high_oRG < 1)
df$depleted_Per <- (df$FT_adj_pval_Per <= pval.cutoff) & (df$FT_odds_Per < 1) & (df$FT_ci_high_Per < 1)
df$depleted_PgG2M <- (df$FT_adj_pval_PgG2M <= pval.cutoff) & (df$FT_odds_PgG2M < 1) & (df$FT_ci_high_PgG2M < 1)
df$depleted_PgS <- (df$FT_adj_pval_PgS <= pval.cutoff) & (df$FT_odds_PgS < 1) & (df$FT_ci_high_PgS < 1)
df$depleted_vRG <- (df$FT_adj_pval_vRG <= pval.cutoff) & (df$FT_odds_vRG < 1) & (df$FT_ci_high_vRG < 1)

write_csv(df, "~/Desktop/DE_all_miRNAs_all_data.csv")
saveRDS(df, df.rds)

# Some Plots ###############
df.mirbase <- dplyr::filter(df)

df.tmp <- df[!is.na(df$enriched_cp_gz),]

df.tmp <- select(df.tmp, 1:23, starts_with("enriched"), starts_with("depleted"))
df.tmp <- select(df.tmp, Name, everything())

write_csv(df.tmp, "~/Desktop/DE_miRBase_miRNAs_Enrichments.csv")

df.tmp <- df[df$source != "miRBase_v22",]

df.tmp <- select(df.tmp, 1:23, starts_with("overlap"))
df.tmp <- select(df.tmp, Name, everything())

write_csv(df.tmp, "~/Desktop/DE_novel_miRNAs_overlaps.csv")

df.odds <- select(df.tmp, Name, starts_with("FT_odds"))

df.odds <- select(df.odds, 1:5, 7:11, 6, 12:21)

df.enriched <- select(df.tmp, Name, starts_with("enriched"))

df.depleted <- select(df.tmp, Name, starts_with("depleted"))

colnames(df.odds) <- c("Name", "gz_cp", "cp_gz", "mat_early", "mat_late", clusters)
colnames(df.depleted) <- c("Name", "gz_cp", "cp_gz", "mat_early", "mat_late", clusters)
colnames(df.enriched) <- c("Name", "gz_cp", "cp_gz", "mat_early", "mat_late", clusters)

mat.odds <- df.odds[,2:21]
mat.depleted <- df.depleted[,2:21]
mat.enriched <- df.enriched[,2:21]

rownames(mat.odds) <- df.odds$Name
rownames(mat.depleted) <- df.depleted$Name
rownames(mat.enriched) <- df.enriched$Name



mat.odds <- log10(mat.odds)
mat.odds[mat.odds == -Inf] <- 0

mat.odds[!(mat.depleted | mat.enriched)] <- 0

mat.odds <- as.matrix(mat.odds)

mat <- (sweep(sweep(mat.odds, 2, colMins(mat.odds), "-"), 2, (colMaxs(mat.odds) - colMins(mat.odds)), "/"))


cols <- colorRampPalette(rev(brewer.pal(10,"RdBu")))(255)

pheatmap(mat = mat,
         color = cols,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_colnames = TRUE,
         show_rownames = FALSE)

