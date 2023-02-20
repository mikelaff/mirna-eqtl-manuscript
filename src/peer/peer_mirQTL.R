# use peer on the mirQTL samples

library(here)
library(dplyr)
library(readr)
library(DESeq2)
library(peer)

# OUTPUT ##########
output.peer.factors.tsv <- here("results/peer/peer_factors_mirQTL.tsv")

# INPUT ###########
# small rna-seq samples used
small.qtl.samples.txt <- here("results/emmax/samples/20191101_mirQTL_RNAID_DonorID_DNAID.txt")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# Load Samples #########
# sample in analysis
samples <- read_delim(small.qtl.samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# Import Expression Data ###############################################################################################
rse <- readRDS(rse.rds)

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd <- varianceStabilizingTransformation(dds)

# subset samples
vsd <- vsd[,samples$RNAID]

rm(rse, dds)

# PEER ########
# normalized expression matrix
expr <- t(assay(vsd))

model <- PEER()

PEER_setPhenoMean(model, as.matrix(expr))

PEER_setNk(model,20)

PEER_getNk(model)

PEER_update(model)

factors <- as.data.frame(PEER_getX(model))


factors$RNAID <- rownames(expr)

write_tsv(factors, output.peer.factors.tsv)




