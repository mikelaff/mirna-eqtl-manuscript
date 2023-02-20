# build gene expression ranged summarized experiment for fetal tissue samples

# this RSE is missing phenotype information (colData)

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(SummarizedExperiment)
library(RMySQL)

# OUTPUT FILES ####################################################################################
output.rse.rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_rse_gene_counts.rds", sep="")

# INPUT FILES #####################################################################################
counts.tsv <- here("results/gene_counts/20180905_fetalTissue_gene_counts.tsv")
# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_totalRNAseq_metadata.tsv")

# # mirna rse for phenotype information
# mirna.rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# # sample data from total rna-seq
# #samples.csv <- here("data/metadata/fetal_tissue_rnaseq.csv")
# samples.tsv <- here("data/metadata/gene_expression_metadata.tsv")

# SteinLab Database ###################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- as_tibble(dbGetQuery(con, "SELECT * FROM `Donors`"))
rm(con)

# Row Data ########################################################################################
# the ensembl human database
edb <- EnsDb.Hsapiens.v86
rowranges <- genes(edb)
rm(edb)

# Column Data #####################################################################################
samples <- read_tsv(samples.metadata.tsv)

# join with donors for gestation week
samples %<>%
    left_join(dplyr::select(donors, DonorID, GestationWeek), by = "DonorID")
rm(donors)

coldata <- as.data.frame(samples)
rownames(coldata) <- coldata$RNAID

# Expression Data #################################################################################
counts <- read_tsv(counts.tsv)
# count matrix
cts <- as.matrix(dplyr::select(counts, -ENSG))
rownames(cts) <- counts$ENSG
rm(counts)

cts <- cts[,samples$RNAID]

# Create RSE ######################################################################################
# order rowranges on rownames of cts
rowranges <- rowranges[rownames(cts)]
stopifnot(all(names(rowranges) == rownames(cts)))

# order count matrix
mode(cts) <- "integer"
cts <- cts[, rownames(coldata)]
stopifnot(all(rownames(coldata) == colnames(cts)))

# build RangedSummarizedExperiment
rse <- SummarizedExperiment(assays = list(counts = cts),
                            rowRanges = rowranges,
                            colData = coldata)

# save rse to rds file
saveRDS(rse, output.rse.rds)
