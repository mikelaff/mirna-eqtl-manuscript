# Create covariate matrix with PCA components and known technical and biological batches

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)
library(RMySQL)
library(DESeq2)
library(ggplot2)

date.prefix <- format(Sys.time(), "%Y%m%d")
date.prefix <- "20200120"

# OUTPUT FILES ########################################################################################################
# Covariate matrix output directory
dir.cov.matrix <- here("results/emmax/covariates/")
dir.create(dir.cov.matrix, recursive = TRUE, showWarnings = FALSE)

# Autosomes covariate matrix
chrAutosomes.cov.mat.output <- paste0(dir.cov.matrix, date.prefix, "_chrAutosomes.mirQTLor.genotypePCA_10.expressionPCA_10.pool.purMethod.rin.sex.gestWeek.edu_index.cov")
chrAutosomes.cov.mat.labeled.output <- paste0(dir.cov.matrix, date.prefix, "_chrAutosomes.mirQTLor.columnLabeled.edu_index.cov")
# Chromosome X covariate matrix
chrX.cov.mat.output <- paste0(dir.cov.matrix, date.prefix, "_chrX.mirQTLor.genotypePCA_10.expressionPCA_10.pool.purMethod.rin.gestWeek.cov")
chrX.cov.mat.labeled.output <- paste0(dir.cov.matrix, date.prefix, "_chrX.mirQTLor.columnLabeled.cov")

# INPUT FILES #########################################################################################################
# tfam file for sample names and ordering
tfam.file <- here("results/emmax/tfiles/mirQTLor.tfam")

# genotype mds components
# mds.comps.file <- here("results/emmax/pca_mds/mirQTLor/chrAutosomes.mirQTLor.mds")

# genotype pca components
genotype.pca.comps.file <- here("results/emmax/pca_mds/mirQTLor/chrAutosomes.mirQTLor.eigenvec")

# small rna-seq samples used
small.qtl.samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# Edu attain index snp at mir 4707
edu.index.genos.raw <- "~/Downloads/edu_attain_index_genotypes.raw"

# SteinLab Database ###################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")
rm(con)
# Import Files ########################################################################################################
tfam <- read_delim(tfam.file, col_names = FALSE, delim = " ")

# mds components
# mds.comps <- read_table(mds.comps.file)

# pca components
geno.pca <- read_delim(genotype.pca.comps.file, col_names = FALSE, delim = " ")

# sample in analysis
samples <- read_delim(small.qtl.samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)

# Format Sample Data ###################################################################################################
# Association design
# Expression ~ Genotype + Kinship + 10 Genotype PCA + 10 Expression PCA + Seq Pool + Purif Method + RIN + Sex + Gest Week
# Expression ~ Genotype + Kinship + 10 Genotype PCA + 10 Expression PCA + Seq Pool + Purif Method + RIN + Gest Week

geno.pca %<>%
    select(DNAID = X2, 3:12)
colnames(geno.pca) <- c("DNAID", paste0("PC", seq(1,10)))

# filter samples
sample.metadata %<>%
    filter(RNAID %in% samples$RNAID)

# select only needed information
sample.metadata %<>%
    select(RNAID,
           DonorID,
           DNAID = VerifiedDNAID,
           Pool,
           PurificationMethod,
           RIN,
           Sex = Sex.by.XIST)

# join with donors for gestation week
sample.metadata %<>%
    left_join(select(donors, DonorID, GestationWeek), by = "DonorID")
rm(donors)

# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################
sample.metadata$PurificationMethod[which(sample.metadata$PurificationMethod == "Trizol w Glycogen")] <- "*Trizol w Glycogen(colum purified)"
# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################

# relabel sequencing_pool factor labels
sample.metadata$Pool <- factor(sample.metadata$Pool)
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool1"] <- "Pool1"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool2"] <- "Pool2"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool3"] <- "Pool3"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool4"] <- "Pool4"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool5"] <- "Pool5"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool6"] <- "Pool6"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool7"] <- "Pool7"
levels(sample.metadata$Pool)[levels(sample.metadata$Pool) == "2015-9087Pool8"] <- "Pool8"

# relabel rna_purification_method factor labels
sample.metadata$PurificationMethod <- factor(sample.metadata$PurificationMethod)
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
#levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "Trizol w Glycogen"] <- "trizol"

# format columns correctly
sample.metadata$Sex <- as.factor(sample.metadata$Sex)

# Import Expression Data ###############################################################################################
rse <- readRDS(rse.rds)

# use only miRBase expression
rse <- rse[mcols(rse)$source == "miRBase_v22",]
# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd <- varianceStabilizingTransformation(dds)

# subset samples
vsd <- vsd[,samples$RNAID]

rm(rse, dds)

# Expression PCA #######################################################################################################
expr.pca <- prcomp(t(assay(vsd)))
expr.pca <- data.frame(expr.pca$x[,1:10])
expr.pca$RNAID <- rownames(expr.pca)

# Import Edu. Attain. Index Genotypes ##################################################################################
df.genos <- read_table2(edu.index.genos.raw)
df.genos %<>%
    select(DNAID = IID, IndexGenos = `chr14:22904777:G:A_G`)

sample.metadata %<>%
    left_join(df.genos, by = "DNAID")

# Create and Export Cov Matrices #######################################################################################

# create combined dataframe for covariates
df.cov <- left_join(geno.pca, sample.metadata, by = "DNAID")
df.cov <- left_join(df.cov, expr.pca, by = "RNAID", suffix = c(".genotype", ".expression"))
df.cov <- as.data.frame(df.cov)
rownames(df.cov) <- df.cov$DNAID

# Autosome covariates
batchVars <- c("Pool", "PurificationMethod", "Sex", "RIN", "GestationWeek", "IndexGenos")

# create model matrix with 10 genotype pca components and 10 expression pca components
genoPCA <- paste0("PC", 1:10, ".genotype")
exprPCA <- paste0("PC", 1:10, ".expression")
formula.string <- paste("~",
                        paste(paste(genoPCA, collapse = " + "),
                              paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.cov))

stopifnot(all(rownames(modelmat) == tfam$X2))
# combine with first 2 columns of tfam file
outputmat <- bind_cols(tfam[,1:2], modelmat)
stopifnot(all(outputmat$X1 == tfam$X1))
stopifnot(all(outputmat$X2 == tfam$X2))
# write covariates matrix
write_delim(outputmat, chrAutosomes.cov.mat.output, col_names = FALSE)
write_delim(outputmat, chrAutosomes.cov.mat.labeled.output, col_names = TRUE)

rm(outputmat)

# chrX covariates
batchVars <- c("Pool", "PurificationMethod", "RIN", "GestationWeek")

# create model matrix with 10 genotype pca components and 10 expression pca components
formula.string <- paste("~",
                        paste(paste(genoPCA, collapse = " + "),
                              paste(exprPCA, collapse = " + "),
                              paste(batchVars, collapse = " + "),
                              sep = " + "))

modelmat <- as.data.frame(model.matrix(as.formula(formula.string), data = df.cov))

stopifnot(all(rownames(modelmat) == tfam$X2))
# combine with first 2 columns of tfam file
outputmat <- bind_cols(tfam[,1:2], modelmat)
stopifnot(all(outputmat$X1 == tfam$X1))
stopifnot(all(outputmat$X2 == tfam$X2))
# write covariates matrix
write_delim(outputmat, chrX.cov.mat.output, col_names = FALSE)
write_delim(outputmat, chrX.cov.mat.labeled.output, col_names = TRUE)
