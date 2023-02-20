# genotype pca and mds for mirna-eqtl samples

library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(mikelaffr)
library(here)

# OUTPUT ##############################################################################################################
dir.pdfs <- here("doc/genotype_pca/pdfs/")

# INPUT ###############################################################################################################
# genotype mds components
mds.comps.file <- here("results/matrixEQTL/pca_mds/chrAll.QC.mirQTL.mds")

# genotype pca components
pca.comps.file <- here("results/matrixEQTL/pca_mds/chrAll.QC.mirQTL.eigenvec")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# GLOBALS #############################################################################################################

# SteinLab Database ###################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")

dbDisconnect(con)
rm(con)

# Import Files ########################################################################################################
tfam <- read_delim(tfam.file, col_names = FALSE, delim = " ")

# mds components
mds.comps <- read_table(mds.comps.file)

# pca components
pca.comps <- read_delim(pca.comps.file, col_names = FALSE, delim = " ")

# sample in analysis
samples <- read_delim(small.qtl.samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)

# Format and Create Covariate Matrix ##################################################################################
# Association design
# Expression ~ Genotype + Kinship + 10 Genotype MDS or PCA + Seq Pool + Purif Method + RIN + Sex + Gest Week

mds.comps %<>%
    select(DNAID = IID, 4:13)

pca.comps %<>%
    select(DNAID = X2, 3:12)
colnames(pca.comps) <- c("DNAID", paste0("PC", seq(1,10)))

# mds.comps %>%
#     ggplot(aes(C1,C2)) + geom_point()
#
# pca.comps %>%
#     ggplot(aes(X3,X4)) + geom_point()

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
levels(sample.metadata$PurificationMethod)[levels(sample.metadata$PurificationMethod) == "Trizol w Glycogen"] <- "trizol"

# format columns correctly
sample.metadata$Sex <- as.factor(sample.metadata$Sex)

