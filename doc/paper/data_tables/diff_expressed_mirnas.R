
library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(RMySQL)
library(DESeq2)

# OUTPUT FILES #########################################################################################################
# supplementary table ?: differentially expressed mirnas
output.csv <- here("doc/paper/data_tables/csv/supplementaryTableX_diff_expressed_mirnas.csv")


# INPUT FILES ##########################################################################################################
# mirQTL samples used in EMMAX association analysis
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GLOBALS ##############################################################################################################
# significance threshold
SIG.THRESHOLD <- 0.1

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

# Load Samples #########################################################################################################
# mirQTL samples
df.samples <- read_table2(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"))

# sample metadata
sample.metadata <- read_tsv(samples.metadata.tsv)
sample.metadata %<>%
    dplyr::filter(RNAID %in% df.samples$RNAID)

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


# Load Expression Data #################################################################################################
rse <- readRDS(rse.rds)

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet
dds <- DESeqDataSet(rse, design = ~1)

# subset samples
dds <- dds[,df.samples$RNAID]

rm(rse)

# Design Equation #################################################################################
# factors for design equation gest. week:
# sequencing_pool
# rna_purification_method
# rin
# rna_concentration
# gestation_week

# relabel sequencing_pool factor labels
dds$sequencing_pool <- factor(dds$sequencing_pool)
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool1"] <- "Pool1"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool2"] <- "Pool2"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool3"] <- "Pool3"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool4"] <- "Pool4"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool5"] <- "Pool5"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool6"] <- "Pool6"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool7"] <- "Pool7"
levels(dds$sequencing_pool)[levels(dds$sequencing_pool) == "2015-9087Pool8"] <- "Pool8"

# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################
dds$rna_purification_method[which(dds$rna_purification_method == "Trizol w Glycogen")] <- "*Trizol w Glycogen(colum purified)"
# THIS CODE RELABELS ONE SAMPLE !!!!!!!!!!!!!!! ##################################

# relabel rna_purification_method factor labels
dds$rna_purification_method <- factor(dds$rna_purification_method)
dds$rna_purification_method_name <- dds$rna_purification_method
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "*Trizol w Glycogen(colum purified)"] <- "trizol_col_pur"
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "miRNeasy w DNase Digestion; DL extracted"] <- "miRNeasy"
levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "miRNeasy-mini w Dnase Digestion; DL extracted"] <- "miRNeasy_mini"
#levels(dds$rna_purification_method_name)[levels(dds$rna_purification_method_name) == "Trizol w Glycogen"] <- "trizol"

# set design equation
design(dds) <- formula(~ sequencing_pool + rna_purification_method_name + rin + rna_concentration + gestation_week)

# Run DESeq2 ###########################################################################################################
dds <- DESeq(dds)

shrunkres <- lfcShrink(dds, coef="gestation_week", type="apeglm")

df.res <- as.data.frame(shrunkres)
df.res$Name <- rownames(df.res)

df.res <- as_tibble(df.res)

df.res$sig <- df.res$padj <= SIG.THRESHOLD


# rough order for plotting
#df <- df[order(df$padj.gw, decreasing = TRUE),]

# add in mcols data
rows <- data.frame(mcols(dds))
rows$Name <- rownames(rows)
rows %<>% dplyr::select(Name, ID, Alias, Derives_from, source, type, score, sequence)

df.res <- left_join(df.res, rows, by = "Name")
rm(rows)

# Format #####

df.res %>%
    select(BASE_MEAN = baseMean,
           LOG2_FOLD_CHANGE = log2FoldChange,
           LFC_SE = lfcSE,
           PVALUE = pvalue,
           PADJ = padj,
           SIGNIFICANT = sig,
           NAME = Name,
           ID,
           ALIAS = Alias,
           DERIVES_FROM = Derives_from,
           SOURCE = source,
           TYPE = type,
           SCORE = score,
           SEQUENCE = sequence) -> df.export


# Export #####
write_csv(df.export, output.csv)

df.res %>%
    filter(sig, log2FoldChange > 0) -> df.pos


