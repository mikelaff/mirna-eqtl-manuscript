
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(DESeq2)
library(ggplot2)
library(ggpubr)
#library(limma)
library(RMySQL)
library(reshape2)
library(mikelaffr)

# FIGURE ###############################################################################################################
# figure 1: mirna target expression at high and low gestation weeks
mir92b.output.pdf <- paste0(here("doc/paper/figure1/pdfs/"), "figure1_target_expression_92b.pdf")

mir124.output.pdf <- paste0(here("doc/paper/figure1/pdfs/"), "figure1_target_expression_124.pdf")

# OUTPUT FILES #########################################################################################################
# output directory for pdf files
dir.pdf <- here("doc/paper/figure1/pdfs/")


# INPUT FILES ##########################################################################################################
# mirQTL samples used in EMMAX association analysis
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# ranged summarized experiment with gene expression counts
rse.rds <- here("results/rdata_files/20190505_rse_gene_counts.rds")

# GLOBALS ##############################################################################################################

mir124genes <- c("ENSG00000067560", "ENSG00000011304", "ENSG00000136518", "ENSG00000185591")
mir9genes <- c("ENSG00000134250", "ENSG00000114315", "ENSG00000115738", "ENSG00000125249", "ENSG00000131711")
mir92bgenes <- c("ENSG00000101665", "ENSG00000165699", "ENSG00000132326", "ENSG00000129757")

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

vsd <- vst(dds)

# subset samples
vsd <- vsd[,vsd$gestation_week == 14 | vsd$gestation_week == 20]
dds <- dds[,dds$gestation_week == 14 | dds$gestation_week == 20]

#rm(rse, dds)


# Plot #################################################################################################################

samples <- as_tibble(colData(dds))

mir92bgenes.names <- sapply(mir92bgenes, function(x) rowData(dds)$symbol[which(rowData(dds)$gene_id == x)])


df.92b.targets <- data.frame(t(assay(vsd)[mir92bgenes,]))
#df.92b.targets <- data.frame(t(assay(normTransform(dds[mir92bgenes,]))))

colnames(df.92b.targets) <- mir92bgenes.names[colnames(df.92b.targets)]
df.92b.targets$rnaid <- rownames(df.92b.targets)

df.92b.targets <- as_tibble(df.92b.targets)



df.92b <- as_tibble(melt(df.92b.targets, measure.vars = mir92bgenes.names, variable.name = "expression"))

df.92b %<>%
    left_join(samples, by = "rnaid")

df.92b %<>%
    mutate(gestation_week = factor(gestation_week))

df.92b %>%
    ggplot(aes(x=expression, y=value, color=gestation_week)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.2) +
    plotTheme("figure") +
    theme(legend.position = "none") +
    scale_color_manual(values = c(paperRed, paperBlue)) +
    labs(y = "mRNA VST Normalized Expression") +
    #ylim(4.5,12.5) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(mir92b.output.pdf, height = 1.2, width = 2.2)



mir124genes.names <- sapply(mir124genes, function(x) rowData(dds)$symbol[which(rowData(dds)$gene_id == x)])


df.124.targets <- data.frame(t(assay(vsd)[mir124genes,]))
#df.92b.targets <- data.frame(t(assay(normTransform(dds[mir92bgenes,]))))

colnames(df.124.targets) <- mir124genes.names[colnames(df.124.targets)]
df.124.targets$rnaid <- rownames(df.124.targets)

df.124.targets <- as_tibble(df.124.targets)



df.124 <- as_tibble(melt(df.124.targets, measure.vars = mir124genes.names, variable.name = "expression"))

df.124 %<>%
    left_join(samples, by = "rnaid")

df.124 %<>%
    mutate(gestation_week = factor(gestation_week))

df.124 %>%
    ggplot(aes(x=expression, y=value, color=gestation_week)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.2) +
    plotTheme("figure") +
    theme(legend.position = "none") +
    scale_color_manual(values = c(paperRed, paperBlue)) +
    labs(y = "mRNA VST Normalized Expression") +
    #ylim(6,12.5) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(mir124.output.pdf, height = 1.2, width = 2.2)



