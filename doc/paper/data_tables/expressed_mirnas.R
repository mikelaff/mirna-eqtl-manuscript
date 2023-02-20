
library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(DESeq2)

# OUTPUT FILES #########################################################################################################
# supplementary table 1: expressed mirnas
output.csv <- here("doc/paper/data_tables/csv/supplementaryTable1_expressed_mirnas.csv")


# INPUT FILES ##########################################################################################################
# mirna expression used in eQTL analysis
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################

# Import #####
vsd <- read_rds(vsd.rds)

df.rowranges <- as_tibble(rowRanges(vsd))

#df.dups <- df.rowranges[duplicated(paste0(df.rowranges$seqnames, df.rowranges$start)) | duplicated(paste0(df.rowranges$seqnames, df.rowranges$start), fromLast = TRUE), ]

#sum(duplicated(df.rowranges$uniqueName))

df.rowranges$meanVSTexpression <- rowMeans(assay(vsd, 1))

# format columns
df.rowranges %>%
    select(UNIQUE_NAME = uniqueName,
           ID,
           ALIAS = Alias,
           NAME = Name,
           DERIVES_FROM = Derives_from,
           DERIVES_FROM_NAME = Derives_from_name,
           SOURCE = source,
           TYPE = type,
           SCORE = score,
           CHR = seqnames,
           START_hg38 = start,
           END_hg38 = end,
           WIDTH = width,
           STRAND = strand,
           SEQUENCE = sequence,
           MEAN_VST_EXPRESSION = meanVSTexpression) -> df.export

df.export %<>%
    arrange(SOURCE, UNIQUE_NAME)

# Export #####
write_csv(df.export, output.csv)




