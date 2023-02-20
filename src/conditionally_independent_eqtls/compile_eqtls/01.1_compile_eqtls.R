
library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
# library(DESeq2)
# library(GenomicRanges)
# library(RColorBrewer)

# OUTPUT ###############################################################################################################
# compilation of conditionally independent eQTLs at different nom P value thresholds
df.compiled.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# index snps list
index.snps.txt <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/index.snps.txt")

# INPUT ################################################################################################################
# eQTLs eigenMT-BH 5%
eqtls.eigenMT.fdr5percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_eQTLs_dataFrame.rds")
nomPvalue.eigenMT.fdr5percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor/compiled/20200120_mirQTLor_conditional_nomPvalue.txt")

# eQTLs eigenMT-BH 10%
eqtls.eigenMT.fdr10percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/20200120_mirQTLor_fdr10percent_conditional_eQTLs_dataFrame.rds")
nomPvalue.eigenMT.fdr10percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/compiled/20200120_mirQTLor_fdr10percent_conditional_nomPvalue.txt")

# eQTLs FDR 5%, no eigenMT
eqtls.fdr5percent.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_conditional_eQTLs_dataFrame.rds")
nomPvalue.fdr5percent.txt <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/compiled/20200120_mirQTLor_fdr5percent_noEigenMT_conditional_nomPvalue.txt")


# GLOBALS ##############################################################################################################


# Import Results ############################################################################################

# import eQTLs
df.eqtls.eigenMT.fdr5 <- read_rds(eqtls.eigenMT.fdr5percent.rds)
df.eqtls.eigenMT.fdr10 <- read_rds(eqtls.eigenMT.fdr10percent.rds)
df.eqtls.fdr5 <- read_rds(eqtls.fdr5percent.rds)

sum(duplicated(df.eqtls.fdr5$eQTL))

all(df.eqtls.eigenMT.fdr5$eQTL %in% df.eqtls.eigenMT.fdr10$eQTL)

all(df.eqtls.eigenMT.fdr5$eQTL %in% df.eqtls.fdr5$eQTL)

all(df.eqtls.eigenMT.fdr10$eQTL %in% df.eqtls.fdr5$eQTL)

# import nom p value thresholds
nomPvalue.eigenMT.fdr5percent <- read_lines(nomPvalue.eigenMT.fdr5percent.txt)
nomPvalue.eigenMT.fdr10percent <- read_lines(nomPvalue.eigenMT.fdr10percent.txt)
nomPvalue.fdr5percent <- read_lines(nomPvalue.fdr5percent.txt)


df.eqtls.eigenMT.fdr5$SIGNIFICANCE <- "eigenMT_fdr5percent"
df.eqtls.eigenMT.fdr5$NOM.P.VALUE.THRESH <- nomPvalue.eigenMT.fdr5percent

df.eqtls.eigenMT.fdr10$SIGNIFICANCE <- "eiganMT_fdr10percent"
df.eqtls.eigenMT.fdr10$NOM.P.VALUE.THRESH <- nomPvalue.eigenMT.fdr10percent

df.eqtls.fdr5$SIGNIFICANCE <- "fdr5percent"
df.eqtls.fdr5$NOM.P.VALUE.THRESH <- nomPvalue.fdr5percent

df.eqtls.eigenMT.fdr5 %>%
    bind_rows(filter(df.eqtls.eigenMT.fdr10, ! eQTL %in% df.eqtls.eigenMT.fdr5$eQTL)) -> df.eqtls

df.eqtls %<>%
    bind_rows(filter(df.eqtls.fdr5, ! eQTL %in% df.eqtls$eQTL))

sum(duplicated(df.eqtls$eQTL))

# Export ###############################################################################
write_rds(df.eqtls, df.compiled.eqtls.rds)

write_lines(unique(df.eqtls$eSNP), index.snps.txt)

# df.eqtls %>%
#     group_by(SIGNIFICANCE, DEGREE) %>%
#     summarise(count = n())


