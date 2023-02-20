
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(AnnotationHub)
#library(Gviz)
#library(gridExtra)
#library(biomaRt)
#library(psych)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/conditional_eqtls/pdfs/")
#dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# primary association results
primary.results.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")
# secondary association results
secondary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/secondary/association_results/compiled/20200120_mirQTLor_secondary_variants_dataFrame.rds")
# tertiary association results
tertiary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/tertiary/association_results/compiled/20200120_mirQTLor_tertiary_variants_dataFrame.rds")
# quarternary association results
quarternary.results.rds <- here("results/conditional_eqtls/20200120_mirQTLor/quarternary/association_results/compiled/20200120_mirQTLor_quarternary_variants_dataFrame.rds")

# primary eqtls
primary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/primary/20200120_mirQTLor_primary_eQTLs_dataFrame.rds")
# secondary eqtls
secondary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/secondary/20200120_mirQTLor_secondary_eQTLs_dataFrame.rds")
# tertiary eqtls
tertiary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor/tertiary/20200120_mirQTLor_tertiary_eQTLs_dataFrame.rds")

# nominal p-value from eigenMT-BH proceedure
nom.p.value.txt <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eigenMT-BH_nomPvalue.txt")

# OLD, clumped eigenMT-BH thresholded eqtls
clumped.eqtls.rds <- here("results/eigenmt/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")

# GLOBALS ##############################################################################################################
FDR.THRESHOLD <- as.numeric(read_lines(nom.p.value.txt))

# Summary #######################################################################################################

df.clumped.eqtls <- read_rds(clumped.eqtls.rds)

print(paste("Number of eQTLs (clumped anlaysis):", nrow(df.clumped.eqtls)))
print(paste("Number of emiRs (clumped analysis):", length(unique(df.clumped.eqtls$emiR))))
print(paste("Number of eSNPs (clumped analysis):", length(unique(df.clumped.eqtls$eSNP))))

df.primary.eqlts <- read_rds(primary.eqtls.rds)
df.secondary.eqlts <- read_rds(secondary.eqtls.rds)
df.tertiary.eqlts <- read_rds(tertiary.eqtls.rds)

df.primary.eqlts$degree <- "primary"
df.secondary.eqlts$degree <- "secondary"
df.tertiary.eqlts$degree <- "tertiary"

df.cond.eqlts <- bind_rows(df.primary.eqlts,
                           df.secondary.eqlts,
                           df.tertiary.eqlts)

print(paste("Number of eQTLs (conditional anlaysis):", nrow(df.cond.eqlts)))
print(paste("Number of emiRs (conditional analysis):", length(unique(df.cond.eqlts$UniName))))
print(paste("Number of eSNPs (conditional analysis):", length(unique(df.cond.eqlts$SNP))))

print(paste("Number of primary eQTLs:", sum(df.cond.eqlts$degree == "primary")))
print(paste("Number of secondary eQTLs:", sum(df.cond.eqlts$degree == "secondary")))
print(paste("Number of tertiary eQTLs:", sum(df.cond.eqlts$degree == "tertiary")))




