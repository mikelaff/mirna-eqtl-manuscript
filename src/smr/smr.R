
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(mikelaffr)

library(pheatmap)

# OUTPUT FILES ########################################################################################################
# output dir for coloc files
dir.output <- here("results/smr/")
dir.create(dir.output, showWarnings = FALSE, recursive = TRUE)

# INPUT FILES ##########################################################################################################
# Conditional eQTLs
conditional.eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# Summarized PRIMARY association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# GWAS Summary Stats: Educational Attain.
edu.sumstats.granges.rds <- here("data/gwas_datasets/educational_attainment/educational_attainment.hg38.GRanges.rds")

# Association results mRNA-eQTLs for chr14
mrna.results.chr14.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/raw/chr14_fetal_brain_mRNA-eQTL_raw.rds")

#1kg EUR bim file to know which rsids can be used in the analysis
eur.bim <- here("data/1000genomes_phase3_hg38/EUR.plink/EUR.chr14_GRCh38.uniqueRSID.bim")

# GLOBALS ##############################################################################################################
# miR-4707-3p eQTL
chr <- "chr14"
eSNP <- "chr14:22953244:A:G"
emiR <- "hsa-mir-4707_hsa-miR-4707-3p"
esnp.bp <- 22953244

# HAUS4 mRNA-eQTL
haus4.ensg <- "ENSG00000092036"

# From Nil and discussion at: https://stats.stackexchange.com/questions/337070/compute-standard-error-from-beta-p-value-sample-size-and-the-number-of-regres

getSE = function(beta, p, n, df){
    beta = as.numeric(beta) # beta estimate
    p = as.numeric(p) # pvalue
    n = as.numeric(n) # sample size
    df = as.numeric(df) # degree of freedom calculated as n - (total number of features(covariates) including the intercept)
    # add 1 for kinship matrix
    t_val <- qt(p/2, df = df) # Calculating the t-value using quantile function
    se = abs(beta)/abs(t_val) # Calculating standard error
    return(se)
}
# Import Data ##################################################################################################

# miRNA-eQTLs
df.eqtls <- readRDS(conditional.eqtls.dataframe.rds)

# Primary association data
df.eqtl.results <- as_tibble(readRDS(summarized.results.dataframe.rds))

# Filter for miR-4707-3p associations
df.eqtl.results %<>%
    filter(UniName == emiR)

max.bp <- max(df.eqtl.results$BP.hg38)
min.bp <- min(df.eqtl.results$BP.hg38)

# Educational Attain sum stats
gr.sumstats.edu <- readRDS(edu.sumstats.granges.rds)

range.overlap <- GRanges(seqnames = chr,
                         ranges = IRanges(start = min.bp, end = max.bp),
                         strand = "*")

gr.sumstats.edu <- subsetByOverlaps(gr.sumstats.edu, range.overlap)

df.sumstats.edu <- as_tibble(gr.sumstats.edu)

rm(gr.sumstats.edu, range.overlap)

df.sumstats.edu %<>%
    mutate(MAF = ifelse(EUR.freq > 0.5, 1 - EUR.freq, EUR.freq))

num.samples.gwas <- 766345

df.sumstats.edu %<>%
    mutate(PVE = (2 * (Beta^2) * MAF * (1 - MAF)) /
               ((2 * (Beta^2) * MAF * (1 - MAF)) + ((SE^2) * 2 * num.samples.gwas * MAF * (1 - MAF))))

# mRNA-eQTL
df.mrna <- readRDS(mrna.results.chr14.rds)
df.mrna %<>%
    filter(ENSG == haus4.ensg)

df.eur.bim <- read_table(eur.bim, col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))
df.eur.bim %<>%
    filter(bp >= min.bp, bp <= max.bp)

# eQTL SE and Z-score ##################################################################################################

# number of samples
num.samples.eqtl <- 212

# degrees of freedom, add 1 for secondary, 2 for tertiary, etc...
# design equation: expression ~ PC1.genotype + PC2.genotype + PC3.genotype + PC4.genotype + PC5.genotype + PC6.genotype +
# PC7.genotype + PC8.genotype + PC9.genotype + PC10.genotype + PC1.expression + PC2.expression + PC3.expression + PC4.expression +
# PC5.expression + PC6.expression + PC7.expression + PC8.expression + PC9.expression + PC10.expression + PoolPool2 + PoolPool3 +
# PoolPool4 + PoolPool5 + PoolPool6 + PoolPool7 + PoolPool8 + PurificationMethodmiRNeasy + PurificationMethodmiRNeasy_mini +
# SexM + RIN + GestationWeek + primaryDosage + secondaryDosage + tertiaryDosage
# add 1 for kinship matrix

primaryDF <- num.samples.eqtl - (32 + 1) - 1

# compute MAF
df.eqtl.results %<>%
    mutate(MAF = ALT_CTS / OBS_CT) %>%
    mutate(MAF = ifelse(MAF > 0.5, 1 - MAF, MAF))


# compute standard error
df.eqtl.results %<>%
    mutate(SE_beta = getSE(BETA, P, num.samples.eqtl, primaryDF))


# compute z-score
df.eqtl.results %<>%
    mutate(Z_score = BETA / SE_beta)

# Format mRNA-eQTL Data ################
df.mrna %<>%
    mutate(SE_beta = getSE(BETA, P, num.samples.eqtl, primaryDF))

df.mrna %<>%
    mutate(Z_score = BETA / SE_beta)

# Format GWAS Data ###########

# compute z-score
df.sumstats.edu %<>%
    mutate(Z_score = Beta / SE)

# Combine Datasets ############

df.eqtl.results %<>%
    select(SNP,
           CHR,
           BP = BP.hg38,
           REF,
           ALT,
           EFFECT.ALLELE,
           MAF,
           BETA,
           P,
           SE_beta,
           Z_score)

df.sumstats.edu %<>%
    select(RSID,
           CHR = seqnames,
           BP = start,
           A1.GWAS = A1.effect,
           A2.GWAS = A2,
           EFFECT.ALLELE = A1.effect,
           MAF,
           BETA = Beta,
           P = Pval,
           SE_beta = SE,
           Z_score)

df.combined <- left_join(df.eqtl.results, df.sumstats.edu, by = c("CHR", "BP"), suffix = c(".eQTL", ".GWAS"))
df.combined %<>%
    filter(!is.na(RSID))

# if effect alleles not the same, flip beta and z
df.combined %<>%
    mutate(BETA.GWAS.fixed = ifelse(EFFECT.ALLELE.eQTL == EFFECT.ALLELE.GWAS, BETA.GWAS, -BETA.GWAS),
           Z_score.GWAS.fixed = ifelse(EFFECT.ALLELE.eQTL == EFFECT.ALLELE.GWAS, Z_score.GWAS, -Z_score.GWAS))

# mRNA data
df.mrna %<>%
    select(CHR,
           BP,
           BETA.mRNA = BETA,
           P.mRNA = P,
           A1.mRNA = ALLELE_MINOR,
           A2.mRNA = ALLELE_MAJOR_EFFECT,
           EFFECT.ALLELE.mRNA = ALLELE_MAJOR_EFFECT,
           SE_beta.mRNA = SE_beta,
           Z_score.mRNA = Z_score)

df.combined %<>%
    left_join(df.mrna, by = c("CHR", "BP"))

df.combined %<>%
    filter(!is.na(BETA.mRNA))

# if effect alleles not the same, flip beta and z
df.combined %<>%
    mutate(BETA.mRNA.fixed = ifelse(EFFECT.ALLELE.eQTL == EFFECT.ALLELE.mRNA, BETA.mRNA, -BETA.mRNA),
           Z_score.mRNA.fixed = ifelse(EFFECT.ALLELE.eQTL == EFFECT.ALLELE.mRNA, Z_score.mRNA, -Z_score.mRNA))

# filter for rsid in the eur bim file in order to make ld matrix
df.combined %<>%
    filter(RSID %in% df.eur.bim$rsid)



# Export #########
# export snp list to make plink fileset and ld matrix

# eQTL snps for LD matrix
write_lines(df.combined$SNP, paste0(dir.output, "eqtl.snps.txt"))
# mRNA eQTL snps for LD matrix
write_lines(df.combined$SNP, paste0(dir.output, "mrna.snps.txt"))
# GWAS snps for LD matrix
write_lines(df.combined$RSID, paste0(dir.output, "gwas.snps.txt"))

# RUN PLINK EXTERNALLY


# Import LD Matrix ############

# plink bim files to get counted allele (a1 allele)
df.eqtl.bim <- read_table(paste0(dir.output, "plink_filesets/eqtl.snps.bim"), col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))
df.gwas.bim <- read_table(paste0(dir.output, "plink_filesets/gwas.snps.bim"), col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))
df.mrna.bim <- read_table(paste0(dir.output, "plink_filesets/mrna.snps.bim"), col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))

df.eqtl.freq <- read_table(paste0(dir.output, "plink_filesets/eqtl.snps.frq"))

# make .flist file
df.eqtls %>%
    mutate(GeneticDistance = 0,
           Orientation = "-",
           PathOfEsd = paste0(dir.output, emiR, ".esd"),
           Chr = "14") %>%
    select(Chr,
           ProbeID = emiR,
           GeneticDistance,
           ProbeBp = SNP.BP.hg38,
           Gene = emiR,
           Orientation,
           PathOfEsd) %>%
    filter(ProbeID == emiR) %>%
    write_tsv(file = paste0(dir.output, "my.flist"))

# make .esd file with eQTL A1 as the effect allele
df.combined %>%
    select(SNP, RSID, Beta = BETA.eQTL, se = SE_beta.eQTL, p = P.eQTL, effect.allele = EFFECT.ALLELE.eQTL, Bp = BP) %>%
    left_join(df.eqtl.freq, by = "SNP") -> df.esd.eqtl.maf

df.esd.eqtl.maf %<>%
    mutate(Beta.fixed = ifelse(effect.allele == A1, Beta, -Beta))

df.esd.eqtl.maf %>%
    select(Chr = CHR,
           SNP,
           Bp,
           A1,
           A2,
           Freq = MAF,
           Beta = Beta.fixed,
           se,
           p) %>%
    write_tsv(file = paste0(dir.output, "hsa-mir-4707_hsa-miR-4707-3p.esd"))

# make .ma file for GWAS data, A1 is the effect allele but taken from the eQTL study

df.combined %>%
    select(SNP, RSID, b = BETA.GWAS, se = SE_beta.GWAS, p = P.GWAS, effect.allele = EFFECT.ALLELE.GWAS) %>%
    mutate(n = num.samples.gwas) %>%
    left_join(df.eqtl.freq, by = "SNP") -> df.ma.eqtl.maf

df.ma.eqtl.maf %<>%
    mutate(Beta.fixed = ifelse(effect.allele == A1, b, -b))

df.ma.eqtl.maf %>%
    select(SNP,
           A1,
           A2,
           freq = MAF,
           b = Beta.fixed,
           se,
           p,
           n) %>%
    write_tsv(file = paste0(dir.output, "mygwas.ma"))


# make .esd file with GWAS A1 as the effect allele

df.gwas.freq <- read_table(paste0(dir.output, "plink_filesets/gwas.snps.frq"))

df.combined %>%
    select(SNP, RSID, Beta = BETA.eQTL, se = SE_beta.eQTL, p = P.eQTL, effect.allele = EFFECT.ALLELE.eQTL, Bp = BP) %>%
    left_join(df.gwas.freq, by = c("RSID" = "SNP")) -> df.esd.eqtl.maf

df.esd.eqtl.maf %<>%
    mutate(Beta.fixed = ifelse(effect.allele == A1, Beta, -Beta))

df.esd.eqtl.maf %>%
    select(Chr = CHR,
           SNP = RSID,
           Bp,
           A1,
           A2,
           Freq = MAF,
           Beta = Beta.fixed,
           se,
           p) %>%
    write_tsv(file = paste0(dir.output, "hsa-mir-4707_hsa-miR-4707-3p.esd"))

# make .ma file for GWAS data, A1 is the effect allele but taken from the GWAS study
df.combined %>%
    select(SNP, RSID, b = BETA.GWAS, se = SE_beta.GWAS, p = P.GWAS, effect.allele = EFFECT.ALLELE.GWAS) %>%
    mutate(n = num.samples.gwas) %>%
    left_join(df.gwas.freq, by = c("RSID" = "SNP")) -> df.ma.eqtl.maf

df.ma.eqtl.maf %<>%
    mutate(Beta.fixed = ifelse(effect.allele == A1, b, -b))

df.ma.eqtl.maf %>%
    select(SNP = RSID,
           A1,
           A2,
           freq = MAF,
           b = Beta.fixed,
           se,
           p,
           n) %>%
    write_tsv(file = paste0(dir.output, "mygwas.ma"))


# Results ###########

df.smr2 <- read_tsv(paste0(dir.output, "mysmr.smr"))



