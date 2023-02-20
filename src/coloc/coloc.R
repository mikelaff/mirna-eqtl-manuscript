
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(tidyr)
library(GenomicRanges)
library(coloc)
library(mikelaffr)

library(pheatmap)

# OUTPUT FILES ########################################################################################################
# output dir for coloc files
dir.output <- here("results/coloc/")
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

idx.center <- which(df.combined$BP > esnp.bp - 100000 & df.combined$BP < esnp.bp + 100000)

# df.combined %>%
#     filter(BP > 22750000, BP < 23250000, P.eQTL < 0.0005) %>%
#     pivot_longer(cols = c(BETA.eQTL, BETA.GWAS, BETA.GWAS.fixed, BETA.mRNA, BETA.mRNA.fixed), names_to = "dataset", values_to = "beta") %>%
#     filter(dataset %in% c("BETA.eQTL", "BETA.GWAS", "BETA.GWAS.fixed")) %>%
#     ggplot(aes(x = BP, y = beta, color = dataset)) +
#     geom_line()

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

# check that the bim files have the same allele order as the combined dataframe
stopifnot(all(df.combined$SNP == df.eqtl.bim$rsid))
stopifnot(all(df.combined$RSID == df.gwas.bim$rsid))
stopifnot(all(df.combined$SNP == df.mrna.bim$rsid))

# correction vector to adjust for different effect allele in each bim file
vec.eqtl.correction <- ifelse(df.combined$EFFECT.ALLELE.eQTL == df.eqtl.bim$a1, 1, -1)
vec.gwas.correction <- ifelse(df.combined$EFFECT.ALLELE.eQTL == df.gwas.bim$a1, 1, -1)
vec.mrna.correction <- ifelse(df.combined$EFFECT.ALLELE.eQTL == df.mrna.bim$a1, 1, -1)

# read in ld matrix
mat.eqtl.ld <- as.matrix(read_table(paste0(dir.output, "plink_filesets/eqtl.snps.ld"), col_names = FALSE))
mat.gwas.ld <- as.matrix(read_table(paste0(dir.output, "plink_filesets/gwas.snps.ld"), col_names = FALSE))
mat.mrna.ld <- as.matrix(read_table(paste0(dir.output, "plink_filesets/mrna.snps.ld"), col_names = FALSE))

rownames(mat.eqtl.ld) <- df.combined$RSID
colnames(mat.eqtl.ld) <- df.combined$RSID

rownames(mat.gwas.ld) <- df.combined$RSID
colnames(mat.gwas.ld) <- df.combined$RSID

rownames(mat.mrna.ld) <- df.combined$RSID
colnames(mat.mrna.ld) <- df.combined$RSID

# corrected ld matrix
mat.eqtl.ld.corrected <- t(t(mat.eqtl.ld) * vec.eqtl.correction) * vec.eqtl.correction
mat.gwas.ld.corrected <- t(t(mat.gwas.ld) * vec.gwas.correction) * vec.gwas.correction
mat.mrna.ld.corrected <- t(t(mat.mrna.ld) * vec.mrna.correction) * vec.mrna.correction

rownames(mat.eqtl.ld.corrected) <- df.combined$RSID
colnames(mat.eqtl.ld.corrected) <- df.combined$RSID

rownames(mat.gwas.ld.corrected) <- df.combined$RSID
colnames(mat.gwas.ld.corrected) <- df.combined$RSID

rownames(mat.mrna.ld.corrected) <- df.combined$RSID
colnames(mat.mrna.ld.corrected) <- df.combined$RSID

# write.table(mat.eqtl.ld.corrected, file = paste0(dir.output, "eqtl.ld.ld.corrected"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
# write.table(mat.gwas.ld.corrected, file = paste0(dir.output, "gwas.ld.ld.corrected"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")

# visualize before and after ld matrix
# pheatmap(1 * (mat.eqtl.ld[1:200,1:200] > 0), cluster_rows = FALSE, cluster_cols = FALSE)
# pheatmap(1 * (mat.gwas.ld[1:200,1:200] > 0), cluster_rows = FALSE, cluster_cols = FALSE)
#
# pheatmap(1 * (mat.eqtl.ld.corrected[1:200,1:200] > 0), cluster_rows = FALSE, cluster_cols = FALSE)
# pheatmap(1 * (mat.gwas.ld.corrected[1:200,1:200] > 0), cluster_rows = FALSE, cluster_cols = FALSE)

# Coloc ##############

# Build Coloc Datasets ###########

# variance of beta
df.combined %<>%
    mutate(varbeta.eQTL = SE_beta.eQTL^2,
           varbeta.GWAS = SE_beta.GWAS^2,
           varbeta.mRNA = SE_beta.mRNA^2)

# miRNA eQTL dataset
coloc.data.eqtl <- list(beta = setNames(df.combined$BETA.eQTL, df.combined$RSID),
                        varbeta = setNames(df.combined$varbeta.eQTL, df.combined$RSID),
                        snp = df.combined$RSID,
                        position = df.combined$BP,
                        type = "quant",
                        N = num.samples.eqtl,
                        MAF = setNames(df.combined$MAF.eQTL, df.combined$RSID),
                        LD = mat.eqtl.ld.corrected)

coloc.data.eqtl <- list(beta = setNames(df.combined$BETA.eQTL[idx.center], df.combined$RSID[idx.center]),
                        varbeta = setNames(df.combined$varbeta.eQTL[idx.center], df.combined$RSID[idx.center]),
                        snp = df.combined$RSID[idx.center],
                        position = df.combined$BP[idx.center],
                        type = "quant",
                        N = num.samples.eqtl,
                        MAF = setNames(df.combined$MAF.eQTL[idx.center], df.combined$RSID[idx.center]),
                        LD = mat.eqtl.ld.corrected[idx.center,idx.center])

str(coloc.data.eqtl)

check_dataset(coloc.data.eqtl)

plot_dataset(coloc.data.eqtl)

check_alignment(coloc.data.eqtl)

# GWAS dataset
coloc.data.gwas <- list(beta = setNames(df.combined$BETA.GWAS.fixed, df.combined$RSID),
                        varbeta = setNames(df.combined$varbeta.GWAS, df.combined$RSID),
                        snp = df.combined$RSID,
                        position = df.combined$BP,
                        type = "quant",
                        N = num.samples.gwas,
                        MAF = setNames(df.combined$MAF.GWAS, df.combined$RSID),
                        LD = mat.gwas.ld.corrected)

coloc.data.gwas <- list(beta = setNames(df.combined$BETA.GWAS.fixed[idx.center], df.combined$RSID[idx.center]),
                        varbeta = setNames(df.combined$varbeta.GWAS[idx.center], df.combined$RSID[idx.center]),
                        snp = df.combined$RSID[idx.center],
                        position = df.combined$BP[idx.center],
                        type = "quant",
                        N = num.samples.gwas,
                        MAF = setNames(df.combined$MAF.GWAS[idx.center], df.combined$RSID[idx.center]),
                        LD = mat.gwas.ld.corrected[idx.center,idx.center])

str(coloc.data.gwas)

check_dataset(coloc.data.gwas)

plot_dataset(coloc.data.gwas)

check_alignment(coloc.data.gwas)


# mRNA eQTL dataset
coloc.data.mrna <- list(beta = setNames(df.combined$BETA.mRNA.fixed, df.combined$RSID),
                        varbeta = setNames(df.combined$varbeta.mRNA, df.combined$RSID),
                        snp = df.combined$RSID,
                        position = df.combined$BP,
                        type = "quant",
                        N = num.samples.eqtl,
                        MAF = setNames(df.combined$MAF.eQTL, df.combined$RSID),
                        LD = mat.mrna.ld.corrected)

coloc.data.mrna <- list(beta = setNames(df.combined$BETA.mRNA.fixed[idx.center], df.combined$RSID[idx.center]),
                        varbeta = setNames(df.combined$varbeta.mRNA[idx.center], df.combined$RSID[idx.center]),
                        snp = df.combined$RSID[idx.center],
                        position = df.combined$BP[idx.center],
                        type = "quant",
                        N = num.samples.eqtl,
                        MAF = setNames(df.combined$MAF.eQTL[idx.center], df.combined$RSID[idx.center]),
                        LD = mat.mrna.ld.corrected[idx.center,idx.center])

str(coloc.data.mrna)

check_dataset(coloc.data.mrna)

plot_dataset(coloc.data.mrna)

check_alignment(coloc.data.mrna)


# data(coloc_test_data)
# attach(coloc_test_data)
#
#
# str(D1)
# str(D1$beta)
#
#
#
#
# str(D1$LD)
# str(coloc.data.eqtl$LD)


# Fine Mapping ########


eqt.res <- finemap.abf(coloc.data.eqtl)

head(eqt.res)
tail(eqt.res)

finemap.signals(coloc.data.eqtl, method = "mask")
finemap.signals(coloc.data.eqtl, method = "cond")
finemap.signals(coloc.data.eqtl, method = "single")

finemap.signals(coloc.data.gwas, method = "mask")
finemap.signals(coloc.data.gwas, method = "cond")
finemap.signals(coloc.data.gwas, method = "single")

# Colocalization

res.eqtl.to.gwas <- coloc.abf(dataset1 = coloc.data.eqtl,
                              dataset2 = coloc.data.gwas)

print(res.eqtl.to.gwas)


res.eqtl.to.mrna <- coloc.abf(dataset1 = coloc.data.eqtl,
                              dataset2 = coloc.data.mrna)

print(res.eqtl.to.mrna)

sensitivity(res.eqtl.to.gwas, rule = "H4 > 0.5")

sensitivity(res.eqtl.to.mrna, rule = "H4 > 0.5")


res.mrna.to.gwas <- coloc.abf(dataset1 = coloc.data.mrna,
                              dataset2 = coloc.data.gwas)

print(res.mrna.to.gwas)

sensitivity(res.mrna.to.gwas, rule = "H4 > 0.5")





# SUSIE ################


susie.eqtl <- runsusie(coloc.data.eqtl)
summary(susie.eqtl)

susie.gwas <- runsusie(coloc.data.gwas)
summary(susie.gwas)

susie.mrna <- runsusie(coloc.data.mrna)
summary(susie.mrna)

susie.res <- coloc.susie(susie.eqtl, susie.gwas)
susie.res$summary

# P-P Plot #################

rsid.index.snp <- "rs4981455"
rsid.index.snp.gwas <- "rs1043209"

df.eqtl.ld.to.index <- read_table(paste0(dir.output, "plink_filesets/eqtl.snps.rs4981455.ld"))
df.eqtl.ld.to.index %<>%
    select(SNP = SNP_B,
           R2_to_rs4981455_eqtl = R2)

df.eqtl.ld.to.index2 <- read_table(paste0(dir.output, "plink_filesets/eqtl.snps.rs1043209.ld"))
df.eqtl.ld.to.index2 %<>%
    select(SNP = SNP_B,
           R2_to_rs1043209_eqtl = R2)

df.gwas.ld.to.index <- read_table(paste0(dir.output, "plink_filesets/gwas.snps.rs4981455.ld"))
df.gwas.ld.to.index %<>%
    select(RSID = SNP_B,
           R2_to_rs4981455_gwas = R2)

df.gwas.ld.to.index2 <- read_table(paste0(dir.output, "plink_filesets/gwas.snps.rs1043209.ld"))
df.gwas.ld.to.index2 %<>%
    select(RSID = SNP_B,
           R2_to_rs1043209_gwas = R2)

df.combined %<>%
    left_join(df.eqtl.ld.to.index, by = "SNP")

df.combined %<>%
    left_join(df.eqtl.ld.to.index2, by = "SNP")

df.combined %<>%
    left_join(df.gwas.ld.to.index, by = "RSID")

df.combined %<>%
    left_join(df.gwas.ld.to.index2, by = "RSID")

# dot color based on LD
df.combined$color_to_rs4981455_eqtl <- "navy"
df.combined$color_to_rs4981455_eqtl[which(df.combined$R2_to_rs4981455_eqtl >= 0.8)] <- "red"
df.combined$color_to_rs4981455_eqtl[which(df.combined$R2_to_rs4981455_eqtl >= 0.6 & df.combined$R2_to_rs4981455_eqtl < 0.8)] <- "orange"
df.combined$color_to_rs4981455_eqtl[which(df.combined$R2_to_rs4981455_eqtl >= 0.4 & df.combined$R2_to_rs4981455_eqtl < 0.6)] <- "green"
df.combined$color_to_rs4981455_eqtl[which(df.combined$R2_to_rs4981455_eqtl >= 0.2 & df.combined$R2_to_rs4981455_eqtl < 0.4)] <- "lightblue"
df.combined$color_to_rs4981455_eqtl[which(df.combined$RSID == rsid.index.snp)] <- "purple"

df.combined$color_to_rs4981455_gwas <- "navy"
df.combined$color_to_rs4981455_gwas[which(df.combined$R2_to_rs4981455_gwas >= 0.8)] <- "red"
df.combined$color_to_rs4981455_gwas[which(df.combined$R2_to_rs4981455_gwas >= 0.6 & df.combined$R2_to_rs4981455_gwas < 0.8)] <- "orange"
df.combined$color_to_rs4981455_gwas[which(df.combined$R2_to_rs4981455_gwas >= 0.4 & df.combined$R2_to_rs4981455_gwas < 0.6)] <- "green"
df.combined$color_to_rs4981455_gwas[which(df.combined$R2_to_rs4981455_gwas >= 0.2 & df.combined$R2_to_rs4981455_gwas < 0.4)] <- "lightblue"
df.combined$color_to_rs4981455_gwas[which(df.combined$RSID == rsid.index.snp)] <- "purple"

df.combined$color_to_rs1043209_eqtl <- "navy"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.8)] <- "red"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.6 & df.combined$R2_to_rs1043209_eqtl < 0.8)] <- "orange"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.4 & df.combined$R2_to_rs1043209_eqtl < 0.6)] <- "green"
df.combined$color_to_rs1043209_eqtl[which(df.combined$R2_to_rs1043209_eqtl >= 0.2 & df.combined$R2_to_rs1043209_eqtl < 0.4)] <- "lightblue"
df.combined$color_to_rs1043209_eqtl[which(df.combined$RSID == rsid.index.snp.gwas)] <- "purple"

df.combined$color_to_rs1043209_gwas <- "navy"
df.combined$color_to_rs1043209_gwas[which(df.combined$R2_to_rs1043209_gwas >= 0.8)] <- "red"
df.combined$color_to_rs1043209_gwas[which(df.combined$R2_to_rs1043209_gwas >= 0.6 & df.combined$R2_to_rs1043209_gwas < 0.8)] <- "orange"
df.combined$color_to_rs1043209_gwas[which(df.combined$R2_to_rs1043209_gwas >= 0.4 & df.combined$R2_to_rs1043209_gwas < 0.6)] <- "green"
df.combined$color_to_rs1043209_gwas[which(df.combined$R2_to_rs1043209_gwas >= 0.2 & df.combined$R2_to_rs1043209_gwas < 0.4)] <- "lightblue"
df.combined$color_to_rs1043209_gwas[which(df.combined$RSID == rsid.index.snp.gwas)] <- "purple"



df.combined %>%
    filter(BP > esnp.bp - 100000 & BP < esnp.bp + 100000) %>%
    mutate(color_to_rs1043209_eqtl = factor(color_to_rs1043209_eqtl, levels = c("purple", "red", "orange", "green", "lightblue", "navy"), ordered = TRUE)) %>%
    arrange(color_to_rs1043209_eqtl) %>%
    ggplot(aes(x = -log10(P.eQTL), y = -log10(P.GWAS), color = color_to_rs1043209_eqtl)) +
    geom_point() +
    scale_color_manual(values = c("purple", "red", "orange", "green", "lightblue", "navy"))









df.combined %>%
    filter(BP > esnp.bp - 100000 & BP < esnp.bp + 100000) -> df.tmp

plot(x = -log10(df.tmp$P.eQTL), y = -log10(df.tmp$P.GWAS), col = df.tmp$color_to_rs4981455_eqtl, pch = 19, main = "eQTL LD")


plot(x = -log10(df.tmp$P.eQTL), y = -log10(df.tmp$P.GWAS), col = df.tmp$color_to_rs4981455_gwas, pch = 19, main = "GWAS LD (1KG EUR)")



plot(x = -log10(df.tmp$P.eQTL), y = -log10(df.tmp$P.GWAS), col = df.tmp$color_to_rs1043209_eqtl, pch = 19, main = "eQTL LD")


plot(x = -log10(df.tmp$P.eQTL), y = -log10(df.tmp$P.GWAS), col = df.tmp$color_to_rs1043209_gwas, pch = 19, main = "GWAS LD (1KG EUR)")



df.results.secondary <- as_tibble(readRDS(here("results/conditional_eqtls/20200120_mirQTLor/secondary/association_results/compiled/20200120_mirQTLor_secondary_variants_dataFrame.rds")))
df.results.secondary %<>%
    filter(UniName == "hsa-mir-4707_hsa-miR-4707-3p")


df.results.secondary %>%
    ggplot(aes(x = BP.hg38, y = -log10(P))) +
    geom_point()


df.eqtl.results %<>%
    left_join(select(df.combined, SNP, color_to_rs4981455_eqtl), by = "SNP")

df.eqtl.results$color_to_rs4981455_eqtl[is.na(df.eqtl.results$color_to_rs4981455_eqtl)] <- "navy"


df.results.secondary %<>%
    left_join(select(df.combined, SNP, color_to_rs4981455_eqtl), by = "SNP")

df.results.secondary$color_to_rs4981455_eqtl[is.na(df.results.secondary$color_to_rs4981455_eqtl)] <- "navy"


df.eqtl.results %>%
    mutate(color_to_rs4981455_eqtl = factor(color_to_rs4981455_eqtl, levels = c("purple", "red", "orange", "green", "lightblue", "navy"), ordered = TRUE)) %>%
    ggplot(aes(x = BP, y = -log10(P), color = color_to_rs4981455_eqtl)) +
    geom_point() +
    scale_color_manual(values = c("purple", "red", "orange", "green", "lightblue", "navy"))

df.results.secondary %>%
    mutate(color_to_rs4981455_eqtl = factor(color_to_rs4981455_eqtl, levels = c("purple", "red", "orange", "green", "lightblue", "navy"), ordered = TRUE)) %>%
    ggplot(aes(x = BP.hg38, y = -log10(P))) +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "navy"), color = "navy") +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "lightblue"), color = "lightblue") +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "green"), color = "green") +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "orange"), color = "orange") +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "red"), color = "red") +
    geom_point(data = function(x) subset(x, color_to_rs4981455_eqtl == "purple"), color = "purple") +
    theme(legend.position = "none") +
    geom_hline(yintercept = -log10(2.804040044e-05), linetype = "dashed") +
    labs(title = "Secondary miRNA-eQTL Signal at miR-4707-3p")

ggsave("secondary_signal.pdf", height = 4, width = 8)












