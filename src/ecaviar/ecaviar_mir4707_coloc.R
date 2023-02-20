# get emiRs and eSNPs to get eQTLs.

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(reshape2)
library(GenomicRanges)
library(mikelaffr)

library(pheatmap)

# OUTPUT FILES ########################################################################################################
# output dir for ecaviar files
dir.output <- here("results/ecaviar/")
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

# SNPs in the 1000G EUR
gwas.bim.file <- here("results/ecaviar/gwas.ld.bim")

# SNPs in the eQTL
eqtl.bim.file <- here("results/ecaviar/eqtl.ld.bim")

# recode AD .raw file for eQTL
eqtl.ad.raw <- here("results/ecaviar/eqtl.ld.raw")

# recode AD .raw file for gwas
gwas.ad.raw <- here("results/ecaviar/gwas.ld.raw")


# GLOBALS ##############################################################################################################
# miR-4707-3p eQTL
chr <- "chr14"
eSNP <- "chr14:22953244:A:G"
emiR <- "hsa-mir-4707_hsa-miR-4707-3p"

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

df.sumstats.edu %<>%
    mutate(PVE = (2 * (Beta^2) * MAF * (1 - MAF)) /
               ((2 * (Beta^2) * MAF * (1 - MAF)) + ((SE^2) * 2 * 1.1e6 * MAF * (1 - MAF)))
    )

# mRNA-eQTL
df.mrna <- readRDS(mrna.results.chr14.rds)
df.mrna %<>%
    filter(ENSG == haus4.ensg)

df.gwas.bim <- read_table(gwas.bim.file, col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))

df.eqtl.bim <- read_table(eqtl.bim.file, col_names = c("chr", "rsid", "cm", "bp", "a1", "a2"))

# eQTL SE and Z-score ##################################################################################################

# number of samples
num.samples <- 212

# degrees of freedom, add 1 for secondary, 2 for tertiary, etc...
# design equation: expression ~ PC1.genotype + PC2.genotype + PC3.genotype + PC4.genotype + PC5.genotype + PC6.genotype +
# PC7.genotype + PC8.genotype + PC9.genotype + PC10.genotype + PC1.expression + PC2.expression + PC3.expression + PC4.expression +
# PC5.expression + PC6.expression + PC7.expression + PC8.expression + PC9.expression + PC10.expression + PoolPool2 + PoolPool3 +
# PoolPool4 + PoolPool5 + PoolPool6 + PoolPool7 + PoolPool8 + PurificationMethodmiRNeasy + PurificationMethodmiRNeasy_mini +
# SexM + RIN + GestationWeek + primaryDosage + secondaryDosage + tertiaryDosage
# add 1 for kinship matrix

primaryDF <- num.samples - (32 + 1) - 1

# compute MAF
df.eqtl.results %<>%
    mutate(MAF = ALT_CTS / OBS_CT) %>%
    mutate(MAF = ifelse(MAF > 0.5, 1 - MAF, MAF))


# compute standard error
df.eqtl.results %<>%
    mutate(SE_beta = getSE(BETA, P, num.samples, primaryDF))


# compute z-score
df.eqtl.results %<>%
    mutate(Z_score = BETA / SE_beta)

# Format mRNA-eQTL Data ################
df.mrna %<>%
    mutate(SE_beta = getSE(BETA, P, num.samples, primaryDF))

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
           BETA,
           P,
           SE_beta,
           Z_score)

df.sumstats.edu %<>%
    select(RSID,
           CHR = seqnames,
           BP = start,
           A1 = A1.effect,
           A2,
           EFFECT.ALLELE = A1.effect,
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



# df.combined %>%
#     filter((P.eQTL < 0.05) | (P.GWAS < 0.05) | (P.mRNA < 0.05)) -> df.combined.pp

df.combined %<>%
    filter(RSID %in% df.gwas.bim$rsid)

# df.combined %<>%
#     filter((P.eQTL < 1.434e-6) | (P.GWAS < 5e-8) | (P.mRNA < 8.17e-4))



df.eqtl.raw <- read_table(eqtl.ad.raw)

df.eqtl.raw %>%
    select(!ends_with("_HET"),
           -FID,
           -IID,
           -PAT,
           -MAT,
           -SEX,
           -PHENOTYPE) -> df.eqtl.raw.variants

df.cols.eqtl <- tibble(var.name = colnames(df.eqtl.raw.variants))
df.cols.eqtl %<>%
    mutate(a1.allele = as.character(lapply(strsplit(var.name, "_"), `[`, 2)),
           SNP = as.character(lapply(strsplit(var.name, "_"), `[`, 1)))

df.cols.eqtl %<>%
    left_join(select(df.combined, SNP, RSID, REF, ALT, EFFECT.ALLELE.eQTL, A1.GWAS = A1, A2.GWAS = A2, EFFECT.ALLELE.GWAS), by = "SNP")

all(df.cols.eqtl$a1.allele == df.eqtl.bim$a1)
all(df.cols.eqtl$SNP == df.combined$SNP)

cor(df.eqtl.raw.variants[,1], df.eqtl.raw.variants[,5])


df.gwas.raw <- read_table(gwas.ad.raw)

df.gwas.raw %>%
    select(!ends_with("_HET"),
           -FID,
           -IID,
           -PAT,
           -MAT,
           -SEX,
           -PHENOTYPE) -> df.gwas.raw.variants

df.cols.gwas <- tibble(var.name = colnames(df.gwas.raw.variants))
df.cols.gwas %<>%
    mutate(a1.allele = as.character(lapply(strsplit(var.name, "_"), `[`, 2)),
           RSID = as.character(lapply(strsplit(var.name, "_"), `[`, 1)))

df.cols.gwas %<>%
    left_join(select(df.combined, SNP, RSID, REF, ALT, EFFECT.ALLELE.eQTL, A1.GWAS = A1, A2.GWAS = A2, EFFECT.ALLELE.GWAS), by = "RSID")

all(df.cols.gwas$a1.allele == df.gwas.bim$a1)
all(df.cols.gwas$RSID == df.combined$RSID)

cor(df.gwas.raw.variants[,1], df.gwas.raw.variants[,5])

# output

# eQTL snps for LD matrix
write_lines(df.combined$SNP, paste0(dir.output, "eqtl.snps.txt"))
# mRNA eQTL snps for LD matrix
write_lines(df.combined$SNP, paste0(dir.output, "mrna.snps.txt"))
# GWAS snps for LD matrix
write_lines(df.combined$RSID, paste0(dir.output, "gwas.snps.txt"))

# z scores
write_delim(select(df.combined, RSID, Z_score.eQTL), col_names = FALSE, file = paste0(dir.output, "eqtl.zscore"))
write_delim(select(df.combined, RSID, Z_score.mRNA.fixed), col_names = FALSE, file = paste0(dir.output, "mrna.zscore"))
write_delim(select(df.combined, RSID, Z_score.GWAS.fixed), col_names = FALSE, file = paste0(dir.output, "gwas.zscore"))


# correct LD matrix
all(df.gwas.bim$bp == df.eqtl.bim$bp)

vec.eqtl.correction <- ifelse(df.combined$EFFECT.ALLELE.eQTL == df.eqtl.bim$a1, 1, -1)
vec.gwas.correction <- ifelse(df.combined$EFFECT.ALLELE.eQTL == df.gwas.bim$a1, 1, -1)

mat.eqtl.ld <- as.matrix(read_table(paste0(dir.output, "eqtl.ld.ld"), col_names = FALSE))
mat.gwas.ld <- as.matrix(read_table(paste0(dir.output, "gwas.ld.ld"), col_names = FALSE))

mat.eqtl.ld.corrected <- t(t(mat.eqtl.ld) * vec.eqtl.correction) * vec.eqtl.correction
mat.gwas.ld.corrected <- t(t(mat.gwas.ld) * vec.gwas.correction) * vec.gwas.correction

write.table(mat.eqtl.ld.corrected, file = paste0(dir.output, "eqtl.ld.ld.corrected"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(mat.gwas.ld.corrected, file = paste0(dir.output, "gwas.ld.ld.corrected"), row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")


pheatmap(mat.eqtl.ld, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(mat.gwas.ld, cluster_rows = FALSE, cluster_cols = FALSE)

pheatmap(mat.eqtl.ld.corrected, cluster_rows = FALSE, cluster_cols = FALSE)
pheatmap(mat.gwas.ld.corrected, cluster_rows = FALSE, cluster_cols = FALSE)

# df.combined %<>%
#     mutate(Z_score.GWAS.fixed.flipped = -Z_score.GWAS.fixed)
#
# write_delim(select(df.combined, RSID, Z_score.GWAS.fixed.flipped), col_names = FALSE, file = paste0(dir.output, "gwas.zscore.flipped"))


#

# #
#
# command <- paste0("module load plink/1.90b6.21 && plink ",
#                   "--bfile ", ,
#                   "--")
#
# command <- paste0('module load plink/1.90b6.21 && plink',
#                   128   ' --bfile ',gen.fn,
#                   129   ' --chr ', sub('chr','',chr),
#                   130   ' --from-bp ', start.pos ,' --to-bp ',end.pos,
#                   131   ' --r square ',
#                   132   ' --out ', ld.fn,'.eQTL')
# 133  print(command)
# 134  system(command,ignore.stdout=T)
# 135 }


# Results ##################


list.files(paste0(dir.output))

df.ecaviar.mirTOgwas.1col <- read_table(paste0(dir.output, "mir4707_to_gwas_1causal_col"))
df.ecaviar.mirTOgwas.1col %<>%
    select(RSID = SNP_ID,
           Prob_in_pCausalSet_mirTOgwas.1causal = Prob_in_pCausalSet,
           CLPP_mirTOgwas.1causal = CLPP)

df.ecaviar.mirTOgwas.2col <- read_table(paste0(dir.output, "mir4707_to_gwas_2causal_col"))
df.ecaviar.mirTOgwas.2col %<>%
    select(RSID = SNP_ID,
           Prob_in_pCausalSet_mirTOgwas.2causal = Prob_in_pCausalSet,
           CLPP_mirTOgwas.2causal = CLPP)

df.ecaviar.mirTOmrna.1col <- read_table(paste0(dir.output, "mir4707_to_mrna_1causal_col"))
df.ecaviar.mirTOmrna.1col %<>%
    select(RSID = SNP_ID,
           Prob_in_pCausalSet_mirTOmrna.1causal = Prob_in_pCausalSet,
           CLPP_mirTOmrna.1causal = CLPP)

df.ecaviar.mirTOmrna.2col <- read_table(paste0(dir.output, "mir4707_to_mrna_2causal_col"))
df.ecaviar.mirTOmrna.2col %<>%
    select(RSID = SNP_ID,
           Prob_in_pCausalSet_mirTOmrna.2causal = Prob_in_pCausalSet,
           CLPP_mirTOmrna.2causal = CLPP)

# df.ecaviar.mirTOgwas.1col.flip <- read_table(paste0(dir.output, "mir4707_to_gwas_1causal_flip_col"))
# df.ecaviar.mirTOgwas.1col.flip %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.1causal.flip = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.1causal.flip = CLPP)
#
# df.ecaviar.mirTOgwas.2col.flip <- read_table(paste0(dir.output, "mir4707_to_gwas_2causal_flip_col"))
# df.ecaviar.mirTOgwas.2col.flip %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.2causal.flip = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.2causal.flip = CLPP)
#
# df.ecaviar.mirTOgwas.1col.eur <- read_table(paste0(dir.output, "mir4707_to_gwas_1causal_EUR_col"))
# df.ecaviar.mirTOgwas.1col.eur %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.1causal.EUR = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.1causal.EUR = CLPP)
#
# df.ecaviar.mirTOgwas.2col.eur <- read_table(paste0(dir.output, "mir4707_to_gwas_2causal_EUR_col"))
# df.ecaviar.mirTOgwas.2col.eur %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.2causal.EUR = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.2causal.EUR = CLPP)
#
# df.ecaviar.mirTOgwas.1col.zflip <- read_table(paste0(dir.output, "mir4707_to_gwas_1causal_Zflip_col"))
# df.ecaviar.mirTOgwas.1col.zflip %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.1causal.zflip = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.1causal.zflip = CLPP)
#
# df.ecaviar.mirTOgwas.2col.zflip <- read_table(paste0(dir.output, "mir4707_to_gwas_2causal_Zflip_col"))
# df.ecaviar.mirTOgwas.2col.zflip %<>%
#     select(RSID = SNP_ID,
#            Prob_in_pCausalSet_mirTOgwas.2causal.zflip = Prob_in_pCausalSet,
#            CLPP_mirTOgwas.2causal.zflip = CLPP)

df.ecaviar.mirTOgwas.1col %>%
    left_join(df.ecaviar.mirTOgwas.2col, by = "RSID") -> df.comb.TOgwas

df.ecaviar.mirTOmrna.1col %>%
    left_join(df.ecaviar.mirTOmrna.2col, by = "RSID") -> df.comb.TOmrna

# df.ecaviar.mirTOgwas.1col.flip %>%
#     left_join(df.ecaviar.mirTOgwas.2col.flip, by = "RSID") -> df.comb.TOgwas.flip
#
# df.ecaviar.mirTOgwas.1col.eur %>%
#     left_join(df.ecaviar.mirTOgwas.2col.eur, by = "RSID") -> df.comb.TOgwas.eur
#
# df.ecaviar.mirTOgwas.1col.zflip %>%
#     left_join(df.ecaviar.mirTOgwas.2col.zflip, by = "RSID") -> df.comb.TOgwas.zflip

df.comb.TOgwas %>%
    left_join(df.comb.TOmrna, by = "RSID") -> df.combined.results

# df.combined.results %<>%
#     left_join(df.comb.TOgwas.flip, by = "RSID")
#
# df.combined.results %<>%
#     left_join(df.comb.TOgwas.eur, by = "RSID")
#
# df.combined.results %<>%
#     left_join(df.comb.TOgwas.zflip, by = "RSID")

df.combined.results %<>%
    left_join(df.combined, by = "RSID")


#melt(df.comb.TOgwas, id.vars = c("RSID"))

df.combined.results %<>%
    mutate(SIG.IN.THREE = (P.eQTL < 1.434e-6) & (P.GWAS < 5e-8) & (P.mRNA < 8.17e-4))

melt(df.combined.results, measure.vars = c("CLPP_mirTOgwas.1causal", "CLPP_mirTOgwas.2causal", "CLPP_mirTOmrna.1causal", "CLPP_mirTOmrna.2causal")) -> df.combined.results.melt

df.combined.results %>%
    filter(BP > 22940000, BP < 22980000) -> df.tmp

melt(df.tmp, measure.vars = c("CLPP_mirTOgwas.1causal", "CLPP_mirTOgwas.2causal", "CLPP_mirTOmrna.1causal", "CLPP_mirTOmrna.2causal")) -> df.tmp.clpp
melt(df.tmp, measure.vars = c("P.eQTL", "P.GWAS", "P.mRNA")) -> df.tmp.pval
melt(df.tmp, measure.vars = c("BETA.eQTL", "BETA.GWAS.fixed", "BETA.mRNA.fixed")) -> df.tmp.beta
melt(df.tmp, measure.vars = c("Z_score.eQTL", "Z_score.GWAS.fixed", "Z_score.mRNA.fixed")) -> df.tmp.z_score
melt(df.tmp, measure.vars = c("Prob_in_pCausalSet_mirTOgwas.1causal", "Prob_in_pCausalSet_mirTOgwas.2causal", "Prob_in_pCausalSet_mirTOmrna.1causal", "Prob_in_pCausalSet_mirTOmrna.2causal")) -> df.tmp.prob

df.tmp.clpp <- as_tibble(df.tmp.clpp)
df.tmp.pval <- as_tibble(df.tmp.pval)
df.tmp.beta <- as_tibble(df.tmp.beta)
df.tmp.prob <- as_tibble(df.tmp.prob)
df.tmp.z_score <- as_tibble(df.tmp.z_score)

# df.combined.results %>%
#     top_n(n = 40, wt = CLPP_mirTOgwas.1causal) %>%
#     ggplot(aes(x = BP, y = CLPP_mirTOgwas.1causal)) +
#     geom_point()
#
# df.combined.results %>%
#     top_n(n = 40, wt = CLPP_mirTOgwas.1causal) %>%
#     ggplot(aes(x = BP, y = CLPP_mirTOgwas.2causal)) +
#     geom_point()
#
# df.combined.results %>%
#     top_n(n = 40, wt = CLPP_mirTOgwas.1causal) %>%
#     ggplot(aes(x = BP, y = CLPP_mirTOmrna.1causal)) +
#     geom_point()
#
# df.combined.results %>%
#     top_n(n = 40, wt = CLPP_mirTOmrna.2causal) %>%
#     ggplot(aes(x = BP, y = CLPP_mirTOmrna.2causal)) +
#     geom_point()

pdf("~/Desktop/ecaviar.pdf", width = 16, height = 8)

df.combined.results.melt %>%
    ggplot(aes(x = BP, y = value, color = variable)) +
    geom_line() +
    geom_hline(yintercept = 0.1, linetype = "dashed") +
    labs(y = "CLPP")

df.tmp.clpp %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable)) +
    geom_point() +
    geom_hline(yintercept = 0.1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(y = "CLPP") +
    geom_point(data = subset(df.tmp.clpp, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.clpp, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1)

ggsave("~/Desktop/ecaviar_sigSNPs.pdf", height = 8, width = 12)

df.tmp.prob %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable)) +
    geom_point() +
    #geom_hline(yintercept = 0.1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(y = "Prob in pCausalSet") +
    geom_point(data = subset(df.tmp.prob, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.prob, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1)

df.tmp.prob %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable)) +
    geom_point() +
    #geom_hline(yintercept = 0.1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Prob in pCausalSet") +
    geom_point(data = subset(df.tmp.prob, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.prob, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1)

df.tmp.pval %>%
    ggplot(aes(x = reorder(RSID, BP), y = -log10(value), color = variable)) +
    geom_point() +
    geom_hline(yintercept = -log10(1.434e-6), linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~variable) +
    labs(y = "-log10(p-val)") +
    geom_point(data = subset(df.tmp.pval, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = -log10(value)), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.pval, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = -log10(value)), size = 4, shape = 1)

df.tmp.pval %>%
    ggplot(aes(x = reorder(RSID, BP), y = -log10(value), color = variable)) +
    geom_point() +
    geom_hline(yintercept = -log10(1.434e-6), linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "-log10(p-val)") +
    geom_point(data = subset(df.tmp.pval, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = -log10(value)), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.pval, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = -log10(value)), size = 4, shape = 1)

df.tmp.beta %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable, shape = SIG.IN.THREE)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_wrap(~variable, scales = "free_y") +
    labs(y = "BETA") +
    geom_point(data = subset(df.tmp.beta, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.beta, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1) +
    scale_shape_manual(values = c(1,16))

df.tmp.beta %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable, shape = SIG.IN.THREE)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "BETA") +
    geom_point(data = subset(df.tmp.beta, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.beta, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1) +
    scale_shape_manual(values = c(1,16))

df.tmp.z_score %>%
    ggplot(aes(x = reorder(RSID, BP), y = value, color = variable, shape = SIG.IN.THREE)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(y = "Z_score") +
    geom_point(data = subset(df.tmp.z_score, RSID == "rs4981455"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 5) +
    geom_point(data = subset(df.tmp.z_score, RSID == "rs2273626"), aes(x = reorder(RSID, BP), y = value), size = 4, shape = 1) +
    scale_shape_manual(values = c(1,16))

dev.off()

df.combined.pp %>%
    filter(P.eQTL < 0.05 & P.GWAS < 0.05) %>%
    ggplot(aes(x = -log10(P.eQTL), y = -log10(P.GWAS))) +
    geom_point() +
    labs(caption = "P-values < 0.05 in eQTL and GWAS")

ggsave("~/Desktop/pp_plot.pdf", width = 10, height = 10)



# Scratch

mat <- matrix(c(1.000,  0.840,  0.050, 0.050,  0.001, -0.010,
                0.840,  1.000,  0.040, 0.040, -0.010, -0.007,
                0.050,  0.040,  1.000, 0.950,  0.060, -0.001,
                0.050,  0.040,  0.950, 1.000,  0.06 ,  0.003,
                0.001, -0.010,  0.060, 0.065,  1.000,  0.018,
                -0.01, -0.007, -0.001, 0.003,  0.018,  1.000), nrow = 6, ncol = 6)

mat

vec <- c(-1, 1, 1, 1, 1, 1)

t(t(mat) * vec) * vec

vec <- c(-1, -1, 1, 1, 1, 1)

t(t(mat) * vec) * vec





