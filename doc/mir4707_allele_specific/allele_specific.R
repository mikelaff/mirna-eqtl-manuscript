
library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# s

# INPUT FILES ##########################################################################################################
# bowtie counts for miR-4707-3p
counts.tsv <- here("results/counts/small_rna_seq/bowtie_mir4707/20210405_mir4707_counts.tsv")

# genotypes directory
dir.genotypes <- here("results/emmax/association_results/20200120_mirQTLor/sample_genotypes/")

# mirna expression rse
mirna.rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")
# vsd for the mirQTL
vsd.mirna.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################
eQTL <- "hsa-mir-4707_hsa-miR-4707-3p+chr14:22953244:A:G"
emiR <- "hsa-mir-4707_hsa-miR-4707-3p"
eSNP <- "chr14:22953244:A:G"
chr <- "chr14"

rsid <- "rs4981455"
ensg <- "ENSG00000092036"

eSNP.seed <- "chr14:22956973:C:A"
rsid.seed <- "rs2273626"

# Import ###############################################################################################################

# import mirQTLor vsd
vsd <- read_rds(vsd.mirna.rds)
# samples used in eqtl analysis
df.samples <- as_tibble(colData(vsd))

# get genotypes for these samples at 4707 index
df.genotypes <- read_table(paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", eSNP, ".raw"))
df.genotypes %<>%
    dplyr::select(donor_id = 1,
                  #DNAID = 2,
                  index_genotypes_number_of_G = 7)
# get genotypes for these samples at 4707 seed
df.genotypes.seed <- read_table(paste0(dir.genotypes, chr, ".hardcall.prefiltered.mirQTLor.", eSNP.seed, ".raw"))
df.genotypes.seed %<>%
    dplyr::select(donor_id = 1,
                  #DNAID = 2,
                  seed_genotypes_number_of_A = 7)

df.samples %<>%
    mutate(donor_id = paste0("D", donor_id))

df.samples %<>%
    left_join(df.genotypes, by = "donor_id")
df.samples %<>%
    left_join(df.genotypes.seed, by = "donor_id")


df.samples %<>%
    mutate(index_genotypes_number_of_G = factor(index_genotypes_number_of_G, levels = c(0,1,2), labels = c("AA", "AG", "GG"), ordered = TRUE))
df.samples %<>%
    mutate(seed_genotypes_number_of_A = factor(seed_genotypes_number_of_A, levels = c(0,1,2), labels = c("CC", "CA", "AA"), ordered = TRUE))

df.samples %<>%
    dplyr::rename(genotype_index_rs4981455 = index_genotypes_number_of_G,
                  genotype_seed_rs2273626 = seed_genotypes_number_of_A)


# import counts data
rse <- read_rds(mirna.rse.rds)
# subset to these samples
rse <- rse[,colnames(vsd)]
dds <- DESeqDataSet(rse, design = ~1)

#counts(dds)
#fpm(dds)





df.counts <- as.data.frame(read_tsv(counts.tsv))
rownames(df.counts) <- df.counts$name
df.counts <- df.counts[,7:246]
df.counts <- as.data.frame(t(df.counts))
colnames(df.counts) <- c("miR4707.SNPC", "miR4707.SNPA", "miR4707.ref")
df.counts$rnaid <- rownames(df.counts)

df.counts <- as_tibble(df.counts)
df.counts %<>%
    filter(rnaid %in% df.samples$rnaid)

df.mirge <- tibble(rnaid = df.counts$rnaid,
                   miR4707.miRge.counts = counts(dds)["hsa-miR-4707-3p", df.counts$rnaid],
                   miR4707.miRge.cpm = fpm(dds)["hsa-miR-4707-3p", df.counts$rnaid])

df.counts %<>%
    left_join(df.mirge, by = "rnaid")

df.counts %<>%
    left_join(select(df.samples, rnaid, genotype_index_rs4981455, genotype_seed_rs2273626), by = "rnaid")

# paired t-test
df.counts %>%
    filter(genotype_seed_rs2273626 == "CA") -> df.hets

t.test(df.hets$miR4707.SNPC, df.hets$miR4707.SNPA, paired = TRUE, alternative = "two.sided")

# Plot ##############################################################################

df.counts %>%
    ggplot(aes(x = genotype, y = miR4707.miRge.cpm)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/mir-4707-3p_miRge_cpm_boxplot.pdf", height = 2, width = 2)


df.counts %>%
    ggplot(aes(x = genotype, y = miR4707.miRge.counts)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.1, size = 0.5) +
    plotTheme("figure")

ggsave("~/Desktop/mir-4707-3p_miRge_counts_boxplot.pdf", height = 2, width = 2)


df.counts %>%
    ggplot(aes(x = genotype, y = miR4707.ref)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(height = 0, width = 0.1, size = 0.5) +
    plotTheme("figure")

df.counts %>%
    select(1,2,4,7) %>%
    melt(measure.vars = c("miR4707.SNPC", "miR4707.SNPA")) %>%
    ggplot(aes(x = genotype, y = value, color = variable)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_jitter(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.5) +
    plotTheme("figure")


df.counts %>%
    filter(genotype == "AG") %>%
    select(1,2,4,7) %>%
    melt(measure.vars = c("miR4707.SNPC", "miR4707.SNPA")) %>%
    ggplot(aes(x = variable, y = value)) +
    geom_line(aes(group = rnaid)) +
    plotTheme("figure")

ggsave("~/Desktop/mir-4707-3p_miRge_ase_boxplot.pdf", height = 2, width = 2)







