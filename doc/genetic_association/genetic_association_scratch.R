# look at genetic association, scratch

library(vcfR)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(here)

# INPUT FILES #####################################################################################
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")

# Build DDS #######################################################################################
# read in summarized experiment
se <- readRDS(mirbase_se_rds)

# label outliers by PCA (after batch effect removal of pool and rna purif. method, see below)
outlier_rnaid <- c("RNAID1710", "RNAID1744", "RNAID1698", "RNAID1814", "RNAID1678")
colData(se)$outlier <- ifelse(colData(se)$rnaid %in% outlier_rnaid, "yes", "no")

# build DESeq data set
dds <- DESeqDataSet(se, design = ~1)
# remove rows with only zero counts
dds <- dds[rowSums(counts(dds)) > 1, ]
# estimate size factors for normalization
dds <- estimateSizeFactors(dds)
# norm. transformed counts (DESeq size factors, log2 transformed)
#ntc <- assay(normTransform(dds))
# variance stabilizing tranformed data
vsd <- varianceStabilizingTransformation(dds)
rm(se)

#####################################################

# sample data
samples.df <- as.data.frame(colData(dds))
samples.df$donor_id <- paste("D", samples.df$donor_id, sep="")

# top 10,000 lines from the latest genotype data from Dan
vcf <- read.vcfR("~/Desktop/allmgergedGenoMafHweQCnewVCFhg38_top10000lines.vcf")

vcf.tidy <- vcfR2tidy(vcf)

# all the variants
variants <- vcf.tidy$fix
# genotypes for all the samples at all the variants
gts <- vcf.tidy$gt
rm(vcf.tidy, vcf)

# sample data from the genotype file
gt.samples.df <- data.frame(individual = as.character(unique(gts$Indiv)), stringsAsFactors = FALSE)
gt.samples.df$donor_id <- sapply(strsplit(gt.samples.df$individual, "_"), `[`, 1)

dups <- gt.samples.df[duplicated(gt.samples.df$donor_id) | duplicated(gt.samples.df$donor_id, fromLast=TRUE),]

# non duplicated samples from genotype file
gt.samples.df <- gt.samples.df[!duplicated(gt.samples.df$donor_id),]

# merge genotype samples with mirna samples
samples.df <- left_join(samples.df, gt.samples.df, by = "donor_id")

# plot genotype information

chrom <- 1
pos <- 598867

i <- 7
chrom <- variants$ChromKey[i]
pos <- variants$POS[i]

gts %>%
  filter(ChromKey == chrom, POS == pos) %>%
  filter(Indiv %in% samples.df$individual) -> tmp

samples.df %>%
  left_join(data.frame(rnaid = colnames(vsd), exprs = assay(vsd)[235,], stringsAsFactors = FALSE), by = "rnaid") %>%
  left_join(tmp, by = c("individual" = "Indiv")) %>%
  filter(!is.na(gt_GT_alleles)) %>%
  ggplot(aes(x=gt_GT_alleles, y=exprs)) +
  geom_boxplot()

samples.df %>%
  left_join(data.frame(rnaid = colnames(vsd), exprs = assay(vsd)[235,], stringsAsFactors = FALSE), by = "rnaid") %>%
  left_join(tmp, by = c("individual" = "Indiv")) %>%
  filter(!is.na(gt_GT_alleles)) -> tmp2

tmp2$gt_GT_alleles <- factor(tmp2$gt_GT_alleles)
tmp2$effect_alleles <- ifelse(tmp2$gt_GT_alleles == "G/G", 0, ifelse(tmp2$gt_GT_alleles == "G/A", 1, 2))

model.matrix(~gt_GT_alleles, data = tmp2)

fit <- lm(mds_c2 ~ gt_GT_alleles, data = tmp2)

summary(fit)

tmp2 %>%
  ggplot(aes(x=gt_GT_alleles, y=mds_c2)) +
  geom_boxplot()

fit2 <- lm(mds_c2 ~ effect_alleles, data = tmp2)
summary(fit2)

tmp2 %>%
  ggplot(aes(x=effect_alleles, y=mds_c2)) +
  geom_point() +
  geom_smooth(method = "lm")

tmp2 %>%
  ggplot(aes(x=effect_alleles, y=mds_c2)) +
  geom_point() +
  geom_abline(slope=0.007747, intercept = -0.0202)


