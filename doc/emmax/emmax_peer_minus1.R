
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(mikelaffr)

# OUTPUT ##############
dir.pdf <- here("doc/emmax/pdfs/")

# INPUT ##############
files.rds <- list.files(here("results/emmax/association_results/20191121_test_covariates_peer_minus1/"), pattern = "*.rds")

# Load Results Summary Tables ##################
# combine all summary tables into one
df.eqtls <- tibble()

for (df.file in files.rds) {
    print(df.file)

    df.tmp <- as_tibble(readRDS(paste0(here("results/emmax/association_results/20191121_test_covariates_peer_minus1/"), df.file)))

    df.eqtls <- bind_rows(df.eqtls, df.tmp)

}

rm(df.tmp)

df.eqtls

df.eqtls$eqtl <- paste(df.eqtls$emir, df.eqtls$esnp, sep = "+")

# df.emirs %>%
#     group_by(covMat) %>%
#     summarise(emirCount = n()) -> df.cov

df.eqtls$genotypePCA <- sapply(strsplit(df.eqtls$covMat, "\\."), `[`, 2)
df.eqtls$genotypeNumComps <- as.numeric(sapply(strsplit(df.eqtls$genotypePCA, "_"), `[`, 2))

df.eqtls$peer <- sapply(strsplit(df.eqtls$covMat, "\\."), `[`, 3)
df.eqtls$peerFactors <- as.numeric(sapply(strsplit(df.eqtls$peer, "_"), `[`, 2))

df.eqtls %>%
    group_by(covMat) %>%
    summarise(emirCount = n_distinct(emir), esnpCount = n_distinct(esnp), eqtlCount = n_distinct(eqtl)) -> df.cov

df.cov$genotypePCA <- sapply(strsplit(df.cov$covMat, "\\."), `[`, 2)
df.cov$genotypeNumComps <- as.numeric(sapply(strsplit(df.cov$genotypePCA, "_"), `[`, 2))

df.cov$peer <- sapply(strsplit(df.cov$covMat, "\\."), `[`, 3)
df.cov$peerFactors <- as.numeric(sapply(strsplit(df.cov$peer, "_"), `[`, 2))

df.cov$batchVars <- grepl("pool", df.cov$covMat)

df.cov$peerFactors <- ifelse(df.cov$peerFactors != 0, df.cov$peerFactors - 1, 0)

df.cov %>%
    ggplot(aes(x=peerFactors, y=emirCount, color=batchVars)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of PEER Factors",
         y = "Number of emiRs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (PEER Factors) + (Seq Pool + Purif Method + RIN + Sex + Gest Week)",
         title = "Inclusions of PEER Factors and Known Tech/Bio Covariates (Excluding PEER1)")

ggsave(paste0(dir.pdf, "emirs_by_peerFactors_excludingPEER1.pdf"), width = 9, height = 5)

df.cov %>%
    ggplot(aes(x=peerFactors, y=esnpCount, color=batchVars)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of PEER Factors",
         y = "Number of eSNPs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (PEER Factors) + (Seq Pool + Purif Method + RIN + Sex + Gest Week)",
         title = "Inclusions of PEER Factors and Known Tech/Bio Covariates (Excluding PEER1)")

df.cov %>%
    ggplot(aes(x=peerFactors, y=eqtlCount, color=batchVars)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of PEER Factors",
         y = "Number of eQTLs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (PEER Factors) + (Seq Pool + Purif Method + RIN + Sex + Gest Week)",
         title = "Inclusions of PEER Factors and Known Tech/Bio Covariates (Excluding PEER1)")

ggsave(paste0(dir.pdf, "eqtls_by_peerFactors_excludingPEER1.pdf"), width = 9, height = 5)


# Comparison to Expression PCs ####################################

df.cov.peer <- df.cov
df.eqtls.peer <- df.eqtls

# INPUT ##############
files.rds <- list.files(here("results/emmax/association_results/20191118_test_covariates_expression_pcs/"), pattern = "*.rds")

# Load Results Summary Tables ##################
# combine all summary tables into one
df.eqtls <- tibble()

for (df.file in files.rds) {
    print(df.file)

    df.tmp <- as_tibble(readRDS(paste0(here("results/emmax/association_results/20191118_test_covariates_expression_pcs/"), df.file)))

    df.eqtls <- bind_rows(df.eqtls, df.tmp)

}

rm(df.tmp)

df.eqtls

df.eqtls$eqtl <- paste(df.eqtls$emir, df.eqtls$esnp, sep = "+")

# df.emirs %>%
#     group_by(covMat) %>%
#     summarise(emirCount = n()) -> df.cov

df.eqtls$genotypePCA <- sapply(strsplit(df.eqtls$covMat, "\\."), `[`, 2)
df.eqtls$genotypeNumComps <- as.numeric(sapply(strsplit(df.eqtls$genotypePCA, "_"), `[`, 2))

df.eqtls$expressionPCA <- sapply(strsplit(df.eqtls$covMat, "\\."), `[`, 3)
df.eqtls$expressionNumComps <- as.numeric(sapply(strsplit(df.eqtls$expressionPCA, "_"), `[`, 2))

df.eqtls %>%
    group_by(covMat) %>%
    summarise(emirCount = n_distinct(emir), esnpCount = n_distinct(esnp), eqtlCount = n_distinct(eqtl)) -> df.cov

df.cov$genotypePCA <- sapply(strsplit(df.cov$covMat, "\\."), `[`, 2)
df.cov$genotypeNumComps <- as.numeric(sapply(strsplit(df.cov$genotypePCA, "_"), `[`, 2))

df.cov$expressionPCA <- sapply(strsplit(df.cov$covMat, "\\."), `[`, 3)
df.cov$expressionNumComps <- as.numeric(sapply(strsplit(df.cov$expressionPCA, "_"), `[`, 2))

df.cov$batchVars <- grepl("pool", df.cov$covMat)

df.cov.pca <- df.cov
df.eqtls.pca <- df.eqtls

rm(df.cov, df.eqtls)

# combine and plot

df.cov.peer$expressionCorrection <- "PEER"
df.cov.pca$expressionCorrection <- "PCA"

df.cov <- bind_rows(select(df.cov.pca,
                           covMat,
                           emirCount,
                           esnpCount,
                           eqtlCount,
                           genotypePCA,
                           genotypeNumComps,
                           correction = expressionPCA,
                           correctionComps = expressionNumComps,
                           batchVars,
                           expressionCorrection),
                    select(df.cov.peer,
                           covMat,
                           emirCount,
                           esnpCount,
                           eqtlCount,
                           genotypePCA,
                           genotypeNumComps,
                           correction = peer,
                           correctionComps = peerFactors,
                           batchVars,
                           expressionCorrection))

df.cov %>%
    filter(batchVars) %>%
    ggplot(aes(x=correctionComps, y=eqtlCount, color=expressionCorrection)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of PEER Factors/Expression PCs",
         y = "Number of eQTLs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (PEER Factors/Expression PCs) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Comparison of Expression PCs vs PEER Factors (Excluding PEER1, Including Batch Vars.)")

ggsave(paste0(dir.pdf, "eqtls_by_peer_or_expressionPCs_excludingPEER1.pdf"), width = 9, height = 5)




