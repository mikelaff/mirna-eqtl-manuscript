
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(mikelaffr)

# OUTPUT ##############
dir.pdf <- here("doc/emmax/pdfs/")

# INPUT ##############
files.rds <- list.files(here("results/emmax/association_results/20191122_mirQTLor_test_covariates_expression_pcs/"), pattern = "*.rds")

# Load Results Summary Tables ##################
# combine all summary tables into one
df.eqtls <- tibble()

for (df.file in files.rds) {
    print(df.file)

    df.tmp <- as_tibble(readRDS(paste0(here("results/emmax/association_results/20191122_mirQTLor_test_covariates_expression_pcs/"), df.file)))

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
df.cov$expressionPCA <- sapply(strsplit(df.cov$expressionPCA, "_"), `[`, 1)
df.cov$expressionPCA <- ifelse(grepl("full", df.cov$expressionPCA), "Full", "miRBase Only")

df.cov$batchVars <- grepl("pool", df.cov$covMat)

df.cov %>%
    ggplot(aes(x=expressionNumComps, y=emirCount, color=expressionPCA)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Expression PCs",
         y = "Number of emiRs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (Expression PCs) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusion of Expression PCs: miRBase Expression Only or All Known and Novel Expression")

df.cov %>%
    ggplot(aes(x=expressionNumComps, y=esnpCount, color=expressionPCA)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Expression PCs",
         y = "Number of eSNPs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (Expression PCs) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusion of Expression PCs: miRBase Expression Only or All Known and Novel Expression")

df.cov %>%
    ggplot(aes(x=expressionNumComps, y=eqtlCount, color=expressionPCA)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Expression PCs",
         y = "Number of eQTLs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (Expression PCs) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusion of Expression PCs: miRBase Expression Only or All Known and Novel Expression")

#ggsave(paste0(dir.pdf, "eqtls_by_mirQTLor_expressionPCs.pdf"), width = 9, height = 5)

# Comparison to mirQTL ##########

df.cov.mirQTLor <- df.cov
df.eqtls.mirQTLor <- df.eqtls

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
df.cov$expressionPCA <- "Full"

df.cov$batchVars <- grepl("pool", df.cov$covMat)

df.cov.mirQTL <- df.cov
df.eqtls.mirQTL <- df.eqtls

rm(df.cov, df.eqtls)

# combine and plot

df.cov.mirQTL$samples <- "mirQTL"
df.cov.mirQTLor$samples <- "mirQTLor"

df.cov <- bind_rows(select(filter(df.cov.mirQTL, batchVars),
                           covMat,
                           emirCount,
                           esnpCount,
                           eqtlCount,
                           genotypePCA,
                           genotypeNumComps,
                           expressionPCA,
                           expressionNumComps,
                           samples),
                    select(df.cov.mirQTLor,
                           covMat,
                           emirCount,
                           esnpCount,
                           eqtlCount,
                           genotypePCA,
                           genotypeNumComps,
                           expressionPCA,
                           expressionNumComps,
                           samples))

df.cov %>%
    ggplot(aes(x=expressionNumComps, y=eqtlCount, color=expressionPCA)) +
    geom_line(mapping = aes(linetype=samples), size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Expression PCs",
         y = "Number of eQTLs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + 10 Genotype PCs + (Expression PCs) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Outlier Removal and Choice of Expression Data")

ggsave(paste0(dir.pdf, "eqtls_by_samples_and_expressionPCA.pdf"), width = 9, height = 5)






