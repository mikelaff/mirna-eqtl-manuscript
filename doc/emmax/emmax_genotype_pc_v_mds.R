
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(mikelaffr)

# OUTPUT ##############
dir.pdf <- here("doc/emmax/pdfs/")

# INPUT ############
files.rds <- list.files(here("results/emmax/association_results/20191115_test_covariates_genotype_pc_v_mds/"), pattern = "*.rds")

# Load Results Summary Tables ##################
# combine all summary tables into one
df.eqtls <- tibble()

for (df.file in files.rds) {
    print(df.file)

    df.tmp <- as_tibble(readRDS(paste0(here("results/emmax/association_results/20191115_test_covariates_genotype_pc_v_mds/"), df.file)))

    df.eqtls <- bind_rows(df.eqtls, df.tmp)

}

rm(df.tmp)

df.eqtls

df.eqtls$eqtl <- paste(df.eqtls$emir, df.eqtls$esnp, sep = "+")

# df.emirs %>%
#     group_by(covMat) %>%
#     summarise(emirCount = n()) -> df.cov

df.eqtls$genotypeType <- sapply(strsplit(df.eqtls$covMat, "\\."), `[`, 2)
df.eqtls$numComps <- as.numeric(sapply(strsplit(df.eqtls$genotypeType, "_"), `[`, 2))
df.eqtls$genotypeType <- sapply(strsplit(df.eqtls$genotypeType, "_"), `[`, 1)

df.eqtls %>%
    group_by(covMat) %>%
    summarise(emirCount = n_distinct(emir), esnpCount = n_distinct(esnp), eqtlCount = n_distinct(eqtl)) -> df.cov

df.cov$genotypeType <- sapply(strsplit(df.cov$covMat, "\\."), `[`, 2)
df.cov$numComps <- as.numeric(sapply(strsplit(df.cov$genotypeType, "_"), `[`, 2))
df.cov$genotypeType <- sapply(strsplit(df.cov$genotypeType, "_"), `[`, 1)

df.cov %>%
    ggplot(aes(x=numComps, y=emirCount, color=genotypeType)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Components",
         y = "Number of emiRs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + (MDS or PCA) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusions of Genotype PCs or MDS Comp.")

df.cov %>%
    ggplot(aes(x=numComps, y=esnpCount, color=genotypeType)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Components",
         y = "Number of eSNPs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + (MDS or PCA) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusions of Genotype PCs or MDS Comp.")

df.cov %>%
    ggplot(aes(x=numComps, y=eqtlCount, color=genotypeType)) +
    geom_line(size = 1) +
    scale_x_continuous(breaks = seq(0,20,1)) +
    scale_color_manual(values=cbPalette) +
    theme(panel.grid.minor.x = element_blank()) +
    labs(x = "Number of Components",
         y = "Number of eQTLs (FDR < 0.05)",
         caption = "Expression ~ Genotype + Kinship + (MDS or PCA) + Seq Pool + Purif Method + RIN + Sex + Gest Week",
         title = "Inclusions of Genotype PCs or MDS Comp.")

ggsave(paste0(dir.pdf, "eqtls_by_genoPCs_or_MDS.pdf"), width = 9, height = 5)







