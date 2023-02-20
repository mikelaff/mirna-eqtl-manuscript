# look at the variance explained by the miRNA-eQTLs

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(ggplot2)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/mirQTL_eda/pdfs/")

# INPUT ################################################################################################################
# mirna-eqtls
eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# GLOBALS ##############################################################################################################
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

# Import eQTL Data ###############################################################################################

df.eqtls <- read_rds(eqtls.rds)

# conditional rank, for adding to degrees of freedom
df.eqtls %<>%
    mutate(RANK = as.integer(factor(DEGREE,
                                    levels = c("primary", "secondary", "tertiary", "quarternary", "quinary"),
                                    labels = c(1,2,3,4,5))))

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

# compute DF, add one because primary (rank 1) doesnt need additonal DF
df.eqtls %<>%
    mutate(DF = primaryDF - RANK + 1)


# compute MAF
df.eqtls %<>%
    mutate(MAF = ALT_CTS / OBS_CT)

df.eqtls %<>%
    mutate(MAF = ifelse(MAF > 0.5, 1 - MAF, MAF))


# compute standard error
df.eqtls %<>%
    mutate(SE_beta = getSE(BETA, P, num.samples, DF))


# compute proportion of variance explained
df.eqtls %<>%
    mutate(PVE = (2 * (BETA^2) * MAF * (1 - MAF)) /
               ((2 * (BETA^2) * MAF * (1 - MAF)) + ((SE_beta^2) * 2 * num.samples * MAF * (1 - MAF)))
    )

df.eqtls %>%
    filter(SIGNIFICANCE == "eigenMT_fdr5percent", DEGREE == "primary") %>%
    ggplot(aes(x = PVE*100)) +
    geom_histogram(bins = 10) +
    stat_bin(bins=10, geom="text", aes(label=..count..), vjust=-1.5) +
    plotTheme("figure") +
    scale_x_continuous(limits = c(0,100))
# labs(x = "Prop. of variance explained (PVE)",
#      title = "Local-miRNA-eQTLs (stringent)")

ggsave("~/Desktop/pve.pdf", height = 2, width = 2)


df.eqtls %>%
    filter(SIGNIFICANCE == "eigenMT_fdr5percent") %>%
    group_by(eSNP) %>%
    summarise(n_emiRs = sum(!duplicated(SEQ))) %>%
    ggplot(aes(x = n_emiRs)) +
    geom_histogram(binwidth = 1) +
    stat_bin(binwidth=1, geom="text", aes(label=..count..), vjust=-1.5) +
    scale_x_continuous(breaks = c(1,2,3,4,5)) +
    # labs(x = "Number of emiRs per eSNP",
    #           title = "Local-miRNA-eQTLs (stringent)") +
    plotTheme("figure")

ggsave("~/Desktop/n_emirs.pdf", height = 2, width = 2)









