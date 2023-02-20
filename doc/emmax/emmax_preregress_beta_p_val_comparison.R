# plot index snps and snps in ld with index

library(here)
library(dplyr)
library(tidyr)
library(readr)
library(magrittr)
library(ggplot2)
#library(Gviz)
#library(gridExtra)
#library(biomaRt)
library(DESeq2)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/")
dir.pngs <- here("doc/emmax/pngs/")

# summary file of computed betas and p-values
dir.output <- here("doc/emmax/rdata/")
dir.create(dir.output, recursive = TRUE, showWarnings = FALSE)

output.rds <- paste0(dir.output, "preregress.eqtl.ld.snps.dataFrame.rds")

# INPUT ################################################################################################################
# eSNPs, emiRs, and eQTLs for mirQTLor association analysis
eqtls.dataframe.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor.eQTLs.dataFrame.rds")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_dataFrame.rds")

# nominal p-value threshold (FDR < 0.05) used for clumping and plotting
nom.p.val.txt <- here("results/emmax/association_results/20200103_mirQTLor_preregress/compiled/20200103_mirQTLor_nomPvalue.txt")

# summarized experiment with vst expression values used in this analysis
vsd.rds <- here("results/emmax/phenotype_files/20200103_mirQTLor_VST_miRNA_expression_preregress/20200103_mirQTLor_VST_miRNA_expression_preregress_rse.rds")

# Directory for sample genotypes at each index SNP
dir.genotypes <- here("results/emmax/association_results/20200103_mirQTLor_preregress/sample_genotypes/")

# Root directory for prefiltered hardcall genotypes per mirna
dir.root.hardcalls <- here("results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/")
# Root directory for prefiltered dosage genotyeps per mirna
dir.root.dosages <- here("results/emmax/tfiles/prefiltered_mirQTLor_dosages/")

# Directory for LD at each index SNP
dir.ld <- here("results/emmax/association_results/20200103_mirQTLor_preregress/ld/")

# GLOBALS ##############################################################################################################

# Import Summarized Results ############################################################################################

# compiled eQTLs
df.eqtls <- as_tibble(readRDS(eqtls.dataframe.rds))
df.eqtls %<>%
    mutate(eqtl = paste(emir, esnp, sep = "+"))

# all variants
df.results <- as_tibble(readRDS(summarized.results.dataframe.rds))
df.results %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))
#filter for only hsa-mir-101-1 and hsa-mir-101-2 results
#df.results %<>%
#    dplyr::filter(!grepl("hsa-mir-101-", UniName))

# expression data
vsd <- readRDS(vsd.rds)

# nominal p-value for plotting
nom.p.val <- as.numeric(read_lines(nom.p.val.txt))

# Investigation Plots ##################################################################################################

# tibble to save comparison results
df.beta.p <- tibble()

for (k in 1:nrow(df.eqtls)) {

    printMessage(paste("Comparing", k, "of", nrow(df.eqtls), ":", df.eqtls$eqtl[k]))

    eqtl.ind <- k

    esnp <- df.eqtls$esnp[eqtl.ind]
    emir <- df.eqtls$emir[eqtl.ind]
    eqtl <- df.eqtls$eqtl[eqtl.ind]

    chr <- strsplit(esnp, ":")[[1]][1]

    # genotype hardcall file: chr1.hardcall.prefiltered.mirQTLor.hsa-mir-101-1_hsa-miR-101-3p.traw
    geno.hardcall.file <- paste0(dir.root.hardcalls, chr, "/", chr, ".hardcall.prefiltered.mirQTLor.", emir, ".traw")
    df.genotype.hardcalls <- read_tsv(geno.hardcall.file, col_types = cols())

    # genotype dosage file: chr1.dosage.prefiltered.mirQTLor.hsa-mir-101-1_hsa-miR-101-3p.traw
    geno.dosage.file <- paste0(dir.root.dosages, chr, "/", chr, ".dosage.prefiltered.mirQTLor.", emir, ".traw")
    df.genotype.dosages <- read_tsv(geno.dosage.file, col_types = cols())

    # subset results for this miRNA
    df.results.emir <- dplyr::filter(df.results, UniName == emir)

    # subset expression data
    df.expression <- data.frame(expression = assays(vsd[emir,])[["VST"]][1,],
                                expression.corrected = assays(vsd[emir,])[["VST.RESIDUAL"]][1,],
                                DonorID = paste0("D", vsd$donor_id),
                                stringsAsFactors = FALSE)
    df.expression$RNAID <- rownames(df.expression)

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")
    # ld for this esnp
    df.ld <- read_table(ld.file, col_types = cols())
    df.ld %<>%
        dplyr::select(SNP = SNP_B, R2)

    # add ld to results table
    df.results.emir %<>%
        left_join(df.ld, by = "SNP")

    # dot color based on LD
    df.results.emir$color <- "navy"
    df.results.emir$color[which(df.results.emir$R2 >= 0.8)] <- "red"
    df.results.emir$color[which(df.results.emir$R2 >= 0.6 & df.results.emir$R2 < 0.8)] <- "orange"
    df.results.emir$color[which(df.results.emir$R2 >= 0.4 & df.results.emir$R2 < 0.6)] <- "green"
    df.results.emir$color[which(df.results.emir$R2 >= 0.2 & df.results.emir$R2 < 0.4)] <- "lightblue"
    df.results.emir$color[which(df.results.emir$SNP == esnp)] <- "purple"

    # loop over variants in high LD with index (red variants)
    df.results.emir %>%
        dplyr::filter(!is.na(R2)) %>%
        arrange(desc(R2)) -> df.vars

    # columns for lm beta and p-value
    df.vars$BETA.LM <- NA
    df.vars$P.LM <- NA

    # if (nrow(df.vars) > 10) {
    #     df.vars <- df.vars[1:10,]
    # }

    for (i in 1:nrow(df.vars)) {

        # this variant
        this.var <- df.vars$SNP[i]

        printMessage(paste("Modeling", i, "of", nrow(df.vars), ":", this.var),
                     fillChar = "-")

        # reference allele
        allele.ref <- df.vars$REF[match(this.var, df.vars$SNP)]
        # alt allele (effect allele)
        allele.alt <- df.vars$ALT[match(this.var, df.vars$SNP)]
        # minor allele
        allele.minor <- df.vars$A1[match(this.var, df.vars$SNP)]
        # major allele
        allele.major <- df.vars$A2[match(this.var, df.vars$SNP)]
        # effect allele (for dosage based associations)
        effect.allele<- df.vars$EFFECT.ALLELE[match(this.var, df.vars$SNP)]

        if (effect.allele != allele.ref) {
            stop("Effect allele != reference allele!")
        }


        # hardcall genotype counted allele
        counted.hardcall.allele <- df.genotype.hardcalls$COUNTED[which(df.genotype.hardcalls$SNP == this.var)]
        # genotypes for this variant
        df.genotype.hardcalls %>%
            dplyr::filter(SNP == this.var) %>%
            dplyr::select(starts_with("D")) %>%
            pivot_longer(cols = everything(),
                         names_to = "DonorID_DNAID",
                         values_to = "Genotype_Hard") -> df.hardcalls

        # hardcall genotype labels
        if (counted.hardcall.allele == effect.allele) {
            hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                    paste(allele.alt, allele.ref, sep = "/"),
                                    paste(allele.ref, allele.ref, sep = "/"))
            df.hardcalls$Genotype_Hard <- factor(df.hardcalls$Genotype_Hard, levels = c(0, 1, 2), ordered = TRUE, labels = hardcall.genotypes)
        } else if (counted.hardcall.allele == allele.alt) {
            hardcall.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                    paste(allele.alt, allele.ref, sep = "/"),
                                    paste(allele.ref, allele.ref, sep = "/"))
            # format genotype hardcalls as factor
            df.hardcalls$Genotype_Hard <- factor(df.hardcalls$Genotype_Hard, levels = rev(c(0, 1, 2)), ordered = TRUE, labels = hardcall.genotypes)
        } else {
            stop("Something is wrong with genotype labels!")
        }

        # dosage genotype counted allele
        counted.dosage.allele <- df.genotype.dosages$COUNTED[which(df.genotype.dosages$SNP == this.var)]
        # genotypes for this variant
        df.genotype.dosages %>%
            dplyr::filter(SNP == this.var) %>%
            dplyr::select(starts_with("D")) %>%
            pivot_longer(cols = everything(),
                         names_to = "DonorID_DNAID",
                         values_to = "Genotype_Dose") -> df.dosages

        # dosage genotype check
        if (counted.dosage.allele == effect.allele) {
            dosage.genotypes <- c(paste(allele.alt, allele.alt, sep = "/"),
                                  paste(allele.alt, allele.ref, sep = "/"),
                                  paste(allele.ref, allele.ref, sep = "/"))
        } else {
            stop("Something is wrong with genotype labels!")
        }

        # combine into single genotype dataframe
        df.genotypes <- left_join(df.dosages, df.hardcalls, by = "DonorID_DNAID")

        rm(df.dosages, df.hardcalls)

        # make donor id column to join with expression table
        df.genotypes$DonorID <- sapply(strsplit(df.genotypes$DonorID_DNAID, "_"), `[`, 1)

        # join genotypes with expression data
        df.expression %>%
            left_join(df.genotypes, by = "DonorID") -> df.expression.this.var

        # linear model: corrected_expression ~ genotype_dose
        model.fit <- lm(formula = expression.corrected ~ Genotype_Dose, data = df.expression.this.var)

        #lm.intercept <- model.fit$coefficients["(Intercept)"]
        df.vars$BETA.LM[i] <- model.fit$coefficients["Genotype_Dose"]
        df.vars$P.LM[i] <- summary(model.fit)$coefficients[,4]["Genotype_Dose"]

    }

    # add df.vars to compilation table
    df.beta.p <- bind_rows(df.beta.p, df.vars)

}


# save data frame
saveRDS(df.beta.p, output.rds)

# load data
df.beta.p <- readRDS(here("doc/emmax/rdata/preregress.eqtl.ld.snps.dataFrame.rds"))

df.beta.p %<>%
    mutate(BETA.diff = abs(BETA - BETA.LM), P.diff = abs(P - P.LM))

summary(df.beta.p$BETA.diff)
summary(df.beta.p$P.diff)

df.beta.p %>%
    dplyr::filter(P <= nom.p.val) -> df.tmp

summary(df.tmp$BETA.diff)
summary(df.tmp$P.diff)

df.beta.p %<>%
    mutate(SIG.EMMAX = P <= nom.p.val, SIG.LM = P.LM <= nom.p.val)

df.beta.p$SIG.EMMAX.ONLY <- df.beta.p$SIG.EMMAX & (! df.beta.p$SIG.LM)
df.beta.p$SIG.LM.ONLY <- df.beta.p$SIG.LM & (! df.beta.p$SIG.EMMAX)
df.beta.p$SIG.BOTH <- df.beta.p$SIG.EMMAX & df.beta.p$SIG.LM

df.beta.p$SIG <- factor(ifelse(df.beta.p$SIG.EMMAX.ONLY, "EMMAX.ONLY",
                        ifelse(df.beta.p$SIG.LM.ONLY, "LM.ONLY",
                               ifelse(df.beta.p$SIG.BOTH, "BOTH", "NEITHER"))))

sum(df.beta.p$SIG.EMMAX)
sum(df.beta.p$SIG.LM)



df.beta.p %>%
    ggplot(aes(x = BETA.diff, fill = SIG)) +
    geom_histogram(bins = 10)

df.beta.p %>%
    filter(!SIG == "NEITHER") %>%
    ggplot(aes(x = BETA.diff, fill = SIG)) +
    geom_histogram(bins = 20) +
    scale_fill_manual(values = cbPalette[c(1,2,4)]) +
    labs(x = "Abs. Difference in Beta Values",
         title = "Significant miR-eQTL Associations using EMMAX or LM",
         fill = "Significant\nAssociations")

ggsave(paste0(dir.pdfs, "emmax_v_lm_beta_values.pdf"), height = 5, width = 7)



df.beta.p %>%
    filter(!SIG == "NEITHER") %>%
    ggplot(aes(x = P.diff, fill = SIG)) +
    geom_histogram(bins = 10) +
    scale_y_log10()

