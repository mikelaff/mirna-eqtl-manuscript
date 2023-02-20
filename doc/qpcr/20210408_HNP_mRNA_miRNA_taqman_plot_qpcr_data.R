
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/qpcr/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# mRNA qPCR data from
#mRNA.data.xlsx <- here("results/qpcr/20210325_HNP_mRNA.xlsx")

# miRNA qPCR TaqMan data from 8 Apr 2021
miRNA.taqman.data.xlsx <- here("results/qpcr/20210408_HNP_miRNA_TaqMan.xlsx")


# GLOBALS ##############################################################################################################


# Import mRNA Data #####################################################################################################

df.data <- read_xlsx(mRNA.data.xlsx, sheet = 3, range = "A48:O384", na = c("", "Undetermined"))

df.data %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 3))

df.data %<>%
    select(Sample,
           Donor,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.data$Expression <- factor(df.data$Expression,
                             levels = c("pLenti-Control", "HAUS4", "pTREX-Control", "4707-C", "4707-A", "Mock"),
                             labels = c("pLenti-Control", "HAUS4", "pTREX-Control", "miR-4707-3p-REF(C)", "miR-4707-3p-ALT(A)", "Mock"),
                             ordered = TRUE)

df.data$Donor <- factor(df.data$Donor)
df.data$Replicate <- factor(df.data$Replicate)

df.data %<>%
    mutate(Name = paste(Donor, Expression, Replicate))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
              CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# Filter out outlier samples (determined after first pass analysis of data)
# df.data %<>%
#     filter(!Sample == "D54_4707-A_3",
#            !Sample == "D88_4707-C_2",
#            !Sample == "D54_Control_1",
#            !is.na(CT_mean))


# loop over samples, calculate delta CT to ACTB
samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp

    # get ACTB CT value for this sample
    ct.actb <- df.tmp$CT_mean[match("ACTB", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_ACTB = CT_mean - ct.actb)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.plenti <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and HAUS4/pLenti-Control samples
    df.data %>%
        filter(Donor == donor, Expression == "HAUS4" | Expression == "pLenti-Control" | Expression == "Mock") -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_ACTB[which("pLenti-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB))

        # combine into all data
        df.plenti <- bind_rows(df.plenti, df.target)
    }

}

rm(df.donor, df.target)

# loop over donors, calculate delta delta CT to pTREX-Control of that donor for each target

df.ptrex <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and 4707-A/C/pTREX-Control samples
    df.data %>%
        filter(Donor == donor, Expression == "miR-4707-3p-ALT(A)" | Expression == "miR-4707-3p-REF(C)" | Expression == "pTREX-Control") -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_ACTB[which("pTREX-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB))

        # combine into all data
        df.ptrex <- bind_rows(df.ptrex, df.target)
    }

}

rm(df.donor, df.target)

df.data <- bind_rows(df.plenti, df.ptrex)
rm(df.plenti, df.ptrex)

# Repeat for EIF4A2

# loop over samples, calculate delta CT to EIF4A2
samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp

    # get ACTB CT value for this sample
    ct.eif4a2 <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_EIF4A2 = CT_mean - ct.eif4a2)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.plenti <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and HAUS4/pLenti-Control samples
    df.data %>%
        filter(Donor == donor, Expression == "HAUS4" | Expression == "pLenti-Control" | Expression == "Mock") -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_EIF4A2[which("pLenti-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

        # combine into all data
        df.plenti <- bind_rows(df.plenti, df.target)
    }

}

rm(df.donor, df.target)

# loop over donors, calculate delta delta CT to pTREX-Control of that donor for each target

df.ptrex <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and 4707-A/C/pTREX-Control samples
    df.data %>%
        filter(Donor == donor, Expression == "miR-4707-3p-ALT(A)" | Expression == "miR-4707-3p-REF(C)" | Expression == "pTREX-Control") -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_EIF4A2[which("pTREX-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

        # combine into all data
        df.ptrex <- bind_rows(df.ptrex, df.target)
    }

}

rm(df.donor, df.target)

df.data <- bind_rows(df.plenti, df.ptrex)
rm(df.plenti, df.ptrex)

# Import TaqMan miRNA Data #############################################################################################
df.taqman <- read_xlsx(miRNA.taqman.data.xlsx, sheet = 3, range = "A45:O109", na = c("", "Undetermined"))

df.taqman %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman %<>%
    mutate(Donor = "293T",
           Expression = sapply(strsplit(Sample, "_"), `[`, 1),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 2))

df.taqman %<>%
    select(Sample,
           Donor,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.taqman$Expression <- factor(df.taqman$Expression,
                              levels = c("pLenti-Control", "pLenti-HAUS4", "pLenti-4707-C", "pLenti-4707-A", "pTREX-Control", "pTREX-4707-C", "pTREX-4707-A", "Mock"),
                              labels = c("pLenti-Control", "pLenti-HAUS4", "pLenti-4707-REF(C)", "pLenti-4707-ALT(A)", "pTREX-Control", "pTREX-4707-REF(C)", "pTREX-4707-ALT(A)", "Mock"),
                              ordered = TRUE)

df.taqman$Donor <- factor(df.taqman$Donor)
df.taqman$Replicate <- factor(df.taqman$Replicate)

df.taqman %<>%
    mutate(Name = paste(Donor, Expression, Replicate))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.taqman %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# loop over samples, calculate delta CT to miR-361
samples <- unique(df.taqman$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.taqman %>%
        filter(Sample == sample) -> df.tmp

    # get miR-361 CT value for this sample
    ct.361 <- df.tmp$CT_mean[match("miR-361", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_miR361 = CT_mean - ct.361)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.taqman <- df.new
rm(df.new, df.tmp)

# HAUS4 / pLenti-Control
# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.taqman$Donor)
targets <- unique(df.taqman$Target)

df.plenti <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and HAUS4/pLenti-Control samples
    df.taqman %>%
        filter(Donor == donor, (Expression == "pLenti-HAUS4" | Expression == "pLenti-Control" | Expression == "pLenti-4707-REF(C)" | Expression == "pLenti-4707-ALT(A)" | Expression == "Mock")) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_miR361[which("pLenti-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))

        # combine into all data
        df.plenti <- bind_rows(df.plenti, df.target)
    }

}

# miR-4707-A/C / pTREX-Control
# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.taqman$Donor)
targets <- unique(df.taqman$Target)

df.ptrex <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor and 4707-A/C/pTREX-Control samples
    df.taqman %>%
        filter(Donor == donor, (Expression == "pTREX-4707-ALT(A)" | Expression == "pTREX-4707-REF(C)" | Expression == "pTREX-Control")) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_miR361[which("pTREX-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))

        # combine into all data
        df.ptrex <- bind_rows(df.ptrex, df.target)
    }

}

df.taqman <- bind_rows(df.plenti, df.ptrex)
rm(df.plenti, df.ptrex, df.donor, df.target)

# Plot #################################################################################################################

# my.theme <- theme(axis.title = element_text(size = 16),
#                   axis.text = element_text(size = 14),
#                   axis.text.x = element_text(size = 12),
#                   title = element_text(size = 18),
#                   legend.title = element_text(size = 16),
#                   legend.text = element_text(size = 14),
#                   strip.text = element_text(size = 16),
#                   plot.caption = element_text(size = 12))

#pdf(paste0(dir.pdfs, "20210217_QIAGEN_v_TaqMan.pdf"))

df.data %>%
    ggplot(aes(y = fold_change_ACTB, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = fold_change_ACTB, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(grepl("HAUS", Target) | grepl("Ki67", Target)) %>%
    ggplot(aes(y = fold_change_ACTB, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(grepl("HAUS", Target) | grepl("Ki67", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")


df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    #filter(grepl("HAUS", Target) | grepl("Ki67", Target)) %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "qPCR after HNP Lentivirus Transductions")

df.taqman %>%
    ggplot(aes(y = fold_change_miR361, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: miR-361",
         title = "TaqMan: qPCR after HNP Lentivirus Transductions")

df.taqman %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "293T expression after 2 days post transfection with transfer plasmids.",
         title = "TaqMan: qPCR after 293T Transfections")

ggsave(paste0(dir.pdfs, "293T_miRNA_expression_taqman.pdf"))
#dev.off()

