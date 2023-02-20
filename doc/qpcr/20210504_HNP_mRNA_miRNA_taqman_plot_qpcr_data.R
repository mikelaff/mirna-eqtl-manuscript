
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
# mRNA qPCR data from 6 May 2021
mRNA.data.xlsx <- here("results/qpcr/20210506_HNP_mRNA.xlsx")

# miRNA qPCR TaqMan data from 4 May 2021
miRNA.taqman.data.xlsx <- here("results/qpcr/20210504_HNP_miRNA_TaqMan.xlsx")


# GLOBALS ##############################################################################################################


# Import mRNA Data #####################################################################################################

df.data <- read_xlsx(mRNA.data.xlsx, sheet = 3, range = "A46:O158", na = c("", "Undetermined"))

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
                             levels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                             labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
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

# loop over donors, calculate delta delta CT to pTRIPZ-Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_ACTB[which("pTRIPZ-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB))

        # combine into all data
        df.tripz <- bind_rows(df.tripz, df.target)
    }

}

rm(df.donor, df.target)
df.data <- df.tripz

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

# loop over donors, calculate delta delta CT to pTRIPZ-Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_EIF4A2[which("pTRIPZ-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

        # combine into all data
        df.tripz <- bind_rows(df.tripz, df.target)
    }

}

rm(df.donor, df.target)
df.data <- df.tripz

rm(df.tripz)

# Import TaqMan miRNA Data #############################################################################################
df.taqman <- read_xlsx(miRNA.taqman.data.xlsx, sheet = 3, range = "A45:O125", na = c("", "Undetermined"))

df.taqman %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 3))

df.taqman %<>%
    select(Sample,
           Donor,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.taqman$Expression <- factor(df.taqman$Expression,
                              levels = c("pLenti-Control", "pLenti-4707-C", "pTREX-Control", "pTREX-4707-C", "pTRIPZ-Control", "pTRIPZ-4707-C"),
                              labels = c("pLenti-Control", "pLenti-4707-C", "pTREX-Control", "pTREX-4707-C", "pTRIPZ-Control", "pTRIPZ-4707-C"),
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

df.taqman.293 <- filter(df.taqman, Donor == "293T")
df.taqman.D88 <- filter(df.taqman, Donor == "D88")


# pTRIPZ for D88
# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.taqman.D88$Donor)
targets <- unique(df.taqman.D88$Target)

df.ptripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.taqman.D88 %>%
        filter(Donor == donor) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_miR361[which("pTRIPZ-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))

        # combine into all data
        df.ptripz <- bind_rows(df.ptripz, df.target)
    }

}

df.taqman.D88 <- df.ptripz

# pTRIPZ for 293
# loop over donors, calculate delta delta CT to pLenti-Control of that donor for each target
donors <- unique(df.taqman.293$Donor)
targets <- unique(df.taqman.293$Target)

df.ptripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.taqman.293 %>%
        filter(Donor == donor, grepl("pTRIPZ", Expression)) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_miR361[which("pTRIPZ-Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))

        # combine into all data
        df.ptripz <- bind_rows(df.ptripz, df.target)
    }

}

df.taqman.293 <- df.ptripz

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

df.taqman.D88 %>%
    ggplot(aes(y = fold_change_miR361, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), mapping = aes(shape = Replicate)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: miR-361",
         title = "TaqMan: qPCR after HNP Lentivirus Transductions")

ggsave("D88_pTRIPZ_fold_change.pdf", height = 4, width = 5)


df.taqman.D88 %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "TaqMan: qPCR after HNP Lentivirus Transductions")

ggsave("D88_pTRIPZ_mean_ct.pdf", height = 4, width = 5)

df.taqman.293 %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "293T expression after 4 days.",
         title = "TaqMan: qPCR after HNP Lentivirus Transduction of 293T")

ggsave("293t_pTRIPZ_mean_ct.pdf", height = 4, width = 6)

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
    facet_wrap(~Target, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("fold_change_actb.pdf", height = 5, width = 8)

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
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("fold_change_eif4a2.pdf", height = 5, width = 8)


df.data %>%
    #filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
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
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 1) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target, nrow = 1, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("ct_mean.pdf", height = 5, width = 8)


#dev.off()

