
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
# mRNA qPCR data from 5 Feb 2021
mRNA.data.xlsx <- here("results/qpcr/20210205_HNP_mRNA.xlsx")
# miRNA qPCR data from 7 Feb 2021
miRNA.data.xlsx <- here("results/qpcr/20210207_HNP_miRNA.xlsx")


# GLOBALS ##############################################################################################################


# Import mRNA Data #####################################################################################################

df.data <- read_xlsx(mRNA.data.xlsx, sheet = 3, range = "A46:O430", na = c("", "Undetermined"))

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
                             levels = c("Control", "HAUS4", "4707-C", "4707-A"),
                             labels = c("Control", "HAUS4", "miR-4707-3p-REF(C)", "miR-4707-3p-ALT(A)"),
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

# loop over donors, calculate delta delta CT to Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.new <- tibble()

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
        delta.ct.control <- mean(df.target$delta_CT_ACTB[which("Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB))

        # combine into all data
        df.new <- bind_rows(df.new, df.target)
    }

}

df.data <- df.new
rm(df.new, df.donor, df.target)

# Repeat for EIF4A2

# loop over samples, calculate delta CT to EIF4A2
samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp

    # get EIF4A2 CT value for this sample
    ct.eif4a2 <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_EIF4A2 = CT_mean - ct.eif4a2)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Control of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.new <- tibble()

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
        delta.ct.control <- mean(df.target$delta_CT_EIF4A2[which("Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

        # combine into all data
        df.new <- bind_rows(df.new, df.target)
    }

}

df.data <- df.new
rm(df.new, df.donor, df.target)

# Import miRNA Data ####################################################################################################
df.mirna <- read_xlsx(miRNA.data.xlsx, sheet = 3, range = "A46:O118", na = c("", "Undetermined"))

df.mirna %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.mirna %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 3))

df.mirna %<>%
    select(Sample,
           Donor,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.mirna$Expression <- factor(df.mirna$Expression,
                             levels = c("Control", "HAUS4", "4707-C", "4707-A"),
                             labels = c("Control", "HAUS4", "miR-4707-3p-REF(C)", "miR-4707-3p-ALT(A)"),
                             ordered = TRUE)

df.mirna$Donor <- factor(df.mirna$Donor)
df.mirna$Replicate <- factor(df.mirna$Replicate)

df.mirna %<>%
    mutate(Name = paste(Donor, Expression, Replicate))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.mirna %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# loop over samples, calculate delta CT to miR-361
samples <- unique(df.mirna$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.mirna %>%
        filter(Sample == sample) -> df.tmp

    # get miR-361 CT value for this sample
    ct.361 <- df.tmp$CT_mean[match("miR-361", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_miR361 = CT_mean - ct.361)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.mirna <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Control of that donor for each target
donors <- unique(df.mirna$Donor)
targets <- unique(df.mirna$Target)

df.new <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.mirna %>%
        filter(Donor == donor) -> df.donor

    # loop over targets
    for (target in targets) {
        df.target <- NULL

        # filter for only this target
        df.donor %>%
            filter(Target == target) -> df.target

        # use the mean control value across the replicates for delta delta ct
        delta.ct.control <- mean(df.target$delta_CT_miR361[which("Control" == df.target$Expression)])

        # calculate delta delta CT within this target and donor
        df.target %<>%
            mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)

        # calculate fold change
        df.target %<>%
            mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))

        # combine into all data
        df.new <- bind_rows(df.new, df.target)
    }

}

df.mirna <- df.new
rm(df.new, df.donor, df.target)


# Plot #################################################################################################################

my.theme <- theme(axis.title = element_text(size = 16),
                  axis.text = element_text(size = 14),
                  axis.text.x = element_text(size = 12),
                  title = element_text(size = 18),
                  legend.title = element_text(size = 16),
                  legend.text = element_text(size = 14),
                  strip.text = element_text(size = 16),
                  plot.caption = element_text(size = 12))

pdf(paste0(dir.pdfs, "20200205_HNP_expression.pdf"))

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

ggsave("~/Desktop/all_with_actb.pdf", width = 12, height = 9)


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
    filter(grepl("HAUS4", Target)) %>%
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
    filter(grepl("HAUS2", Target)) %>%
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
    filter(grepl("HAUS7", Target)) %>%
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
    filter(grepl("Ki67", Target)) %>%
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
    filter(grepl("HAUS4", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/haus4.pdf", width = 12, height = 6)


df.data %>%
    filter(grepl("HAUS2", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/haus2.pdf", width = 8, height = 6)


df.data %>%
    filter(grepl("HAUS7", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/haus7.pdf", width = 12, height = 6)


df.data %>%
    filter(grepl("Ki67", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/ki67.pdf", width = 12, height = 6)


df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/mean_ct_controls.pdf", width = 12, height = 6)

df.mirna %>%
    ggplot(aes(y = fold_change_miR361, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 4 days. Endogenous Control: miR-361",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/mirna.pdf", width = 12, height = 6)


df.mirna %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    my.theme +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 4 days.",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave("~/Desktop/mean_ct_mirna.pdf", width = 12, height = 6)


dev.off()
ggsave("~/Desktop/miR4707_expression_qPCR.pdf", height = 5, width = 8)

