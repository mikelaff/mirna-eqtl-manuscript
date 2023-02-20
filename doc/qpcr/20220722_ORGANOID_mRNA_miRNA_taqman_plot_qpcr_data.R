
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/qpcr/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# mRNA qPCR from 28 July 2022
mRNA.data1.xlsx <- here("results/qpcr/20220728_MIR4707_ORGANOID_QPCR_MRNA.xlsx")
# miRNA qPCR TaqMan data from 22 July 2022
miRNA.taqman.data1.xlsx <- here("results/qpcr/20220722_MIR4707_ORGANOID_QPCR.xlsx")

# mRNA data for the 13 diff expressed genes from the HNP experiment
mRNA.data2.xlsx <- here("results/qpcr/20220921_MIR4707_QPCR_13_GENES.xlsx")

# GLOBALS ##############################################################################################################

# Import Sample Data ###################################################################################################
df.data1 <- read_xlsx(mRNA.data1.xlsx, sheet = "Results", range = "A36:O228", na = c("", "Undetermined"))

df.data1 %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data1 %<>%
    mutate(Expression = NA)

df.data1$Expression[df.data1$Sample %in% c("1", "2", "3", "4")] <- "pTRIPZ-Control"
df.data1$Expression[df.data1$Sample %in% c("5", "6", "7", "8")] <- "pTRIPZ-4707"

df.data1$Expression <- factor(df.data1$Expression,
                              levels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                              ordered = TRUE)

df.data1$Target <- factor(df.data1$Target)


# df.data %<>%
#     mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))


# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data1 %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# loop over samples, calculate delta CT to ACTB, EIF4A2
samples <- unique(df.data1$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data1 %>%
        filter(Sample == sample) -> df.tmp

    # get CT value for this sample
    ct.actb <- df.tmp$CT_mean[match("ACT-B", df.tmp$Target)]
    ct.eif <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_ACTB = CT_mean - ct.actb,
               delta_CT_EIF4A2 = CT_mean - ct.eif)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)
}

df.data1 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Control of that donor for each target
targets <- unique(df.data1$Target)

df.new <- tibble()


# loop over targets
for (target in targets) {
    df.target <- NULL

    # filter for only this target
    df.data1 %>%
        filter(Target == target) -> df.target

    # use the mean -DOX value across the replicates for delta delta ct
    delta.ct.control.actb <- mean(df.target$delta_CT_ACTB[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
    delta.ct.control.eif4a2 <- mean(df.target$delta_CT_EIF4A2[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)

    # calculate delta delta CT within this target and donor
    df.target %<>%
        mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control.actb,
               delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control.eif4a2)

    # calculate fold change
    df.target %<>%
        mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB),
               fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

    # combine into all data
    df.new <- bind_rows(df.new, df.target)
}



df.data1 <- df.new
rm(df.target, df.new)

# 13 genes
df.data2 <- read_xlsx(mRNA.data2.xlsx, sheet = "Results", range = "A37:O261", na = c("", "Undetermined"))

df.data2 %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data2 %<>%
    mutate(Expression = NA)

df.data2$Expression[df.data2$Sample %in% c("1", "2", "3", "4")] <- "pTRIPZ-Control"
df.data2$Expression[df.data2$Sample %in% c("5", "6", "7", "8")] <- "pTRIPZ-4707"

df.data2$Expression <- factor(df.data2$Expression,
                              levels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                              ordered = TRUE)

df.data2$Target <- factor(df.data2$Target)


# df.data %<>%
#     mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))


# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data2 %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# loop over samples, calculate delta CT to ACTB, EIF4A2
samples <- unique(df.data2$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data2 %>%
        filter(Sample == sample) -> df.tmp

    # get CT value for this sample
    #ct.actb <- df.tmp$CT_mean[match("ACT-B", df.tmp$Target)]
    ct.eif <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(#delta_CT_ACTB = CT_mean - ct.actb,
               delta_CT_EIF4A2 = CT_mean - ct.eif)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)
}

df.data2 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Control of that donor for each target
targets <- unique(df.data2$Target)

df.new <- tibble()


# loop over targets
for (target in targets) {
    df.target <- NULL

    # filter for only this target
    df.data2 %>%
        filter(Target == target) -> df.target

    # use the mean -DOX value across the replicates for delta delta ct
    #delta.ct.control.actb <- mean(df.target$delta_CT_ACTB[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
    delta.ct.control.eif4a2 <- mean(df.target$delta_CT_EIF4A2[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)

    # calculate delta delta CT within this target and donor
    df.target %<>%
        mutate(#delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control.actb,
               delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control.eif4a2)

    # calculate fold change
    df.target %<>%
        mutate(#fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB),
               fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

    # combine into all data
    df.new <- bind_rows(df.new, df.target)
}



df.data2 <- df.new
rm(df.target, df.new)


# Import TaqMan miRNA Data, More Controls ##############################################################################
df.taqman <- read_xlsx(miRNA.taqman.data1.xlsx, sheet = "Results", range = "A35:O163", na = c("", "Undetermined"))

df.taqman %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman %<>%
    mutate(Expression = NA)

df.taqman$Expression[df.taqman$Sample %in% c("1", "2", "3", "4")] <- "pTRIPZ-Control"
df.taqman$Expression[df.taqman$Sample %in% c("5", "6", "7", "8")] <- "pTRIPZ-4707"

df.taqman$Expression <- factor(df.taqman$Expression,
                               levels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                               ordered = TRUE)

df.taqman$Target <- factor(df.taqman$Target,
                           levels = c("197", "324", "361", "4707"),
                           labels = c("miR-197-3p", "miR-324-3p", "miR-361-5p", "miR-4707-3p"),
                           ordered = TRUE)

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.taqman %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

df.taqman %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = .7, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = .7) +
    facet_wrap(~Target, scales = "free_y", nrow = 2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(size = 5, method = "t.test", label = "p.format")



# loop over samples, calculate delta CT to controls
samples <- unique(df.taqman$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.taqman %>%
        filter(Sample == sample) -> df.tmp

    # get miR-361 CT value for this sample
    ct.361 <- df.tmp$CT_mean[match("miR-361-5p", df.tmp$Target)]
    ct.197 <- df.tmp$CT_mean[match("miR-197-3p", df.tmp$Target)]
    ct.324 <- df.tmp$CT_mean[match("miR-324-3p", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_miR361 = CT_mean - ct.361,
               delta_CT_miR197 = CT_mean - ct.197,
               delta_CT_miR324 = CT_mean - ct.324)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.taqman <- df.new
rm(df.new, df.tmp)

# calculate delta delta CT to Control of that sample for each target
targets <- unique(df.taqman$Target)

df.new <- tibble()

# loop over targets
for (target in targets) {
    df.target <- NULL

    # filter for only this target
    df.taqman %>%
        filter(Target == target) -> df.target

    # use the mean control value across the replicates for delta delta ct
    delta.ct.control.361 <- mean(df.target$delta_CT_miR361[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
    delta.ct.control.197 <- mean(df.target$delta_CT_miR197[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
    delta.ct.control.324 <- mean(df.target$delta_CT_miR324[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)

    # calculate delta delta CT within this target and donor
    df.target %<>%
        mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control.361,
               delta_delta_CT_miR197 = delta_CT_miR197 - delta.ct.control.197,
               delta_delta_CT_miR324 = delta_CT_miR324 - delta.ct.control.324)

    # calculate fold change
    df.target %<>%
        mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361),
               fold_change_miR197 = 2 ^ (-delta_delta_CT_miR197),
               fold_change_miR324 = 2 ^ (-delta_delta_CT_miR324))

    # combine into all data
    df.new <- bind_rows(df.new, df.target)
}



df.taqman <- df.new
rm(df.new, df.target)


# Plot #################################################################################################################

df.data1 %>%
    ggplot(aes(y = fold_change_ACTB, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.data1 %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.data1 %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")


pdf(paste0(dir.pdfs, "20220722_Organoid_miRNA_expression.pdf"), height = 7, width = 12)

df.taqman %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR197, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.taqman %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR324, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.taqman %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = log2(fold_change_miR361), x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")


df.taqman %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

dev.off()

df.data2 %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.data2 %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    facet_wrap(~Target, nrow = 1) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          legend.position = "bottom") +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Expression",
         title = "Organoid Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")



