
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
# mRNA qPCR data from 22 Nov 2021, set 1
mRNA.data1.xlsx <- here("results/qpcr/20211122_HNP_mRNA-1.xlsx")

# mRNA qPCR data from 22 Nov 2021, set 2
mRNA.data2.xlsx <- here("results/qpcr/20211122_HNP_mRNA-2.xlsx")

# miRNA qPCR TaqMan data from 19 Nov 2021
miRNA.taqman.data.xlsx <- here("results/qpcr/20211119_HNP_miRNA_TaqMan.xlsx")

# GLOBALS ##############################################################################################################

# Import mRNA Data  Set 1 ##############################################################################################

df.data <- read_xlsx(mRNA.data1.xlsx, sheet = 3, range = "A46:O430", na = c("", "Undetermined"))

df.data %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           DOX = sapply(strsplit(Sample, "_"), `[`, 3),
           Timepoint = sapply(strsplit(Sample, "_"), `[`, 4),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 5))

df.data %<>%
    select(Sample,
           Donor,
           Timepoint,
           Expression,
           DOX,
           Replicate,
           Well,
           Target,
           CT)

df.data$Expression <- factor(df.data$Expression,
                             levels = c("Control", "4707"),
                             labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                             ordered = TRUE)

df.data$DOX <- factor(df.data$DOX,
                      levels = c("+DOX", "-DOX"),
                      labels = c(TRUE, FALSE),
                      ordered = FALSE)

df.data$Timepoint <- factor(df.data$Timepoint,
                            levels = c("1week", "2week"),
                            labels = c("Week 1", "Week 2"),
                            ordered = TRUE)

df.data$Donor <- factor(df.data$Donor)
df.data$Replicate <- factor(df.data$Replicate)

df.data %<>%
    mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))


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

# loop over donors, calculate delta delta CT to -DOX of that donor for each target
donors <- unique(df.data$Donor)
timepoints <- unique(df.data$Timepoint)
targets <- unique(df.data$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.donor

    for (timepoint in timepoints) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Timepoint == timepoint) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
                filter(Target == target) -> df.target

            # use the mean -DOX value across the replicates for delta delta ct
            delta.ct.control.actb <- mean(df.target$delta_CT_ACTB[which(df.target$DOX == FALSE)])
            delta.ct.control.eif4a2 <- mean(df.target$delta_CT_EIF4A2[which(df.target$DOX == FALSE)])

            # calculate delta delta CT within this target and donor
            df.target %<>%
                mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control.actb,
                       delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control.eif4a2)

            # calculate fold change
            df.target %<>%
                mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB),
                       fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

            # combine into all data
            df.tripz <- bind_rows(df.tripz, df.target)
        }
    }
}

df.data <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

# Import mRNA Data  Set 2 ##############################################################################################

df.data2 <- read_xlsx(mRNA.data2.xlsx, sheet = 3, range = "A46:O430", na = c("", "Undetermined"))

df.data2 %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data2 %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           DOX = sapply(strsplit(Sample, "_"), `[`, 3),
           Timepoint = sapply(strsplit(Sample, "_"), `[`, 4),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 5))

df.data2 %<>%
    select(Sample,
           Donor,
           Timepoint,
           Expression,
           DOX,
           Replicate,
           Well,
           Target,
           CT)

df.data2$Expression <- factor(df.data2$Expression,
                              levels = c("Control", "4707"),
                              labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                              ordered = TRUE)

df.data2$DOX <- factor(df.data2$DOX,
                       levels = c("+DOX", "-DOX"),
                       labels = c(TRUE, FALSE),
                       ordered = FALSE)

df.data2$Timepoint <- factor(df.data2$Timepoint,
                             levels = c("1week", "2week"),
                             labels = c("Week 1", "Week 2"),
                             ordered = TRUE)

df.data2$Donor <- factor(df.data2$Donor)
df.data2$Replicate <- factor(df.data2$Replicate)

df.data2 %<>%
    mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))


# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data2 %<>%
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
samples <- unique(df.data2$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data2 %>%
        filter(Sample == sample) -> df.tmp

    # get ACTB CT value for this sample
    ct.actb <- df.tmp$CT_mean[match("ACTB", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_ACTB = CT_mean - ct.actb)

    # get EIF4A2 CT value for this sample
    ct.eif4a2 <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_EIF4A2 = CT_mean - ct.eif4a2)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)
}

df.data2 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to -DOX of that donor for each target
donors <- unique(df.data2$Donor)
timepoints <- unique(df.data2$Timepoint)
targets <- unique(df.data2$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data2 %>%
        filter(Donor == donor) -> df.donor

    for (timepoint in timepoints) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Timepoint == timepoint) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
                filter(Target == target) -> df.target

            # use the mean -DOX value across the replicates for delta delta ct
            delta.ct.control.actb <- mean(df.target$delta_CT_ACTB[which(df.target$DOX == FALSE)])
            delta.ct.control.eif4a2 <- mean(df.target$delta_CT_EIF4A2[which(df.target$DOX == FALSE)])

            # calculate delta delta CT within this target and donor
            df.target %<>%
                mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control.actb,
                       delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control.eif4a2)

            # calculate fold change
            df.target %<>%
                mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB),
                       fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2))

            # combine into all data
            df.tripz <- bind_rows(df.tripz, df.target)
        }
    }
}

df.data2 <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

# Import TaqMan miRNA Data #############################################################################################
df.taqman <- read_xlsx(miRNA.taqman.data.xlsx, sheet = 3, range = "A45:O141", na = c("", "Undetermined"))

df.taqman %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           DOX = sapply(strsplit(Sample, "_"), `[`, 3),
           Timepoint = sapply(strsplit(Sample, "_"), `[`, 4),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 5))

df.taqman %<>%
    select(Sample,
           Donor,
           Timepoint,
           Expression,
           DOX,
           Replicate,
           Well,
           Target,
           CT)

df.taqman$Expression <- factor(df.taqman$Expression,
                               levels = c("Control", "4707"),
                               labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                               ordered = TRUE)

df.taqman$DOX <- factor(df.taqman$DOX,
                        levels = c("+DOX", "-DOX"),
                        labels = c(TRUE, FALSE),
                        ordered = FALSE)

df.taqman$Timepoint <- factor(df.taqman$Timepoint,
                              levels = c("1week", "2week"),
                              labels = c("Week 1", "Week 2"),
                              ordered = TRUE)

df.taqman$Donor <- factor(df.taqman$Donor)
df.taqman$Replicate <- factor(df.taqman$Replicate)

df.taqman %<>%
    mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))

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
    ct.361 <- df.tmp$CT_mean[match("miR-361-5p", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_miR361 = CT_mean - ct.361)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.taqman <- df.new
rm(df.new, df.tmp)

# pTRIPZ
# loop over donors, calculate delta delta CT to Control of that donor for each target
donors <- unique(df.taqman$Donor)
timepoints <- unique(df.taqman$Timepoint)
targets <- unique(df.taqman$Target)

df.ptripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.taqman %>%
        filter(Donor == donor) -> df.donor

    for (timepoint in timepoints) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Timepoint == timepoint) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
                filter(Target == target) -> df.target

            # use the mean control value across the replicates for delta delta ct
            delta.ct.control <- mean(df.target$delta_CT_miR361[which(df.target$DOX == FALSE)])

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
}

df.taqman <- df.ptripz
rm(df.donor, df.day, df.target, df.ptripz)


# Plot #################################################################################################################

# my.theme <- theme(axis.title = element_text(size = 16),
#                   axis.text = element_text(size = 14),
#                   axis.text.x = element_text(size = 12),
#                   title = element_text(size = 18),
#                   legend.title = element_text(size = 16),
#                   legend.text = element_text(size = 14),
#                   strip.text = element_text(size = 16),
#                   plot.caption = element_text(size = 12))

pdf(paste0(dir.pdfs, "20211123_Diff_1_2_week_D54_D88_qPCR_foldChange.pdf"), height = 8, width = 8)

df.taqman %>%
    filter(Target == "miR-4707-3p-C") %>%
    mutate(DOX = factor(DOX, levels = c(FALSE, TRUE), ordered = TRUE)) %>%
    ggplot(aes(y = fold_change_miR361, x = Timepoint, color = DOX, label = Replicate)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Donor + Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme("figure") +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Week") +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "Diff_1-2week_miRNA_fold_change_withSig.pdf"), height = 2, width = 2.5)

df.data %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data2 %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = fold_change_ACTB, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data2 %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = fold_change_ACTB, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")


df.data %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data2 %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = fold_change_ACTB, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.data2 %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = fold_change_ACTB, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    # labs(y = "Fold Change",
    #      x = "Day",
    #      caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
    #      title = "qPCR after HNP Lentivirus Transductions")
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")


dev.off()


pdf(paste0(dir.pdfs, "20211123_Diff_1_2_week_D54_D88_qPCR_meanCT.pdf"), height = 8, width = 8)

df.taqman %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target + Donor) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    #my.theme +
    scale_color_manual(values = c("darkblue", "darkorange", "gray"))

#ggsave(paste0(dir.pdfs, "D88_miRNA_mean_ct.pdf"), height = 4, width = 7)

df.data %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.data2 %>%
    filter(Donor == "D54") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.data %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.data2 %>%
    filter(Donor == "D88") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

dev.off()








df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = fold_change_ACTB, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(grepl("HAUS", Target) | grepl("Ki67", Target) | grepl("CCND1", Target)) %>%
    ggplot(aes(y = fold_change_ACTB, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: ACTB",
         title = "qPCR after HNP Lentivirus Transductions")

#ggsave(paste0(dir.pdfs, "D88_fold_change_actb.pdf"), height = 5, width = 9)

df.data %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

df.data2 %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.6) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    filter(grepl("HAUS", Target) | grepl("Ki67", Target) | grepl("CCND1", Target)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

#ggsave(paste0(dir.pdfs, "D88_fold_change_eif4a2_repeat.pdf"), height = 5, width = 9)

df.data2 %>%
    filter(grepl("Ki67", Target) | grepl("SOX2", Target) | grepl("DCX", Target) | grepl("TUJ1", Target)) %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "Ki67", "SOX2", "DCX", "TUJ1", "PAX6", "TBR2"))) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y", nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")

ggsave(paste0(dir.pdfs, "D88_fold_change_eif4a2_additional.pdf"), height = 5, width = 9)

df.data %>%
    #filter(Target %in% c("ACTB", "EGFP", "EIF4A2")) %>%
    ggplot(aes(y = CT_mean, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3) +
    geom_point(position = position_dodge(width = 0.85)) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days.",
         title = "qPCR after HNP Lentivirus Transductions")

df.data %>%
    #filter(grepl("HAUS", Target) | grepl("Ki67", Target)) %>%
    ggplot(aes(y = CT_mean, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 1) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days.",
         title = "qPCR after HNP Lentivirus Transductions")

df.data2 %>%
    #filter(grepl("HAUS", Target) | grepl("Ki67", Target)) %>%
    ggplot(aes(y = CT_mean, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), lwd = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 1) +
    #geom_errorbar(aes(ymin = CT_mean - CT_sd, ymax = CT_mean + CT_sd), position = position_dodge(width = 0.85), width = 0.3, lwd = 0.3) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days.",
         title = "qPCR after HNP Lentivirus Transductions")

#ggsave(paste0(dir.pdfs, "D88_ct_mean_repeat.pdf"), height = 8, width = 8)


#dev.off()

# Presentation ####################

library(ggpubr)

df.data2 %>%
    filter(Target == "HAUS4_1" | Target == "CCND1") %>%
    bind_rows(df.data) -> df.data3

df.taqman %>%
    select(Sample, Donor, Timepoint, Expression, DOX, Replicate, Target, Name, CT_mean, CT_sd,
           fold_change_EIF4A2 = fold_change_miR361) %>%
    bind_rows(df.data3) -> df.data4

df.data4 %<>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "miR-361-5p", "miR-4707-3p-C", "EGFP", "Ki67", "CCND1", "HAUS4_1", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE))

df.data4 %<>%
    mutate(DOX = factor(DOX, levels = c(FALSE, TRUE), ordered = TRUE))

df.data4 %<>%
    filter(!is.na(fold_change_EIF4A2))

df.data4 %>%
    filter(grepl("Ki67", Target) | grepl("SOX2", Target) | grepl("DCX", Target) | grepl("TUJ1", Target) | grepl("miR-4707-3p-C", Target) | grepl("PAX6", Target) | grepl("HAUS4_1", Target) | grepl("CCND1", Target)) %>%
    #filter(Day == "Day8") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    facet_wrap(~Donor + Target, nrow = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = paperCatEight) +
    plotTheme("figure") +
    labs(y = "Fold Change",
         x = element_blank()) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D54_D88_fold_change_diff.pdf"), height = 4, width = 7)

# Scratch ###########################

# nanodrop values
#df.nano <- read_xlsx(here("results/qpcr/20210611_RNA_extractions.xlsx"), sheet = 1)

