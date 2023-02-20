
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
# mRNA qPCR data from from 24 June 2022, RNA1-RNA24
mRNA.data1.xlsx <- here("results/qpcr/20220624_HNP_mRNA.xlsx")

# mRNA qPCR data from
#mRNA.data2.xlsx <- here("results/qpcr/20211122_HNP_mRNA-2.xlsx")

# miRNA qPCR TaqMan data from 23 June 2022, RNA1-RNA48
miRNA.taqman.data1.xlsx <- here("results/qpcr/20220623_HNP_miRNA_TaqMan.xlsx")
# miRNA qPCR TaqMan data from 23 June 2022, RNA49-RNA72
miRNA.taqman.data2.xlsx <- here("results/qpcr/20220623_HNP_miRNA_TaqMan_2.xlsx")
# miRNA qPCR TaqMan data from 8 July 2022, RNA1-RNA48
miRNA.taqman.data3.xlsx <- here("results/qpcr/20220708_HNP_miRNA_TaqMan.xlsx")

# RNA extractions, nanodrop data
sample.data.xlsx <- here("results/qpcr/20220621_RNA_extractions.xlsx")

# GLOBALS ##############################################################################################################

# Import Sample Data ###################################################################################################
df.samples <- read_xlsx(sample.data.xlsx, sheet = 1)

df.nanodrop <- read_xlsx(sample.data.xlsx, sheet = 2)

df.samples %<>%
    left_join(df.nanodrop, by = c("RNAID" = "Sample ID"))

rm(df.nanodrop)

# Import mRNA Data  Set 1 ##############################################################################################

df.data1 <- read_xlsx(mRNA.data1.xlsx, sheet = 3, range = "A36:O420", na = c("", "Undetermined"))

df.data1 %<>%
    select(RNAID = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data1 %<>%
    mutate(RNAID = paste0("RNA", sapply(strsplit(RNAID, " "), `[`, 2)))

df.data1 %<>%
    left_join(select(df.samples, RNAID, Donor, Expression, Replicate, Media, Timepoint), by = "RNAID")

df.data1$Expression <- factor(df.data1$Expression,
                             levels = c("Control", "4707"),
                             labels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                             ordered = TRUE)

df.data1$Timepoint <- factor(df.data1$Timepoint,
                            levels = c("Week 1", "Week 2"),
                            labels = c("Week 1", "Week 2"),
                            ordered = TRUE)

df.data1$Donor <- factor(df.data1$Donor)
df.data1$Replicate <- factor(df.data1$Replicate)

# df.data %<>%
#     mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))


# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data1 %<>%
    group_by(RNAID, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# filter out high standard dev samples
df.data1 %<>%
    filter(CT_sd < 0.8) %>%
    filter(!(RNAID == "RNA15"))

# loop over samples, calculate delta CT to ACTB, EIF4A2, GAPDH
samples <- unique(df.data1$RNAID)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data1 %>%
        filter(RNAID == sample) -> df.tmp

    # get CT value for this sample
    ct.actb <- df.tmp$CT_mean[match("ACTB", df.tmp$Target)]
    ct.eif <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]
    ct.gapdh <- df.tmp$CT_mean[match("GAPDH", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_ACTB = CT_mean - ct.actb,
               delta_CT_EIF4A2 = CT_mean - ct.eif,
               delta_CT_GAPDH = CT_mean - ct.gapdh)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)
}

df.data1 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Control of that donor for each target
donors <- unique(df.data1$Donor)
timepoints <- unique(df.data1$Timepoint)
targets <- unique(df.data1$Target)

df.new <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data1 %>%
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
            delta.ct.control.actb <- mean(df.target$delta_CT_ACTB[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
            delta.ct.control.eif4a2 <- mean(df.target$delta_CT_EIF4A2[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)
            delta.ct.control.gapdh <- mean(df.target$delta_CT_GAPDH[which(df.target$Expression == "pTRIPZ-Control")], na.rm = TRUE)

            # calculate delta delta CT within this target and donor
            df.target %<>%
                mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control.actb,
                       delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control.eif4a2,
                       delta_delta_CT_GAPDH = delta_CT_GAPDH - delta.ct.control.gapdh)

            # calculate fold change
            df.target %<>%
                mutate(fold_change_ACTB = 2 ^ (-delta_delta_CT_ACTB),
                       fold_change_EIF4A2 = 2 ^ (-delta_delta_CT_EIF4A2),
                       fold_change_GAPDH = 2 ^ (-delta_delta_CT_GAPDH))

            # combine into all data
            df.new <- bind_rows(df.new, df.target)
        }
    }
}

df.data1 <- df.new
rm(df.donor, df.target, df.day, df.new)

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
df.taqman1 <- read_xlsx(miRNA.taqman.data1.xlsx, sheet = 3, range = "A35:O419", na = c("", "Undetermined"))

df.taqman1 %<>%
    select(RNAID = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

# remove space from RNAID
df.taqman1 %<>%
    mutate(RNAID = paste0("RNA", sapply(strsplit(RNAID, " "), `[`, 2)))

df.taqman2 <- read_xlsx(miRNA.taqman.data2.xlsx, sheet = 3, range = "A35:O227", na = c("", "Undetermined"))

df.taqman2 %<>%
    select(RNAID = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman2 %<>%
    mutate(RNAID = paste0("RNA", sapply(strsplit(RNAID, " "), `[`, 2)))


df.taqman <- bind_rows(df.taqman1, df.taqman2)

rm(df.taqman1, df.taqman2)

df.taqman %<>%
    left_join(select(df.samples, RNAID, Donor, Expression, Replicate, Media, Timepoint), by = "RNAID")


df.taqman$Expression <- factor(df.taqman$Expression,
                               levels = c("Control", "4707"),
                               labels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                               ordered = TRUE)


df.taqman$Timepoint <- factor(df.taqman$Timepoint,
                              levels = c("Week 1", "Week 2", "Day 8"),
                              labels = c("Week 1", "Week 2", "Day 8"),
                              ordered = TRUE)

df.taqman$Donor <- factor(df.taqman$Donor)
df.taqman$Replicate <- factor(df.taqman$Replicate)

# df.taqman %<>%
#     mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.taqman %<>%
    group_by(RNAID, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

# loop over samples, calculate delta CT to miR-361
samples <- unique(df.taqman$RNAID)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.taqman %>%
        filter(RNAID == sample) -> df.tmp

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

# Diff Experiment
# select diff samples, loop over donors, calculate delta delta CT to Control of that donor for each target
df.taqman %>%
    filter(Media == "Diff") -> df.diff

donors <- unique(df.diff$Donor)
timepoints <- unique(df.diff$Timepoint)
targets <- unique(df.diff$Target)

df.new <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.diff %>%
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
            delta.ct.control <- mean(df.target$delta_CT_miR361[which(df.target$Expression == "pTRIPZ-Control")])

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
}

df.diff <- df.new
rm(df.new, df.donor, df.day, df.target)

# Prolif Experiment
# select diff samples, loop over donors, calculate delta delta CT to Control of that donor for each target
df.taqman %>%
    filter(Media == "Prolif") -> df.prolif

donors <- unique(df.prolif$Donor)
timepoints <- unique(df.prolif$Timepoint)
targets <- unique(df.prolif$Target)

df.new <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.prolif %>%
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
            delta.ct.control <- mean(df.target$delta_CT_miR361[which(df.target$Expression == "pTRIPZ-Control")])

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
}

df.prolif <- df.new
rm(df.new, df.donor, df.day, df.target)




# Import TaqMan miRNA Data, More Controls ##############################################################################
df.taqman3 <- read_xlsx(miRNA.taqman.data3.xlsx, sheet = 3, range = "A35:O419", na = c("", "Undetermined"))

df.taqman3 %<>%
    select(RNAID = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

# remove space from RNAID
df.taqman3 %<>%
    mutate(RNAID = paste0("RNA", sapply(strsplit(RNAID, " "), `[`, 2)))

df.taqman3 %<>%
    left_join(select(df.samples, RNAID, Donor, Expression, Replicate, Media, Timepoint), by = "RNAID")


df.taqman3$Expression <- factor(df.taqman3$Expression,
                               levels = c("Control", "4707"),
                               labels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                               ordered = TRUE)


df.taqman3$Timepoint <- factor(df.taqman3$Timepoint,
                              levels = c("Week 1", "Week 2", "Day 8"),
                              labels = c("Week 1", "Week 2", "Day 8"),
                              ordered = TRUE)

df.taqman3$Donor <- factor(df.taqman3$Donor)
df.taqman3$Replicate <- factor(df.taqman3$Replicate)

# df.taqman %<>%
#     mutate(Name = paste(Donor, Expression, DOX, Timepoint, Replicate, sep = "_"))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.taqman3 %<>%
    group_by(RNAID, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

df.taqman3 %<>%
    filter(CT_sd < 1)

df.taqman3 %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = .7, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = .7) +
    facet_wrap(~Donor + Target, scales = "free_y", nrow = 2) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(label.x=2, size = 5, method = "t.test", label = "p.format")



# loop over samples, calculate delta CT to miR-361
samples <- unique(df.taqman3$RNAID)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.taqman3 %>%
        filter(RNAID == sample) -> df.tmp

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

df.taqman3 <- df.new
rm(df.new, df.tmp)

# Diff Experiment
# select diff samples, loop over donors, calculate delta delta CT to Control of that donor for each target
df.taqman3 %>%
    filter(Media == "Diff") -> df.diff

donors <- unique(df.diff$Donor)
timepoints <- unique(df.diff$Timepoint)
targets <- unique(df.diff$Target)

df.new <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.diff %>%
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
    }
}

df.diff <- df.new
rm(df.new, df.donor, df.day, df.target)

# Prolif Experiment
# select diff samples, loop over donors, calculate delta delta CT to Control of that donor for each target
# df.taqman %>%
#     filter(Media == "Prolif") -> df.prolif
#
# donors <- unique(df.prolif$Donor)
# timepoints <- unique(df.prolif$Timepoint)
# targets <- unique(df.prolif$Target)
#
# df.new <- tibble()
#
# for (donor in donors) {
#     df.donor <- NULL
#
#     # filter for only this donor
#     df.prolif %>%
#         filter(Donor == donor) -> df.donor
#
#     for (timepoint in timepoints) {
#         df.day <- NULL
#
#         # filter for only this day
#         df.donor %>%
#             filter(Timepoint == timepoint) -> df.day
#
#         # loop over targets
#         for (target in targets) {
#             df.target <- NULL
#
#             # filter for only this target
#             df.day %>%
#                 filter(Target == target) -> df.target
#
#             # use the mean control value across the replicates for delta delta ct
#             delta.ct.control <- mean(df.target$delta_CT_miR361[which(df.target$Expression == "pTRIPZ-Control")])
#
#             # calculate delta delta CT within this target and donor
#             df.target %<>%
#                 mutate(delta_delta_CT_miR361 = delta_CT_miR361 - delta.ct.control)
#
#             # calculate fold change
#             df.target %<>%
#                 mutate(fold_change_miR361 = 2 ^ (-delta_delta_CT_miR361))
#
#             # combine into all data
#             df.new <- bind_rows(df.new, df.target)
#         }
#     }
# }
#
# df.prolif <- df.new
# rm(df.new, df.donor, df.day, df.target)




# Plot #################################################################################################################

# my.theme <- theme(axis.title = element_text(size = 16),
#                   axis.text = element_text(size = 14),
#                   axis.text.x = element_text(size = 12),
#                   title = element_text(size = 18),
#                   legend.title = element_text(size = 16),
#                   legend.text = element_text(size = 14),
#                   strip.text = element_text(size = 16),
#                   plot.caption = element_text(size = 12))

#pdf(paste0(dir.pdfs, "20211123_Diff_1_2_week_D54_D88_qPCR_foldChange.pdf"), height = 8, width = 8)

df.diff %>%
    filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR361, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Donor + Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme("figure") +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(label.x=2, label.y=.25, size = 2, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "20220624_Diff_1-2week_miRNA_fold_change_withSig.pdf"), height = 3.5, width = 5)

df.diff %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    facet_wrap(~Donor + Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme("figure") +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Week",
         title = "Differentiation Experiment")

ggsave(paste0(dir.pdfs, "20220624_Diff_1-2week_miRNA_meanCT.pdf"), height = 4, width = 5)

df.prolif %>%
    filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR361, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme("figure") +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Donor",
         title = "Proliferation Experiment") +
    stat_compare_means(label.x=2, label.y=.25, size = 2, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "20220624_Prolif_miRNA_fold_change_withSig.pdf"), height = 3.5, width = 3.5)

df.prolif %>%
    ggplot(aes(y = CT_mean, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 0.3) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme("figure") +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Proliferation Experiment")

ggsave(paste0(dir.pdfs, "20220624_Prolif_miRNA_meanCT.pdf"), height = 2, width = 5)

df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `260/230`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.3), size = 4) +
    facet_wrap(~Donor + Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Week",
         title = "Differentiation Experiment")

#ggsave(paste0(dir.pdfs, "20220624_Diff_1-2week_miRNA_meanCT.pdf"), height = 4, width = 5)

pdf(paste0(dir.pdfs, "20220628_Diff_1_week_D54_D88_qPCR.pdf"), height = 8, width = 10)

df.data1 %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    ggplot(aes(y = fold_change_ACTB, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, Endog.Control: ACTB") +
    stat_compare_means(label.x=2, label.y=.25, method = "t.test", label = "p.format")

df.data1 %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, Endog.Control: EIF4A2") +
    stat_compare_means(label.x=2, label.y=.25, method = "t.test", label = "p.format")

df.data1 %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    ggplot(aes(y = fold_change_GAPDH, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85)) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, Endog.Control: GAPDH") +
    stat_compare_means(label.x=2, label.y=.25, method = "t.test", label = "p.format")

df.data1 %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    ggplot(aes(y = CT_mean, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Differentiation Experiment, Week 1") +
    stat_compare_means(method = "t.test", label = "p.format")

dev.off()


pdf(paste0(dir.pdfs, "20220628_Diff_1_week_D54_D88_raw_CT_rna_quality.pdf"), height = 8, width = 10)

df.data1 %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    filter(Target %in% c("ACTB", "EIF4A2", "GAPDH")) %>%
    ggplot(aes(y = CT_mean, x = Donor, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `ng/ul`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.4), size = 4) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, RNA conc.") +
    stat_compare_means(method = "t.test", label = "p.format")

df.data1 %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    filter(Target %in% c("ACTB", "EIF4A2", "GAPDH")) %>%
    ggplot(aes(y = CT_mean, x = Donor, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `260/280`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.4), size = 4) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, 260/280") +
    stat_compare_means(method = "t.test", label = "p.format")

df.data1 %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "GAPDH", "KI67", "PAX6", "SOX2", "DCX", "TUJ1"), ordered = TRUE)) %>%
    filter(Target %in% c("ACTB", "EIF4A2", "GAPDH")) %>%
    ggplot(aes(y = CT_mean, x = Donor, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `260/230`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.4), size = 4) +
    facet_wrap(~Target, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Differentiation Experiment, Week 1, 260/230") +
    stat_compare_means(method = "t.test", label = "p.format")

dev.off()

pdf(paste0(dir.pdfs, "20220628_Diff_1_week_ng_ul.pdf"), height = 5, width = 6)

df.samples %>%
    filter(Timepoint == "Week 1") %>%
    ggplot(aes(y = `ng/ul`, x = Donor, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `ng/ul`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.3), size = 4) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    stat_compare_means(method = "t.test", label = "p.format")


df.samples %>%
    filter(Timepoint == "Week 1") %>%
    ggplot(aes(y = `ng/ul`, x = Expression, shape = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `ng/ul`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.3), size = 4) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    stat_compare_means(method = "t.test", label = "p.format")

dev.off()

df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, shape = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), outlier.shape = NA) +
    geom_point(mapping = aes(color = `ng/ul`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.4), size = 4) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Donor",
         title = "Differentiation Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = `ng/ul`, shape = Expression)) +
    geom_point(mapping = aes(color = `ng/ul`), position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.4), size = 4) +
    facet_wrap(~Target + Donor, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    #scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "ng/ul",
         title = "Differentiation Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

pdf(paste0(dir.pdfs, "20220624_Diff_miRNA_meanCT_by_RNA.pdf"), height = 7, width = 9)


df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = `ng/ul`)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    facet_wrap(~Target + Donor, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "ng/ul",
         title = "Differentiation Experiment")

df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = `260/280`)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    facet_wrap(~Target + Donor, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "260/280",
         title = "Differentiation Experiment")

df.diff %>%
    left_join(select(df.samples, RNAID, `ng/ul`, A260, A280, `260/280`, `260/230`), by = "RNAID") %>%
    ggplot(aes(y = CT_mean, x = `260/230`)) +
    geom_point(size = 2) +
    geom_smooth(method = "lm") +
    facet_wrap(~Target + Donor, scales = "free") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    #scale_color_gradient(low = "blue", high = "red") +
    #theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "260/230",
         title = "Differentiation Experiment")


dev.off()
#ggsave(paste0(dir.pdfs, "20220624_Diff_miRNA_meanCT_by_conc.pdf"), height = 7, width = 9)


pdf(paste0(dir.pdfs, "20220719_Diff_miRNA_controls.pdf"), height = 7, width = 12)

df.diff %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR197, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Donor + Target, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(label.x=2, label.y=.25, size = 2, method = "t.test", label = "p.format")

df.diff %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR324, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Donor + Target, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(label.x=2, label.y=.25, size = 2, method = "t.test", label = "p.format")

df.diff %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = fold_change_miR361, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    facet_wrap(~Donor + Target, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "Fold Change",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(label.x=2, label.y=.25, size = 2, method = "t.test", label = "p.format")

df.diff %>%
    #filter(Target == "miR-4707-3p") %>%
    ggplot(aes(y = CT_mean, x = Timepoint, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_dodge(width = 0.85), size = 2) +
    facet_wrap(~Donor + Target, nrow = 2, scales = "free_y") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    plotTheme() +
    theme(strip.text.x = element_text(size = 7, margin = margin(t = 1.5, r = 0, b = 1.5, l = 0, unit = "pt"))) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    labs(y = "mean(CT)",
         x = "Week",
         title = "Differentiation Experiment") +
    stat_compare_means(method = "t.test", label = "p.format")

dev.off()
