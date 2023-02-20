
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/paper/figure6/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# mRNA qPCR data from 13 June 2021 (prolif)
prolif.mRNA.data1.xlsx <- here("results/qpcr/20210702_HNP_mRNA_repeat.xlsx")
# miRNA qPCR TaqMan data from 13 June 2021 (prolif)
prolif.miRNA.taqman.data.xlsx <- here("results/qpcr/20210613_HNP_miRNA_TaqMan.xlsx")
# additional mRNA qPCR data from 30 July 2021 (prolif)
prolif.mRNA.data2.xlsx <- here("results/qpcr/20210730_HNP_mRNA.xlsx")

# mRNA qPCR data from 22 Nov 2021, set 1 (diff)
diff.mRNA.data1.xlsx <- here("results/qpcr/20211122_HNP_mRNA-1.xlsx")
# mRNA qPCR data from 22 Nov 2021, set 2 (diff)
diff.mRNA.data2.xlsx <- here("results/qpcr/20211122_HNP_mRNA-2.xlsx")
# miRNA qPCR TaqMan data from 19 Nov 2021 (diff)
diff.miRNA.taqman.data.xlsx <- here("results/qpcr/20211119_HNP_miRNA_TaqMan.xlsx")


# GLOBALS ##############################################################################################################

# Import mRNA Data 1 Prolif ############################################################################################

df.data <- read_xlsx(prolif.mRNA.data1.xlsx, sheet = 3, range = "A46:O334", na = c("", "Undetermined"))

df.data %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Day = sapply(strsplit(Sample, "_"), `[`, 2),
           Expression = sapply(strsplit(Sample, "_"), `[`, 3),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 4))

df.data %<>%
    select(Sample,
           Donor,
           Day,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.data$Expression <- factor(df.data$Expression,
                             levels = c("Control", "4707"),
                             labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                             ordered = TRUE)

df.data$Day <- factor(df.data$Day)
df.data$Donor <- factor(df.data$Donor)
df.data$Replicate <- factor(df.data$Replicate)

df.data %<>%
    mutate(Name = paste(Donor, Expression, Day, Replicate))

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
days <- unique(df.data$Day)
targets <- unique(df.data$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.donor

    for (day in days) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Day == day) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
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
}

df.data <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

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
days <- unique(df.data$Day)
targets <- unique(df.data$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.donor

    for (day in days) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Day == day) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
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
}

df.data <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

df.prolif.mRNA1 <- df.data
rm(df.data)

# Import mRNA Data 2 Prolif ############################################################################################

df.data2 <- read_xlsx(prolif.mRNA.data2.xlsx, sheet = 3, range = "A46:O334", na = c("", "Undetermined"))

df.data2 %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data2 %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Day = sapply(strsplit(Sample, "_"), `[`, 2),
           Expression = sapply(strsplit(Sample, "_"), `[`, 3),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 4))

df.data2 %<>%
    select(Sample,
           Donor,
           Day,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.data2$Expression <- factor(df.data2$Expression,
                              levels = c("Control", "4707"),
                              labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                              ordered = TRUE)

df.data2$Day <- factor(df.data2$Day)
df.data2$Donor <- factor(df.data2$Donor)
df.data2$Replicate <- factor(df.data2$Replicate)

df.data2 %<>%
    mutate(Name = paste(Donor, Expression, Day, Replicate))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data2 %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
           CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()

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

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data2 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to pTRIPZ-Control of that donor for each target
donors <- unique(df.data2$Donor)
days <- unique(df.data2$Day)
targets <- unique(df.data2$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data2 %>%
        filter(Donor == donor) -> df.donor

    for (day in days) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Day == day) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
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
}

df.data2 <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

# Repeat for EIF4A2

# loop over samples, calculate delta CT to EIF4A2
samples <- unique(df.data2$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data2 %>%
        filter(Sample == sample) -> df.tmp

    # get ACTB CT value for this sample
    ct.eif4a2 <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_EIF4A2 = CT_mean - ct.eif4a2)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data2 <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to pTRIPZ-Control of that donor for each target
donors <- unique(df.data2$Donor)
days <- unique(df.data2$Day)
targets <- unique(df.data2$Target)

df.tripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.data2 %>%
        filter(Donor == donor) -> df.donor

    for (day in days) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Day == day) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
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
}

df.data2 <- df.tripz
rm(df.donor, df.target, df.day, df.tripz)

df.prolif.mRNA2 <- df.data2
rm(df.data2)

# Import TaqMan miRNA Data Prolif ######################################################################################
df.taqman <- read_xlsx(prolif.miRNA.taqman.data.xlsx, sheet = 3, range = "A45:O117", na = c("", "Undetermined"))

df.taqman %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.taqman %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Day = sapply(strsplit(Sample, "_"), `[`, 2),
           Expression = sapply(strsplit(Sample, "_"), `[`, 3),
           Replicate = sapply(strsplit(Sample, "_"), `[`, 4))

df.taqman %<>%
    select(Sample,
           Donor,
           Day,
           Expression,
           Replicate,
           Well,
           Target,
           CT)

df.taqman$Expression <- factor(df.taqman$Expression,
                               levels = c("Control", "4707"),
                               labels = c("pTRIPZ-Control", "pTRIPZ-4707-C"),
                               ordered = TRUE)

df.taqman$Day <- factor(df.taqman$Day)
df.taqman$Donor <- factor(df.taqman$Donor)
df.taqman$Replicate <- factor(df.taqman$Replicate)

df.taqman %<>%
    mutate(Name = paste(Donor, Expression, Day, Replicate))

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
days <- unique(df.taqman$Day)
targets <- unique(df.taqman$Target)

df.ptripz <- tibble()

for (donor in donors) {
    df.donor <- NULL

    # filter for only this donor
    df.taqman %>%
        filter(Donor == donor) -> df.donor

    for (day in days) {
        df.day <- NULL

        # filter for only this day
        df.donor %>%
            filter(Day == day) -> df.day

        # loop over targets
        for (target in targets) {
            df.target <- NULL

            # filter for only this target
            df.day %>%
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
}

df.taqman <- df.ptripz
rm(df.donor, df.day, df.target, df.ptripz)

df.prolif.taqman <- df.taqman
rm(df.taqman)

# Import mRNA Data 1 Diff ##############################################################################################

df.data <- read_xlsx(diff.mRNA.data1.xlsx, sheet = 3, range = "A46:O430", na = c("", "Undetermined"))

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

df.diff.mRNA1 <- df.data
rm(df.data)

# Import mRNA Data 2 Diff ##############################################################################################

df.data2 <- read_xlsx(diff.mRNA.data2.xlsx, sheet = 3, range = "A46:O430", na = c("", "Undetermined"))

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

df.diff.mRNA2 <- df.data2
rm(df.data2)

# Import TaqMan miRNA Data Diff ########################################################################################
df.taqman <- read_xlsx(diff.miRNA.taqman.data.xlsx, sheet = 3, range = "A45:O141", na = c("", "Undetermined"))

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

df.diff.taqman <- df.taqman
rm(df.taqman)

# Plots ################################################################################################################

df.prolif.mRNA1 %>%
    filter(Target == "HAUS4_1" | Target == "CCND1") %>%
    bind_rows(df.prolif.mRNA2) -> df.prolif

df.prolif.taqman %>%
    select(Sample, Donor, Day, Expression, Replicate, Target, CT_mean, CT_sd,
           fold_change_EIF4A2 = fold_change_miR361) %>%
    bind_rows(df.prolif) -> df.prolif

df.prolif %<>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "miR-361-5p", "miR-4707-3p-C", "Ki67", "CCND1", "HAUS4_1", "PAX6", "SOX2", "DCX", "TBR2", "TUJ1")))

df.prolif %>%
    filter(grepl("Ki67", Target) | grepl("SOX2", Target) | grepl("DCX", Target) | grepl("TUJ1", Target) | grepl("miR-4707-3p-C", Target) | grepl("PAX6", Target) | grepl("HAUS4_1", Target) | grepl("CCND1", Target)) %>%
    #filter(Day == "Day8") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Day, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
    #geom_point(position = position_dodge(width = 0.85), size = 0.2) +
    facet_wrap(~Target, nrow = 1) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")
    #stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D88_prolif_fold_change_figure.pdf"), height = 2, width = 6.5)


df.diff.mRNA2 %>%
    filter(Target == "HAUS4_1" | Target == "CCND1") %>%
    bind_rows(df.diff.mRNA1) -> df.diff

df.diff.taqman %>%
    select(Sample, Donor, Timepoint, Expression, DOX, Replicate, Target, CT_mean, CT_sd,
           fold_change_EIF4A2 = fold_change_miR361) %>%
    bind_rows(df.diff) -> df.diff

df.diff %<>%
    mutate(Target = factor(Target, levels = c("ACTB", "EIF4A2", "miR-361-5p", "miR-4707-3p-C", "Ki67", "CCND1", "HAUS4_1", "PAX6", "SOX2", "DCX", "TBR2", "TUJ1")))

df.diff %<>%
    mutate(DOX = factor(DOX, levels = c(FALSE, TRUE), ordered = TRUE))

df.diff %>%
    filter(grepl("Ki67", Target) | grepl("SOX2", Target) | grepl("DCX", Target) | grepl("TUJ1", Target) | grepl("miR-4707-3p-C", Target) | grepl("PAX6", Target) | grepl("HAUS4_1", Target) | grepl("CCND1", Target)) %>%
    #filter(Day == "Day8") %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Timepoint, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
    #geom_point(position = position_dodge(width = 0.85), size = 0.2) +
    facet_wrap(~Donor + Target, nrow = 2) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Day",
         caption = "Donor 88. HNP expression after 4-8 days. Endogenous Control: EIF4A2",
         title = "qPCR after HNP Lentivirus Transductions")
#stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D54_D88_diff_fold_change_figure.pdf"), height = 4, width = 6.5)

