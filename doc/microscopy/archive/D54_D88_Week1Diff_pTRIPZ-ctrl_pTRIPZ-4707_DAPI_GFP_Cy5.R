
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/microscopy/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# D54/D88 Week 1 Diff, pTRIPZ-Ctrl/4707, DAPI/GFP/Cy5 (Tuj1 or GFAP): Nuclei and Dilated Nuclei
nuclei.data.csv <- here("results/microscopy/output/D54_D88_Week1Diff_DAPI_GFP_Cy5_Nuclei.csv")
dilated.data.csv <- here("results/microscopy/output/D54_D88_Week1Diff_DAPI_GFP_Cy5_DilatedNuclei.csv")

# GLOBALS ##############################################################################################################


# Import Data #####################################################################################################

# nuclei measurements
df.nuclei <- read_csv(nuclei.data.csv)

df.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           File = Metadata_FileLocation,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity"))

df.nuclei$Donor <- NA
df.nuclei$Donor[df.nuclei$Row %in% c("C", "D")] <- "D54"
df.nuclei$Donor[df.nuclei$Row %in% c("E", "F")] <- "D88"

df.nuclei$Column <- as.numeric(df.nuclei$Column)
df.nuclei$Cy5_Ab <- NA
df.nuclei$Cy5_Ab[df.nuclei$Column %in% c(1,2,3,7,8,9)] <- "Tuj1"
df.nuclei$Cy5_Ab[df.nuclei$Column %in% c(4,5,6,10,11,12)] <- "GFAP"

df.nuclei$Expression <- NA
df.nuclei$Expression[df.nuclei$Column %in% c(1,6,7,12)] <- "Untransduced"
df.nuclei$Expression[df.nuclei$Column %in% c(2,3,4,5)] <- "Control"
df.nuclei$Expression[df.nuclei$Column %in% c(8,9,10,11)] <- "4707"

df.nuclei %>%
    #filter(Expression == "Control") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_GFP)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Donor + Expression, scales = "free")

df.nuclei %>%
    #filter(Expression == "Control") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_Cy5)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Cy5_Ab + Donor + Expression, scales = "free", nrow = 2)

df.nuclei %>%
    filter(Cy5_Ab == "Tuj1") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_Cy5)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Donor + Expression, scales = "fixed", nrow = 2)

df.nuclei %>%
    filter(Cy5_Ab == "GFAP") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_Cy5)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Donor + Expression, scales = "fixed", nrow = 2)

df.nuclei %<>%
    mutate(GFP_positive = Intensity_IntegratedIntensity_GFP > 75,
           Cy5_positive = Intensity_IntegratedIntensity_Cy5 > 100)

df.nuclei %>%
    group_by(ImageNumber, Donor, Expression, Cy5_Ab) %>%
    summarise(num_nuclei = n(),
              GFP_total_intensity = sum(Intensity_IntegratedIntensity_GFP),
              Cy5_total_intensity = sum(Intensity_IntegratedIntensity_Cy5),
              GFP_pos_nuclei = sum(GFP_positive),
              Cy5_pos_nuclei = sum(Cy5_positive),
              GFP_Cy5_pos_nuclei = sum(GFP_positive & Cy5_positive)) %>%
    mutate(GFP_over_nuclei = GFP_total_intensity / num_nuclei,
           Cy5_over_nuclei = Cy5_total_intensity / num_nuclei) %>%
    mutate(Cy5_over_GFP = Cy5_over_nuclei / GFP_over_nuclei,
           GFP_pos_frac = GFP_pos_nuclei / num_nuclei,
           Cy5_pos_frac = Cy5_pos_nuclei / num_nuclei,
           GFP_Cy5_pos_frac = GFP_Cy5_pos_nuclei / num_nuclei) -> df.summary

df.summary %>%
    ggplot(mapping = aes(x = Expression, y = GFP_pos_frac, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.summary %>%
    ggplot(mapping = aes(x = Expression, y = Cy5_pos_frac, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    facet_wrap(~Cy5_Ab)

df.summary %>%
    ggplot(mapping = aes(x = Expression, y = GFP_Cy5_pos_frac, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    facet_wrap(~Cy5_Ab)

df.summary %>%
    ggplot(mapping = aes(x = Expression, y = GFP_over_nuclei, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.summary %>%
    ggplot(mapping = aes(x = Expression, y = Cy5_over_nuclei, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Cy5_Ab) +
    scale_color_manual(values = c("darkblue", "darkorange"))




# dilated nuclei measurements
df.dilated <- read_csv(dilated.data.csv)

df.dilated %<>%
    select(ImageNumber,
           ObjectNumber,
           File = Metadata_FileLocation,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity"))

df.dilated$Donor <- NA
df.dilated$Donor[df.dilated$Row %in% c("C", "D")] <- "D54"
df.dilated$Donor[df.dilated$Row %in% c("E", "F")] <- "D88"

df.dilated$Column <- as.numeric(df.dilated$Column)
df.dilated$Cy5_Ab <- NA
df.dilated$Cy5_Ab[df.dilated$Column %in% c(1,2,3,7,8,9)] <- "Tuj1"
df.dilated$Cy5_Ab[df.dilated$Column %in% c(4,5,6,10,11,12)] <- "GFAP"

df.dilated$Expression <- NA
df.dilated$Expression[df.dilated$Column %in% c(1,6,7,12)] <- "Untransduced"
df.dilated$Expression[df.dilated$Column %in% c(2,3,4,5)] <- "Control"
df.dilated$Expression[df.dilated$Column %in% c(8,9,10,11)] <- "4707"

df.dilated %>%
    #filter(Expression == "Control") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_GFP)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Donor + Expression, scales = "free")


df.dilated %>%
    group_by(ImageNumber, Donor, Expression, Cy5_Ab) %>%
    summarise(num_nuclei = n(),
              GFP_total_intensity = sum(Intensity_IntegratedIntensity_GFP),
              Cy5_total_intensity = sum(Intensity_IntegratedIntensity_Cy5)) %>%
    mutate(GFP_over_nuclei = GFP_total_intensity / num_nuclei,
           Cy5_over_nuclei = Cy5_total_intensity / num_nuclei) %>%
    mutate(Cy5_over_GFP = Cy5_over_nuclei / GFP_over_nuclei) -> df.dilated.summary


df.dilated.summary %>%
    ggplot(mapping = aes(x = Expression, y = GFP_over_nuclei, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    scale_color_manual(values = c("darkblue", "darkorange"))

df.dilated.summary %>%
    ggplot(mapping = aes(x = Expression, y = Cy5_over_GFP, color = Donor)) +
    geom_boxplot(position = position_dodge(width = 0.85)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Cy5_Ab) +
    scale_color_manual(values = c("darkblue", "darkorange"))


























df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_DAPI, color = Expression)) +
    geom_density()

df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFPcorrected)) +
    geom_histogram()

df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_Tuj1corrected)) +
    geom_histogram()

df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_EdUcorrected)) +
    geom_histogram()

# classify nuclei
df.nuclei %<>%
    mutate(GFP_positive = Intensity_IntegratedIntensity_GFPcorrected > 1,
           Tuj1_positive = Intensity_IntegratedIntensity_Tuj1corrected > 1,
           EdU_positive = Intensity_IntegratedIntensity_EdUcorrected > 1)

df.nuclei %>%
    group_by(ImageNumber) %>%
    summarise(num_nuclei = n(),
              GFP_positive_nuclei = sum(GFP_positive),
              Tuj1_positive_nuclei = sum(Tuj1_positive),
              EdU_positive_nuclei = sum(EdU_positive),
              GFP_Tuj1_positive_nuclei = sum(GFP_positive & Tuj1_positive),
              GFP_EdU_positive_nuclei = sum(GFP_positive & EdU_positive)) %>%
    mutate(GFP_positive_frac = GFP_positive_nuclei / num_nuclei,
           Tuj1_positive_frac = Tuj1_positive_nuclei / num_nuclei,
           EdU_positive_frac = EdU_positive_nuclei / num_nuclei,
           GFP_Tuj1_positive_frac = GFP_Tuj1_positive_nuclei / num_nuclei,
           GFP_EdU_positive_frac = GFP_EdU_positive_nuclei / num_nuclei) %>%
    left_join(df.images, by = "ImageNumber") -> df.summary


df.summary %>%
    ggplot(aes(x = Expression, y = GFP_positive_frac)) +
    geom_boxplot()

df.summary %>%
    ggplot(aes(x = Expression, y = Tuj1_positive_frac)) +
    geom_boxplot()

df.summary %>%
    ggplot(aes(x = Expression, y = EdU_positive_frac)) +
    geom_boxplot()

df.summary %>%
    ggplot(aes(x = Expression, y = GFP_Tuj1_positive_frac)) +
    geom_boxplot()

df.summary %>%
    ggplot(aes(x = Expression, y = GFP_EdU_positive_frac)) +
    geom_boxplot()




