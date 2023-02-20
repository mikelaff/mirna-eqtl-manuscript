
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
# Day 8 Prolif experiment with Donor54, stained DAPI, EdU, GFP, and Tuj1
data.csv <- here("results/microscopy/D54_Day8_Prolif_10xby4images/output/Day8Prolif_D54_pTRIPZ_Image.csv")

# Day 8 Prolif experiment with Donor54, stained DAPI, EdU, GFP, and Tuj1, NUCLEI Measurements
nuclei.csv <- here("results/microscopy/D54_Day8_Prolif_10xby4images_nd2/output/Day8Prolif_D54_pTRIPZ_Nuclei.csv")
images.csv <- here("results/microscopy/D54_Day8_Prolif_10xby4images_nd2/output/Day8Prolif_D54_pTRIPZ_Image.csv")

# GLOBALS ##############################################################################################################


# Import Data #####################################################################################################

# df.data <- read_csv(data.csv)
#
# df.data %>%
#     ggplot(aes(x = Expression, y = Classify_Tuj1_positive_PctObjectsPerBin)) +
#     geom_boxplot()
#
# df.data %>%
#     ggplot(aes(x = Expression, y = Classify_EdU_positive_PctObjectsPerBin)) +
#     geom_boxplot()
#
# df.data %>%
#     ggplot(aes(x = Expression, y = Classify_GFP_pos_PctObjectsPerBin)) +
#     geom_boxplot()
#
# df.data %>%
#     ggplot(aes(x = Expression, y = Count_Nuclei)) +
#     geom_boxplot()

# nuclei measurements
df.nuclei <- read_csv(nuclei.csv)

# image metadata
df.images <- read_csv(images.csv)
df.images %<>%
    select(ImageNumber, Expression)

df.nuclei %<>%
    left_join(df.images, by = "ImageNumber")

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




