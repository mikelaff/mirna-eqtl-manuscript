
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)
library(ggpubr)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/microscopy/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
# D54/D88 Week 1 Diff, pTRIPZ-ctrl or pTRIPZ-4707 low MOI, +DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
manual.counts.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/manual_counting.csv")
image.key.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/counting_key.csv")

# D54/D88 Week 2 Diff, pTRIPZ-ctrl or pTRIPZ-4707 low MOI, +DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)


# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################

df.key <- read_csv(image.key.csv)

df.counts <- read_csv(manual.counts.csv)

df.counts %<>%
    left_join(df.key, by = "Counting_Filenames")

df.counts %>%
    filter(!is.na(image_counted)) %>%
    mutate(fraction_tuj1 = gfp_cy5_cells / gfp_cells) %>%
    ggplot(aes(x = Expression, y = fraction_tuj1)) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(width = 0.1) +
    facet_wrap(~Donor)




    #stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")









