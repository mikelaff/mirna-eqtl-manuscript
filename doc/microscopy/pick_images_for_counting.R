
library(here)
library(readr)
library(dplyr)
library(stringr)
library(magrittr)
library(ggplot2)
library(readxl)
library(mikelaffr)

# OUTPUT ###############################################################################################################

# INPUT ################################################################################################################
# D54/D88 Week 1 Diff, pTRIPZ-Ctrl or pTRIPZ-4707 or untransduced, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week1.ctrl.4707.dir <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/")
# D54/D88 Week 1 Diff, pTRIPZ-4707 stable line, +/- DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week1.4707stbl.dir <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-4707stbl_GFAP_Tuj1/")
# D54/D88 Week 2 Diff, pTRIPZ-Ctrl or pTRIPZ-4707 or untransduced, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week2.ctrl.4707.dir <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/")
# D54/D88 Week 2 Diff, pTRIPZ-4707 stable line, +/- DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week2.4707stbl.dir <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/")

# GLOBALS ##############################################################################################################

# Week 1 Ctrl / 4707 ###################################################################################################

# image output dir
output.dir <- paste0(week1.ctrl.4707.dir, "images_for_counting/")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
# image filename prefix
filename.prefix <- "Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_"
# image input dir
input.dir <- paste0(week1.ctrl.4707.dir, "images/")
# input plate layout csv
layout.csv <- paste0(week1.ctrl.4707.dir, "plate_layout.csv")
# output csv for counting
counting.csv <- paste0(week1.ctrl.4707.dir, "manual_counting.csv")
# output key
key.csv <- paste0(week1.ctrl.4707.dir, "counting_key.csv")


# import plate layout
df.layout <- read_csv(layout.csv)
df.layout %<>%
    mutate(Well = paste0(Row, str_pad(Column, 2, pad = "0")))

# get images file list
df.files <- tibble(filenames = list.files(input.dir))

# wells, rows, columns
df.files %<>%
    mutate(Well = sapply(strsplit(filenames, split = "_"), `[`, 2)) %>%
    mutate(Well = sapply(strsplit(Well, split = "Well"), `[`, 2))

df.files %<>%
    left_join(df.layout, by = c("Well"))

df.files %<>%
    mutate(Donor = factor(Donor),
           Expression = factor(Expression),
           Primary_AB = factor(Primary_AB))

# for each combination of Donor, Expression, and Primary_AB select one well at random, relable, copy image, and save key
df.files %<>%
    mutate(Label = paste(Donor, Expression, DOX, Primary_AB, sep = "_"))

# create key dataframe (select only the transduced ones: ctrl and 4707)
df.files %>%
    filter(!Expression == "Untransduced") %>%
    group_by(Label) %>%
    sample_n(8) -> df.key

# new set of filenames with random number
new.filenames <- paste0(filename.prefix, str_pad(sample(1:nrow(df.key), nrow(df.key)), 2, pad = "0"), "_", df.key$Primary_AB, ".nd2")

# assign new filenames to key
df.key$Counting_Filenames <- new.filenames

# create counting dataframe for export
df.key %>%
    ungroup() %>%
    select(Counting_Filenames) %>%
    arrange(Counting_Filenames) -> df.counting

# loop over images, copy with new filename
for (i in 1:nrow(df.key)) {
    file.copy(paste0(input.dir, df.key$filenames[i]), paste0(output.dir, df.key$Counting_Filenames[i]), overwrite = TRUE)
}

# export key and counting dataframe
write_csv(df.key, key.csv)
write_csv(df.counting, counting.csv)


# Week 2 Ctrl / 4707 ###################################################################################################

# image output dir
output.dir <- paste0(week2.ctrl.4707.dir, "images_for_counting/")
dir.create(output.dir, showWarnings = FALSE, recursive = TRUE)
# image filename prefix
filename.prefix <- "Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_"
# image input dir
input.dir <- paste0(week2.ctrl.4707.dir, "images/")
# input plate layout csv
layout.csv <- paste0(week2.ctrl.4707.dir, "plate_layout.csv")
# output csv for counting
counting.csv <- paste0(week2.ctrl.4707.dir, "manual_counting.csv")
# output key
key.csv <- paste0(week2.ctrl.4707.dir, "counting_key.csv")


# import plate layout
df.layout <- read_csv(layout.csv)
df.layout %<>%
    mutate(Well = paste0(Row, str_pad(Column, 2, pad = "0")))

# get images file list
df.files <- tibble(filenames = list.files(input.dir))

# wells, rows, columns
df.files %<>%
    mutate(Well = sapply(strsplit(filenames, split = "_"), `[`, 2)) %>%
    mutate(Well = sapply(strsplit(Well, split = "Well"), `[`, 2))

df.files %<>%
    left_join(df.layout, by = c("Well"))

df.files %<>%
    mutate(Donor = factor(Donor),
           Expression = factor(Expression),
           Primary_AB = factor(Primary_AB))

# for each combination of Donor, Expression, and Primary_AB select one well at random, relable, copy image, and save key
df.files %<>%
    mutate(Label = paste(Donor, Expression, DOX, Primary_AB, sep = "_"))

# create key dataframe
df.files %>%
    filter(!Expression == "Untransduced") %>%
    group_by(Label) %>%
    sample_n(8) -> df.key

# new set of filenames with random number
new.filenames <- paste0(filename.prefix, str_pad(sample(1:nrow(df.key), nrow(df.key)), 2, pad = "0"), "_", df.key$Primary_AB, ".nd2")

# assign new filenames to key
df.key$Counting_Filenames <- new.filenames

# create counting dataframe for export
df.key %>%
    ungroup() %>%
    select(Counting_Filenames) %>%
    arrange(Counting_Filenames) -> df.counting

# loop over images, copy with new filename
for (i in 1:nrow(df.key)) {
    file.copy(paste0(input.dir, df.key$filenames[i]), paste0(output.dir, df.key$Counting_Filenames[i]), overwrite = TRUE)
}

# export key and counting dataframe
write_csv(df.key, key.csv)
write_csv(df.counting, counting.csv)

