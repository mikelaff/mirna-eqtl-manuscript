
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
# D54 Week 1 Diff, pTRIPZ-Control or pTRIPZ-4707, +DOX, DAPI/GFP/TxRed(NES)/Cy7(Tuj1)
d54.week1.nuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_redoGFP_Nuclei.csv")
d54.week1.Cy7cells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_redoGFP_Cy7_cells.csv")
d54.week1.TxRedcells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_redoGFP_TxRed_cells.csv")
d54.week1.GFPcells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_redoGFP_GFP_cells.csv")
d54.week1.images.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_redoGFP_Image.csv")
d54.week1.plate.layout.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/D54_week1_plate_layout.csv")

# D54 Week 2 Diff, pTRIPZ-Control or pTRIPZ-4707, +DOX, DAPI/GFP/TxRed(NES)/Cy7(Tuj1)
d54.week2.nuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_redoGFP_Nuclei.csv")
d54.week2.Cy7cells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_redoGFP_Cy7_cells.csv")
d54.week2.TxRedcells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_redoGFP_TxRed_cells.csv")
d54.week2.GFPcells.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_redoGFP_GFP_cells.csv")
d54.week2.images.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_redoGFP_Image.csv")
d54.week2.plate.layout.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/D54_week2_plate_layout.csv")

# Donut Analysis: D54 Week 1 Diff, pTRIPZ-Control or pTRIPZ-4707, +DOX, DAPI/GFP/TxRed(NES)/Cy7(Tuj1)
#d54.week1.donut.nuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_donuts_Nuclei.csv")
#d54.week1.donut.expandednuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_donuts_ExpandedNuclei.csv")
#d54.week1.donut.donuts.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week1/output/D54_week1_donuts_Donuts.csv")

# Donut Analysis: D54 Week 2 Diff, pTRIPZ-Control or pTRIPZ-4707, +DOX, DAPI/GFP/TxRed(NES)/Cy7(Tuj1)
#d54.week2.donut.nuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_donuts_Nuclei.csv")
#d54.week2.donut.expandednuclei.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_donuts_ExpandedNuclei.csv")
#d54.week2.donut.donuts.csv <- here("results/microscopy/2022_Diff_Experiment/D54_week2/output/D54_week2_donuts_Donuts.csv")


# GLOBALS ##############################################################################################################


# Plate Layout #############################################################################################

# plate layout metadata
df.d54.w1.plate <- read_csv(d54.week1.plate.layout.csv)
df.d54.w2.plate <- read_csv(d54.week2.plate.layout.csv)

df.d54.w1.plate %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))
df.d54.w2.plate %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))

df.d54.w1.plate %<>%
    mutate(Week = factor(Week,
                         levels = c(1, 2),
                         labels = c("Week 1", "Week 2"),
                         ordered = TRUE))
df.d54.w2.plate %<>%
    mutate(Week = factor(Week,
                         levels = c(1, 2),
                         labels = c("Week 1", "Week 2"),
                         ordered = TRUE))

df.d54.w1.plate %<>%
    mutate(Expression = factor(Expression,
                         levels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                         ordered = TRUE))
df.d54.w2.plate %<>%
    mutate(Expression = factor(Expression,
                         levels = c("pTRIPZ-Control", "pTRIPZ-4707"),
                         ordered = TRUE))

# Images ##############
# images
df.d54.w1.images <- read_csv(d54.week1.images.csv)
df.d54.w2.images <- read_csv(d54.week2.images.csv)


# make row and column
df.d54.w1.images %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.images %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Count_"),
           starts_with("Intensity_"))
df.d54.w2.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Count_"),
           starts_with("Intensity_"))

df.d54.w1.images %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.images %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))

df.d54.images <- bind_rows(df.d54.w1.images, df.d54.w2.images)

df.d54.images %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))

df.d54.images %<>%
    mutate(Week = factor(Week))

df.d54.images %>%
    filter(! ICC_control) %>%
    group_by(Well, Donor, Expression, Week) %>%
    summarise(count_nuclei = sum(Count_Nuclei)) -> df.d54.images.by.well



df.d54.images.by.well %>%
    ggplot(aes(x = Week, y = count_nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.2), size = 1) +
    plotTheme("figure") +
    stat_compare_means(method = "t.test", label = "p.format") +
    scale_color_manual(values = paperCatEight) +
    theme(legend.position = "none")

ggsave("~/Desktop/d54_nuclei_count.pdf", height = 2, width = 2)





df.d54.images %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Week, y = Count_Nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")

df.d54.images %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Week, y = Intensity_TotalIntensity_GFP_raw / Count_Nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")

df.d54.images %>%
    filter(! ICC_control, (Intensity_TotalIntensity_Cy7_raw / Count_Nuclei)< 20000) %>%
    ggplot(aes(x = Week, y = Intensity_TotalIntensity_TxRed_raw / Count_Nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")

df.d54.images %>%
    filter(! ICC_control, (Intensity_TotalIntensity_Cy7_raw / Count_Nuclei)< 20000) %>%
    ggplot(aes(x = Week, y = Intensity_TotalIntensity_Cy7_raw / Count_Nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")


df.d54.images %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Week, y = Intensity_TotalIntensity_TxRed_raw / Intensity_TotalIntensity_GFP_raw, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")

df.d54.images %>%
    filter(! ICC_control, (Intensity_TotalIntensity_Cy7_raw / Count_Nuclei)< 20000) %>%
    ggplot(aes(x = Week, y = Intensity_TotalIntensity_Cy7_raw / (Intensity_TotalIntensity_GFP_raw ), color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    stat_compare_means(method = "t.test", label = "p.format")


# Nuclei ##############
# Nuclei
df.d54.w1.nuclei <- read_csv(d54.week1.nuclei.csv)
df.d54.w2.nuclei <- read_csv(d54.week2.nuclei.csv)

df.d54.w1.nuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.nuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.nuclei %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.nuclei %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.nuclei) == colnames(df.d54.w2.nuclei))

df.d54.nuclei <- bind_rows(df.d54.w1.nuclei, df.d54.w2.nuclei)

df.d54.nuclei %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))

df.d54.nuclei %<>%
    mutate(Week = factor(Week))


# Cells ##################

# GFP cells
# Cy7 cells
df.d54.w1.GFPcells <- read_csv(d54.week1.GFPcells.csv)
df.d54.w2.GFPcells <- read_csv(d54.week2.GFPcells.csv)

df.d54.w1.GFPcells %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.GFPcells %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.GFPcells %<>%
    select(ImageNumber,
           ObjectNumber,
           Parent_Nuclei,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.GFPcells %<>%
    select(ImageNumber,
           ObjectNumber,
           Parent_Nuclei,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.GFPcells %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.GFPcells %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.GFPcells) == colnames(df.d54.w2.GFPcells))

df.d54.GFPcells <- bind_rows(df.d54.w1.GFPcells, df.d54.w2.GFPcells)

df.d54.GFPcells %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))

df.d54.GFPcells %<>%
    mutate(Week = factor(Week))

colnames(df.d54.GFPcells)



df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_TxRed_raw / Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %<>%
    mutate(Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_GFP_raw,
           Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw = Intensity_MeanIntensity_TxRed_raw / Intensity_MeanIntensity_GFP_raw)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_point(size = 1, shape = 1) +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw / Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x =  AreaShape_FormFactor, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.GFPcells %<>%
    mutate(GFP_positive = AreaShape_Area)

# combine with nuclei
df.d54.nuclei %>%
    filter(! ICC_control) %>%
    select(ImageNumber, ObjectNumber, AreaShape_Area.nuclei = AreaShape_Area, Week) %>%
    mutate(object = paste(Week, ImageNumber, ObjectNumber, sep = "_")) %>%
    select(object, AreaShape_Area.nuclei) -> df.nuc.tmp

df.d54.GFPcells %<>%
    filter(! ICC_control) %>%
    mutate(object = paste(Week, ImageNumber, Parent_Nuclei, sep = "_"))

df.d54.GFPcells %<>%
    left_join(df.nuc.tmp, by = "object")

df.d54.GFPcells %<>%
    mutate(GFP_size = AreaShape_Area - AreaShape_Area.nuclei,
           GFP_size_ratio = AreaShape_Area / AreaShape_Area.nuclei)

df.d54.GFPcells %<>%
    mutate(#GFP_positive =  GFP_size_ratio >= 1.1,
        GFP_positive = GFP_size > 0)

# classify cells
df.d54.GFPcells %>%
    filter(Expression == "pTRIPZ-Control", Week == "Week 2", GFP_positive) %>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed_GFPnormalized = Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw / Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw) %>%
    slice_max(Intensity_MeanIntensity_Cy7_over_TxRed_GFPnormalized, prop = 0.2) -> df.tmp

RATIO_CUTOFF <- min(df.tmp$Intensity_MeanIntensity_Cy7_over_TxRed_GFPnormalized)

df.d54.GFPcells %<>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed_GFPnormalized = Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw / Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw) %>%
    mutate(Neuron = Intensity_MeanIntensity_Cy7_over_TxRed_GFPnormalized >= RATIO_CUTOFF)

pdf("~/Desktop/gfp_cells.pdf")

df.d54.GFPcells %>%
    filter(GFP_positive) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw, color = Neuron)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "GFP cells")

df.d54.GFPcells %>%
    filter(GFP_positive) %>%
    group_by(Well, Row, Column, Donor, Expression, Week, ICC_control) %>%
    summarise(Count_Neuron = sum(Neuron),
              Count_Projenitor = sum(!Neuron),
              Count_Cells = n(),
              Mean_GFP = mean(Intensity_MeanIntensity_GFP_raw),
              Mean_TxRed = mean(Intensity_MeanIntensity_TxRed_raw),
              Mean_Cy7 = mean(Intensity_MeanIntensity_Cy7_raw),
              Mean_TxRed_gfp = mean(Intensity_MeanIntensity_TxRed_raw_over_Intensity_MeanIntensity_GFP_raw),
              Mean_Cy7_gfp = mean(Intensity_MeanIntensity_Cy7_raw_over_Intensity_MeanIntensity_GFP_raw),
              Mean_Area = mean(AreaShape_Area),
              Mean_ConvexArea = mean(AreaShape_ConvexArea)) -> df.wells.GFPcells

df.wells.GFPcells %>%
    left_join(df.d54.images.by.well) -> df.wells



df.wells %>%
    ggplot(aes(x = Expression, y = Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells %>%
    ggplot(aes(x = Expression, y = count_nuclei, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Nuclei, by well")

df.wells %>%
    ggplot(aes(x = Expression, y = Count_Cells / count_nuclei, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells %>%
    ggplot(aes(x = Expression, y = Count_Neuron, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells %>%
    ggplot(aes(x = Expression, y = Count_Neuron / Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

# for paper
df.wells %>%
    ggplot(aes(x = Expression, y = Count_Neuron / count_nuclei, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells %>%
    ggplot(aes(x = Week, y = Count_Neuron / count_nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.2), size = 1) +
    plotTheme("figure") +
    stat_compare_means(method = "t.test", label = "p.format") +
    scale_color_manual(values = paperCatEight)

ggsave("~/Desktop/d54_neuron_by_nuclei_legend.pdf", height = 2, width = 2)


df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_Area, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_GFP, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_TxRed, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_Cy7, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_TxRed_gfp, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

df.wells.GFPcells %>%
    ggplot(aes(x = Expression, y = Mean_Cy7_gfp, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "GFP cells, by well")

dev.off()

# Cy7 cells
df.d54.w1.Cy7cells <- read_csv(d54.week1.Cy7cells.csv)
df.d54.w2.Cy7cells <- read_csv(d54.week2.Cy7cells.csv)

df.d54.w1.Cy7cells %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.Cy7cells %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.Cy7cells %<>%
    select(ImageNumber,
           ObjectNumber,
           Parent_Nuclei,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.Cy7cells %<>%
    select(ImageNumber,
           ObjectNumber,
           Parent_Nuclei,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.Cy7cells %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.Cy7cells %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.Cy7cells) == colnames(df.d54.w2.Cy7cells))

df.d54.Cy7cells <- bind_rows(df.d54.w1.Cy7cells, df.d54.w2.Cy7cells)

df.d54.Cy7cells %<>%
    mutate(ICC_control = ! (Primary_GFP_Chicken & Primary_NES_Rabbit & Primary_TUJ1_Mouse & Secondary_Chicken_AF488 & Secondary_Rabbit_AF568 & Secondary_Mouse_AF750))

df.d54.Cy7cells %<>%
    mutate(Week = factor(Week))

colnames(df.d54.Cy7cells)

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_TxRed_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    geom_histogram(position = "dodge") +
    facet_wrap(~Week)

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_Cy7_raw, y = Intensity_IntegratedIntensity_TxRed_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFP_raw, y = Intensity_IntegratedIntensity_Cy7_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFP_raw, y = Intensity_IntegratedIntensity_TxRed_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_Cy7_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    geom_point(size = 0.5) +
    scale_y_log10() +
    scale_x_log10()



# scratch
df.d54.nuclei %>%
    filter(! ICC_control) %>%
    select(ImageNumber, ObjectNumber, AreaShape_Area.nuclei = AreaShape_Area, Week) %>%
    mutate(object = paste(Week, ImageNumber, ObjectNumber, sep = "_")) %>%
    select(object, AreaShape_Area.nuclei) -> df.nuc.tmp

df.d54.Cy7cells %>%
    filter(! ICC_control) %>%
    mutate(object = paste(Week, ImageNumber, Parent_Nuclei, sep = "_")) -> df.cy7.tmp

df.tmp <- left_join(df.cy7.tmp, df.nuc.tmp, by = "object")

df.tmp %<>%
    mutate(cy7_size = AreaShape_Area - AreaShape_Area.nuclei)

df.tmp %<>%
    filter(cy7_size >= 0)

df.tmp %>%
    ggplot(aes(x = cy7_size, color = Expression)) +
    geom_density()

df.tmp %>%
    ggplot(aes(x = AreaShape_FormFactor, color = Expression)) +
    geom_density()

df.tmp %>%
    ggplot(aes(x = AreaShape_FormFactor, y = AreaShape_Area, color = Expression)) +
    geom_point()

df.tmp %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, color = Expression)) +
    geom_density()

df.tmp %>%
    filter(AreaShape_FormFactor < 0.6) %>%
    ggplot(aes(x = Intensity_MedianIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    facet_wrap(~Week) +
    geom_point()

df.tmp %>%
    filter(AreaShape_FormFactor < 0.5) %>%
    group_by(Well, Donor, Expression, Week) %>%
    summarise(count_cy7_cells = n()) %>%
    left_join(df.d54.images.by.well) %>%
    ggplot(aes(x = Week, y = count_cy7_cells / count_nuclei, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.75, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 1) +
    plotTheme() +
    stat_compare_means(method = "t.test", label = "p.format")

df.tmp %>%
    filter(AreaShape_FormFactor < 0.5) %>%
    ggplot(aes(x = Intensity_MeanIntensity_TxRed_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

df.tmp %>%
    filter(AreaShape_FormFactor < 0.5) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, color = Expression)) +
    geom_density() +
    facet_wrap(~Week)

# Donuts ###############

# donuts
df.d54.w1.donuts <- read_csv(d54.week1.donut.donuts.csv)
df.d54.w2.donuts <- read_csv(d54.week2.donut.donuts.csv)

df.d54.w1.donuts %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.donuts %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.donuts %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.donuts %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.donuts %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.donuts %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.donuts) == colnames(df.d54.w2.donuts))

df.d54.donuts <- bind_rows(df.d54.w1.donuts, df.d54.w2.donuts)
rm(df.d54.w1.donuts, df.d54.w2.donuts)

# expanded nuclei
df.d54.w1.expnuclei <- read_csv(d54.week1.donut.expandednuclei.csv)
df.d54.w2.expnuclei <- read_csv(d54.week2.donut.expandednuclei.csv)

df.d54.w1.expnuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.expnuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.expnuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.expnuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.expnuclei %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.expnuclei %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.expnuclei) == colnames(df.d54.w2.expnuclei))

df.d54.expnuclei <- bind_rows(df.d54.w1.expnuclei, df.d54.w2.expnuclei)
rm(df.d54.w1.expnuclei, df.d54.w2.expnuclei)


# nuclei
df.d54.w1.nuclei <- read_csv(d54.week1.donut.nuclei.csv)
df.d54.w2.nuclei <- read_csv(d54.week2.donut.nuclei.csv)

df.d54.w1.nuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))
df.d54.w2.nuclei %<>%
    mutate(Row = sapply(strsplit(Metadata_Well, "[0-9]"), `[`, 1),
           Column = as.numeric(sapply(strsplit(Metadata_Well, "[A-Z]"), `[`, 2)))

df.d54.w1.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))
df.d54.w2.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row,
           Column,
           starts_with("Intensity"),
           starts_with("AreaShape"))

df.d54.w1.nuclei %<>%
    left_join(df.d54.w1.plate, by = c("Row", "Column"))
df.d54.w2.nuclei %<>%
    left_join(df.d54.w2.plate, by = c("Row", "Column"))



all(colnames(df.d54.w1.nuclei) == colnames(df.d54.w2.nuclei))

df.d54.nuclei <- bind_rows(df.d54.w1.nuclei, df.d54.w2.nuclei)
rm(df.d54.w1.nuclei, df.d54.w2.nuclei)

pdf(paste0(dir.pdfs, "20220815_donut_analysis.pdf"))
# donuts
# classify cells
df.d54.donuts %>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    filter(Week == "Week 2") %>%
    slice_max(Intensity_MeanIntensity_Cy7_over_TxRed, prop = 0.2) -> df.tmp

RATIO_CUTOFF <- 4.722982

df.d54.donuts %<>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    mutate(Neuron = Intensity_MeanIntensity_Cy7_over_TxRed >= RATIO_CUTOFF)


df.d54.donuts %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed)) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Neuron)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Donuts")

df.d54.donuts %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed), ! is.nan(AreaShape_Area)) %>%
    group_by(Well, Row, Column, Donor, Expression, Week, ICC_control) %>%
    summarise(Count_Neuron = sum(Neuron),
              Count_Projenitor = sum(!Neuron),
              Count_Cells = n(),
              Mean_DAPI = mean(Intensity_MeanIntensity_DAPI_raw),
              Mean_GFP = mean(Intensity_MeanIntensity_GFP_raw),
              Mean_TxRed = mean(Intensity_MeanIntensity_TxRed_raw),
              Mean_Cy7 = mean(Intensity_MeanIntensity_Cy7_raw),
              Mean_Area = mean(AreaShape_Area),
              Mean_ConvexArea = mean(AreaShape_ConvexArea)) -> df.wells.donuts

df.wells.donuts %>%
    ggplot(aes(x = Expression, y = Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Donuts, by well")

df.wells.donuts %>%
    ggplot(aes(x = Expression, y = Count_Neuron / Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Donuts, by well")

df.wells.donuts %>%
    ggplot(aes(x = Expression, y = Mean_Area, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Donuts, by well")



# ext nuclei
# classify cells
df.d54.expnuclei %>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    filter(Week == "Week 2") %>%
    slice_max(Intensity_MeanIntensity_Cy7_over_TxRed, prop = 0.2) -> df.tmp

RATIO_CUTOFF <- 5.078791

df.d54.expnuclei %<>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    mutate(Neuron = Intensity_MeanIntensity_Cy7_over_TxRed >= RATIO_CUTOFF)

df.d54.expnuclei %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed)) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Neuron)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Expanded Nuclei")

df.d54.expnuclei %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed), ! is.nan(AreaShape_Area)) %>%
    group_by(Well, Row, Column, Donor, Expression, Week, ICC_control) %>%
    summarise(Count_Neuron = sum(Neuron),
              Count_Projenitor = sum(!Neuron),
              Count_Cells = n(),
              Mean_DAPI = mean(Intensity_MeanIntensity_DAPI_raw),
              Mean_GFP = mean(Intensity_MeanIntensity_GFP_raw),
              Mean_TxRed = mean(Intensity_MeanIntensity_TxRed_raw),
              Mean_Cy7 = mean(Intensity_MeanIntensity_Cy7_raw),
              Mean_Area = mean(AreaShape_Area),
              Mean_ConvexArea = mean(AreaShape_ConvexArea)) -> df.wells.expanded

df.wells.expanded %>%
    ggplot(aes(x = Expression, y = Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Expanded Nuclei, by well")

df.wells.expanded %>%
    ggplot(aes(x = Expression, y = Count_Neuron / Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Expanded Nuclei, by well")

df.wells.expanded %>%
    ggplot(aes(x = Expression, y = Mean_Area, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Expanded Nuclei, by well")


#  nuclei
# classify cells
df.d54.nuclei %>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    filter(Week == "Week 2") %>%
    slice_max(Intensity_MeanIntensity_Cy7_over_TxRed, prop = 0.2) -> df.tmp

RATIO_CUTOFF <- 5.286366

df.d54.nuclei %<>%
    mutate(Intensity_MeanIntensity_Cy7_over_TxRed = Intensity_MeanIntensity_Cy7_raw / Intensity_MeanIntensity_TxRed_raw) %>%
    mutate(Neuron = Intensity_MeanIntensity_Cy7_over_TxRed >= RATIO_CUTOFF)

df.d54.nuclei %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed)) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = Neuron)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Nuclei")

df.d54.nuclei %>%
    filter(! ICC_control, ! is.na(Intensity_MeanIntensity_Cy7_over_TxRed), ! is.nan(AreaShape_Area)) %>%
    group_by(Well, Row, Column, Donor, Expression, Week, ICC_control) %>%
    summarise(Count_Neuron = sum(Neuron),
              Count_Projenitor = sum(!Neuron),
              Count_Cells = n(),
              Mean_DAPI = mean(Intensity_MeanIntensity_DAPI_raw),
              Mean_GFP = mean(Intensity_MeanIntensity_GFP_raw),
              Mean_TxRed = mean(Intensity_MeanIntensity_TxRed_raw),
              Mean_Cy7 = mean(Intensity_MeanIntensity_Cy7_raw),
              Mean_Area = mean(AreaShape_Area),
              Mean_ConvexArea = mean(AreaShape_ConvexArea)) -> df.wells.nuclei

df.wells.nuclei %>%
    ggplot(aes(x = Expression, y = Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Nuclei, by well")

df.wells.nuclei %>%
    ggplot(aes(x = Expression, y = Count_Neuron / Count_Cells, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Nuclei, by well")

df.wells.nuclei %>%
    ggplot(aes(x = Expression, y = Mean_Area, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Week, scales = "free_y") +
    stat_compare_means(method = "t.test", label = "p.format") +
    labs(title = "Nuclei, by well")



#dev.off()

df.d54.donuts %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Donuts")

df.d54.donuts %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Donuts")

df.d54.donuts %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_Cy7_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Donuts")


df.d54.nuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Nuclei")

df.d54.nuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Nuclei")

df.d54.nuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_Cy7_raw)) +
    geom_point(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Nuclei")


df.d54.expnuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_hex(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Expanded Nuclei")

df.d54.expnuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_TxRed_raw)) +
    geom_hex(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Expanded Nuclei")

df.d54.expnuclei %>%
    filter(! ICC_control) %>%
    ggplot(aes(x = Intensity_MeanIntensity_GFP_raw, y = Intensity_MeanIntensity_Cy7_raw)) +
    geom_hex(size = 0.1) +
    facet_wrap(~Expression + Week) +
    labs(title = "Expanded Nuclei")


dev.off()



# df.d54.GFPcells %<>%
#     mutate(TxRed_only = Intensity_MeanIntensity_Cy7_raw < .022 & Intensity_MeanIntensity_TxRed_raw > .015,
#            Cy7_only = Intensity_MeanIntensity_Cy7_raw > .025 & Intensity_MeanIntensity_TxRed_raw < .015)
#
# df.d54.GFPcells %>%
#     ggplot(aes(x = Intensity_MeanIntensity_Cy7_raw, y = Intensity_MeanIntensity_TxRed_raw, color = TxRed_only)) +
#     geom_point(size = 0.5) +
#     scale_y_log10() +
#     scale_x_log10()
#
# df.d54.GFPcells %>%
#     filter(! ICC_control) %>%
#     filter(TxRed_only | Cy7_only) -> df.d54.GFPcells.positive
#
# df.d54.GFPcells.positive %>%
#     group_by(Well, Row, Column, Donor, Expression, Week)

# df.w2.images <- read_csv(week2.images.csv)
#
# df.w2.images %<>%
#     select(ImageNumber,
#            Well = Metadata_Well,
#            Row = Metadata_Row,
#            Column = Metadata_Column,
#            Count_Cy5objects,
#            Count_GFPobjects,
#            Count_Nuclei) %>%
#     mutate(Column = as.numeric(Column))
#
# df.w2.images %<>%
#     left_join(df.w2.plate, by = c("Row", "Column"))
#
# all(colnames(df.w1.images) == colnames(df.w2.images))
#
# df.w1.images %>%
#     bind_rows(df.w2.images) -> df.images
#
# rm(df.w1.images, df.w2.images)
#
# # summarize by well
# df.images %>%
#     group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
#     summarise(Count_Cy5objects = sum(Count_Cy5objects),
#               Count_GFPobjects = sum(Count_GFPobjects),
#               Count_Nuclei = sum(Count_Nuclei)) -> df.wells
#
# df.wells$Count_GFPobjects[df.wells$Expression == "Untransduced"] <- NA
#
# df.wells %<>%
#     mutate(Cy5objects_over_nuclei = Count_Cy5objects / Count_Nuclei,
#            GFPobjects_over_nuclei = Count_GFPobjects / Count_Nuclei,
#            Week = factor(Week))
#
# df.wells %>%
#     filter(!Expression == "Untransduced") %>%
#     ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
# stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
#
# df.wells %>%
#     filter(!Expression == "Untransduced") %>%
#     ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray"))
#
# df.wells %>%
#     filter(Expression == "Untransduced" | Count_GFPobjects < 1200) %>%
#     ggplot(mapping = aes(y = Count_GFPobjects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray"))
#
# df.wells %>%
#     filter(Expression == "Untransduced" | Count_GFPobjects < 1200) %>%
#     ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
#
#
# # Cy5 objects
# df.w1.Cy5objects <- read_csv(week1.Cy5objects.csv)
#
# df.w1.Cy5objects %<>%
#     select(ImageNumber,
#            ObjectNumber,
#            Well = Metadata_Well,
#            Row = Metadata_Row,
#            Column = Metadata_Column,
#            starts_with("Intensity")) %>%
#     mutate(Column = as.numeric(Column))
#
# df.w1.Cy5objects %<>%
#     left_join(df.w1.plate, by = c("Row", "Column"))
#
# df.w2.Cy5objects <- read_csv(week2.Cy5objects.csv)
#
# df.w2.Cy5objects %<>%
#     select(ImageNumber,
#            ObjectNumber,
#            Well = Metadata_Well,
#            Row = Metadata_Row,
#            Column = Metadata_Column,
#            starts_with("Intensity")) %>%
#     mutate(Column = as.numeric(Column))
#
# df.w2.Cy5objects %<>%
#     left_join(df.w2.plate, by = c("Row", "Column"))
#
# all(colnames(df.w1.Cy5objects) == colnames(df.w2.Cy5objects))
#
# df.w1.Cy5objects %>%
#     bind_rows(df.w2.Cy5objects) -> df.Cy5objects
#
# rm(df.w1.Cy5objects, df.w2.Cy5objects)
#






# Plot #################################################################################################################

df.d54.w1.plate %>%
    ggplot(mapping = aes(x = Column, y = Row.plot)) +
    geom_point(data = expand.grid(seq(1, 12), seq(1, 8)), mapping = aes(x = Var1, y = Var2), color="grey90", fill="white", shape=21, size=20) +
    geom_point(aes(color = Primary_NES_Rabbit)) +
    coord_fixed(ratio=(13/12)/(9/8), xlim=c(0.5, 12.5), ylim=c(8.5, 0.5), expand = FALSE) +
    scale_y_reverse(breaks=seq(1, 8), labels=LETTERS[1:8]) +
    scale_x_continuous(breaks=seq(1, 12), position = "top") +
    labs(title="Plate Layout for My Experiment")






