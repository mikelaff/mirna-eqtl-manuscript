
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
# D54/D88 Week 1 Diff, pTRIPZ-4707 stable line, +/- DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week1.nuclei.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Nuclei.csv")
week1.Cy5objects.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Cy5objects.csv")
week1.images.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Image.csv")
week1.plate.layout.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-4707stbl_GFAP_Tuj1/plate_layout.csv")

# D54/D88 Week 2 Diff, pTRIPZ-4707 stable line, +/- DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week2.nuclei.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Nuclei.csv")
week2.Cy5objects.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Cy5objects.csv")
week2.images.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Image.csv")
week2.plate.layout.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/plate_layout.csv")

# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################

# plate layout metadata
df.w1.plate <- read_csv(week1.plate.layout.csv)
df.w2.plate <- read_csv(week2.plate.layout.csv)

# images
df.w1.images <- read_csv(week1.images.csv)

df.w1.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           Count_Cy5objects,
           Count_Nuclei) %>%
    mutate(Column = as.numeric(Column))

df.w1.images %<>%
    left_join(df.w1.plate, by = c("Row", "Column"))

df.w2.images <- read_csv(week2.images.csv)

df.w2.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           Count_Cy5objects,
           Count_Nuclei) %>%
    mutate(Column = as.numeric(Column))

df.w2.images %<>%
    left_join(df.w2.plate, by = c("Row", "Column"))

all(colnames(df.w1.images) == colnames(df.w2.images))

df.w1.images %>%
    bind_rows(df.w2.images) -> df.images

rm(df.w1.images, df.w2.images)

# summarize by well
df.images %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
    summarise(Count_Cy5objects = sum(Count_Cy5objects),
              Count_Nuclei = sum(Count_Nuclei)) -> df.wells

df.wells %<>%
    mutate(Cy5objects_over_nuclei = Count_Cy5objects / Count_Nuclei,
           Week = factor(Week))

df.wells %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.wells %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.wells %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")



# Cy5 objects
df.w1.Cy5objects <- read_csv(week1.Cy5objects.csv)

df.w1.Cy5objects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.w1.Cy5objects %<>%
    left_join(df.w1.plate, by = c("Row", "Column"))

df.w2.Cy5objects <- read_csv(week2.Cy5objects.csv)

df.w2.Cy5objects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.w2.Cy5objects %<>%
    left_join(df.w2.plate, by = c("Row", "Column"))

all(colnames(df.w1.Cy5objects) == colnames(df.w2.Cy5objects))

df.w1.Cy5objects %>%
    bind_rows(df.w2.Cy5objects) -> df.Cy5objects

rm(df.w1.Cy5objects, df.w2.Cy5objects)

# only cy5 objects with a max intensity > 0.8
df.Cy5objects %<>%
    mutate(Cy5positive_.5 = Intensity_MaxIntensity_Cy5 > 0.5,
           Cy5positive_.8 = Intensity_MaxIntensity_Cy5 > 0.8,
           Cy5positive_.9 = Intensity_MaxIntensity_Cy5 > 0.9)

df.Cy5objects %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
    summarise(Count_Cy5objects = n(),
              Cy5positive_.5 = sum(Cy5positive_.5),
              Cy5positive_.8 = sum(Cy5positive_.8),
              Cy5positive_.9 = sum(Cy5positive_.9)) -> df.Cy5objects.summary

df.Cy5objects.summary %<>%
    mutate(Week = factor(Week))

df.Cy5objects.summary %<>%
    left_join(df.wells, by = c("Well", "Row", "Column", "Donor", "Expression", "DOX", "Week", "Primary_AB", "Count_Cy5objects"))

df.Cy5objects.summary %<>%
    mutate(Cy5positive_.5_over_nuclei = Cy5positive_.5 / Count_Nuclei,
           Cy5positive_.8_over_nuclei = Cy5positive_.8 / Count_Nuclei,
           Cy5positive_.9_over_nuclei = Cy5positive_.9 / Count_Nuclei)

pdf("~/Desktop/d54_d88_diff_4707stable.pdf")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.5, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.8, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.9, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.5_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.8_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.9_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-4707 stable",
         caption = "24 wells/condition")

dev.off()


# Figure ###########################

dir.pdfs <- here("doc/paper/figure6/pdfs/")

df.Cy5objects.summary %>%
    ggplot(aes(y = Count_Nuclei, x = Week, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
    facet_wrap(~Donor, nrow = 1, scales = "fixed") +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c(paperBlue, paperRed)) +
    labs(y = "Nuclei Count",
         x = "Week") +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D88_D54_diff_nuclei_figure_withSig.pdf"), height = 2, width = 3)

df.Cy5objects.summary %>%
    ggplot(aes(y = Cy5objects_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
    facet_wrap(~Primary_AB + Donor, nrow = 1, scales = "fixed") +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Cy5 Positive / Nuclei",
         x = "Week") +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D88_D54_diff_cy5_figure_withSig.pdf"), height = 2, width = 4)


