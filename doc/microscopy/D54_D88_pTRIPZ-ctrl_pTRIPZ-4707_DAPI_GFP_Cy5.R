
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
week1.nuclei.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Nuclei.csv")
week1.Cy5objects.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Cy5objects.csv")
week1.GFPobjects.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_GFPobjects.csv")
week1.images.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Image.csv")
week1.plate.layout.csv <- here("results/microscopy/D54_D88_Diff1week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/plate_layout.csv")

# D54/D88 Week 2 Diff, pTRIPZ-ctrl or pTRIPZ-4707 low MOI, +DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
week2.nuclei.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Nuclei.csv")
week2.Cy5objects.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Cy5objects.csv")
week2.GFPobjects.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_GFPobjects.csv")
week2.images.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Image.csv")
week2.plate.layout.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-ctrl_pTRIPZ-4707_GFAP_Tuj1/plate_layout.csv")

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
           Count_GFPobjects,
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
           Count_GFPobjects,
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
              Count_GFPobjects = sum(Count_GFPobjects),
              Count_Nuclei = sum(Count_Nuclei)) -> df.wells

df.wells$Count_GFPobjects[df.wells$Expression == "Untransduced"] <- NA

df.wells %<>%
    mutate(Cy5objects_over_nuclei = Count_Cy5objects / Count_Nuclei,
           GFPobjects_over_nuclei = Count_GFPobjects / Count_Nuclei,
           Week = factor(Week))

df.wells %>%
    filter(!Expression == "Untransduced") %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")


df.wells %>%
    filter(!Expression == "Untransduced") %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray"))

df.wells %>%
    filter(Expression == "Untransduced" | Count_GFPobjects < 1200) %>%
    ggplot(mapping = aes(y = Count_GFPobjects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray"))

df.wells %>%
    filter(Expression == "Untransduced" | Count_GFPobjects < 1200) %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
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

# GFP objects
df.w1.GFPobjects <- read_csv(week1.GFPobjects.csv)

df.w1.GFPobjects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.w1.GFPobjects %<>%
    left_join(df.w1.plate, by = c("Row", "Column"))

df.w2.GFPobjects <- read_csv(week2.GFPobjects.csv)

df.w2.GFPobjects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.w2.GFPobjects %<>%
    left_join(df.w2.plate, by = c("Row", "Column"))

all(colnames(df.w1.GFPobjects) == colnames(df.w2.GFPobjects))

df.w1.GFPobjects %>%
    bind_rows(df.w2.GFPobjects) -> df.GFPobjects

rm(df.w1.GFPobjects, df.w2.GFPobjects)


df.wells %<>%
    filter(! (Week == 2 & Well == "E11")) %>%
    filter(! (Week == 2 & Well == "E09")) %>%
    filter(! (Week == 2 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "E03")) %>%
    filter(! (Week == 1 & Well == "B02")) %>%
    filter(! (Week == 1 & Well == "E11")) %>%
    filter(! (Week == 1 & Well == "B11")) %>%
    filter(! (Week == 1 & Well == "B09")) %>%
    filter(! (Week == 1 & Well == "E10")) %>%
    filter(! (Week == 1 & Well == "B04")) %>%
    filter(! (Week == 1 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "B08"))

df.GFPobjects %<>%
    filter(! (Week == 2 & Well == "E11")) %>%
    filter(! (Week == 2 & Well == "E09")) %>%
    filter(! (Week == 2 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "E03")) %>%
    filter(! (Week == 1 & Well == "B02")) %>%
    filter(! (Week == 1 & Well == "E11")) %>%
    filter(! (Week == 1 & Well == "B11")) %>%
    filter(! (Week == 1 & Well == "B09")) %>%
    filter(! (Week == 1 & Well == "E10")) %>%
    filter(! (Week == 1 & Well == "B04")) %>%
    filter(! (Week == 1 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "B08"))

df.Cy5objects %<>%
    filter(! (Week == 2 & Well == "E11")) %>%
    filter(! (Week == 2 & Well == "E09")) %>%
    filter(! (Week == 2 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "E03")) %>%
    filter(! (Week == 1 & Well == "B02")) %>%
    filter(! (Week == 1 & Well == "E11")) %>%
    filter(! (Week == 1 & Well == "B11")) %>%
    filter(! (Week == 1 & Well == "B09")) %>%
    filter(! (Week == 1 & Well == "E10")) %>%
    filter(! (Week == 1 & Well == "B04")) %>%
    filter(! (Week == 1 & Well == "E02")) %>%
    filter(! (Week == 1 & Well == "B08"))

df.GFPobjects %<>%
    filter(! (Expression == "Untransduced"))

df.GFPobjects %>%
    ggplot(aes(x = Intensity_MaxIntensity_GFP)) +
    geom_histogram(bins = 50)

df.GFPobjects %>%
    ggplot(aes(x = Intensity_MaxIntensity_Cy5)) +
    geom_histogram(bins = 50)

df.GFPobjects %<>%
    mutate(Cy5positive_GFPobjects = Intensity_MaxIntensity_Cy5 > .9)

df.Cy5objects %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
    summarise(Count_Cy5objects = n()) -> df.Cy5objects.summary

df.GFPobjects %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
    summarise(Count_GFPobjects = n(),
              Count_Cy5positive_GFPobjects = sum(Cy5positive_GFPobjects)) -> df.GFPobjects.summary

df.Cy5objects.summary %<>%
    mutate(Week = factor(Week))

df.GFPobjects.summary %<>%
    mutate(Week = factor(Week))

df.Cy5objects.summary %<>%
    left_join(df.wells, by = c("Well", "Row", "Column", "Donor", "Expression", "DOX", "Week", "Primary_AB", "Count_Cy5objects"))

df.GFPobjects.summary %<>%
    left_join(df.wells, by = c("Well", "Row", "Column", "Donor", "Expression", "DOX", "Week", "Primary_AB", "Count_GFPobjects"))

df.GFPobjects.summary %<>%
    mutate(Cy5positive_GFPobjects_over_nuclei = Count_Cy5positive_GFPobjects / Count_Nuclei)

df.GFPobjects.summary %<>%
    mutate(Cy5positive_GFPobjects_over_GFPobjects = Count_Cy5positive_GFPobjects / Count_GFPobjects)

df.GFPobjects.summary %<>%
    mutate(Expression = factor(Expression, levels = c("pTRIPZ-Ctrl", "pTRIPZ-4707"), ordered = TRUE))



pdf("~/Desktop/d54_d88_diff_ctrl_4707.pdf")

# df.Cy5objects.summary %>%
#     filter(! Expression == "Untransduced") %>%
#     ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "16 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "8 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Count_GFPobjects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "8 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = GFPobjects_over_nuclei, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "16 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Count_Cy5positive_GFPobjects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "8 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_GFPobjects_over_GFPobjects, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "8 wells/condition")

df.GFPobjects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_GFPobjects_over_nuclei, x = Week, color = Expression)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Diff. pTRIPZ-ctrl/4707 Low MOI",
         caption = "8 wells/condition")

dev.off()









df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.5, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.8, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.9, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.5_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.8_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5positive_.9_over_nuclei, x = Week, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB + Donor) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

dev.off()
