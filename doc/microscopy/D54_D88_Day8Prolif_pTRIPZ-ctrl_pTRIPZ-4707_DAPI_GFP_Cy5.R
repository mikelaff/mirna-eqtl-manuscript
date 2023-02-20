
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
# D54 Day 8 prolif, pTRIPZ-ctrl or pTRIPZ-4707, +DOX, DAPI/GFP/Cy5 (EdU)
d54.nuclei.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D54_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Nuclei.csv")
d54.Cy5objects.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D54_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Cy5objects.csv")
d54.GFPobjects.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D54_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_GFPobjects.csv")
d54.images.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D54_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Image.csv")
d54.plate.layout.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/D54_plate_layout.csv")

# D88 Day 8 prolif, pTRIPZ-ctrl or pTRIPZ-4707, +DOX, DAPI/GFP/Cy5 (EdU)
d88.nuclei.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D88_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Nuclei.csv")
d88.Cy5objects.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D88_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Cy5objects.csv")
d88.GFPobjects.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D88_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_GFPobjects.csv")
d88.images.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/output/D88_pTRIPZ-ctrl_pTRIPZ-4707_DAPI_GFP_Cy5_Image.csv")
d88.plate.layout.csv <- here("results/microscopy/D54_D88_Prolif8days_pTRIPZ-ctrl_pTRIPZ-4707_EdU/D88_plate_layout.csv")

# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################

# plate layout metadata
df.d54.plate <- read_csv(d54.plate.layout.csv)
df.d88.plate <- read_csv(d88.plate.layout.csv)

# images
df.d54.images <- read_csv(d54.images.csv)

df.d54.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           Count_Cy5objects,
           Count_GFPobjects,
           Count_Nuclei) %>%
    mutate(Column = as.numeric(Column))

df.d54.images %<>%
    left_join(df.d54.plate, by = c("Row", "Column"))

df.d88.images <- read_csv(d88.images.csv)

df.d88.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           Count_Cy5objects,
           Count_GFPobjects,
           Count_Nuclei) %>%
    mutate(Column = as.numeric(Column))

df.d88.images %<>%
    left_join(df.d88.plate, by = c("Row", "Column"))

all(colnames(df.d54.images) == colnames(df.d88.images))

df.d54.images %>%
    bind_rows(df.d88.images) -> df.images

rm(df.d54.images, df.d88.images)

# summarize by well
df.images %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Day, Primary_AB) %>%
    summarise(Count_Cy5objects = sum(Count_Cy5objects),
              Count_GFPobjects = sum(Count_GFPobjects),
              Count_Nuclei = sum(Count_Nuclei)) -> df.wells

df.wells %<>%
    mutate(Cy5objects_over_nuclei = Count_Cy5objects / Count_Nuclei,
           GFPobjects_over_nuclei = Count_GFPobjects / Count_Nuclei,
           Expression = factor(Expression, levels = c("pTRIPZ-Ctrl", "pTRIPZ-4707"), ordered = TRUE))

pdf("~/Desktop/d54_d88_prolif_ctrl_4707.pdf")

df.wells %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.wells %>%
    ggplot(mapping = aes(y = Count_GFPobjects, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition",
         y = "GFP positive cells")

df.wells %>%
    ggplot(mapping = aes(y = Count_Cy5objects, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition",
         y = "Edu positive nuclei")

df.wells %>%
    ggplot(mapping = aes(y = GFPobjects_over_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition",
         y = "GFP positive / nuclei")

df.wells %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition",
         y = "EdU positive / nuclei")

dev.off()

# Classify Nuclei by GFP and Cy5 Positive ##################
# donor 54 nuclei
df.d54.nuclei <- read_csv(d54.nuclei.csv)

df.d54.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity_")) %>%
    mutate(Column = as.numeric(Column))

df.d54.nuclei %<>%
    left_join(df.d54.plate, by = c("Row", "Column"))

# donor 88 nuclei
df.d88.nuclei <- read_csv(d88.nuclei.csv)

df.d88.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity_")) %>%
    mutate(Column = as.numeric(Column))

df.d88.nuclei %<>%
    left_join(df.d88.plate, by = c("Row", "Column"))

all(colnames(df.d54.nuclei) == colnames(df.d88.nuclei))

df.d54.nuclei %>%
    bind_rows(df.d88.nuclei) -> df.nuclei

rm(df.d54.nuclei, df.d88.nuclei)


df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_Cy5enhance)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = 100)

df.nuclei %<>%
    mutate(EdU_positive = Intensity_IntegratedIntensity_Cy5enhance >= 100)

df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_Cy5enhance, fill = EdU_positive)) +
    geom_histogram(bins = 60) +
    geom_vline(xintercept = 100)

df.nuclei %>%
    ggplot(aes(x = Intensity_IntegratedIntensity_GFPenhance)) +
    geom_histogram(bins = 120) +
    geom_vline(xintercept = 30)

df.nuclei %<>%
    mutate(GFP_positive = Intensity_IntegratedIntensity_GFPenhance >= 20)

df.nuclei %<>%
    mutate(GFP_EdU_positive = GFP_positive & EdU_positive)





# summarize by well
df.nuclei %>%
    group_by(Well, Row, Column, Donor, Expression, DOX, Day, Primary_AB) %>%
    summarise(Count_EdU_positive_nuclei = sum(EdU_positive),
              Count_GFP_positive_nuclei = sum(GFP_positive),
              Count_GFP_EdU_positive_nuclei = sum(GFP_EdU_positive),
              Count_Nuclei = n()) -> df.nuclei.by.well

df.nuclei.by.well %<>%
    mutate(nuclei_over_gfp_positive = Count_Nuclei / Count_GFP_positive_nuclei,
           edu_positive_over_gfp_positive = Count_EdU_positive_nuclei / Count_GFP_positive_nuclei,
           gfp_positive_over_nuclei = Count_GFP_positive_nuclei / Count_Nuclei,
           gfp_edu_positive_over_gfp_positive = Count_GFP_EdU_positive_nuclei / Count_GFP_positive_nuclei,
           Expression = factor(Expression, levels = c("pTRIPZ-Ctrl", "pTRIPZ-4707"), ordered = TRUE))



df.nuclei.by.well %>%
    ggplot(mapping = aes(y = Count_GFP_positive_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = Count_Nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = gfp_positive_over_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = nuclei_over_gfp_positive, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = Count_EdU_positive_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = edu_positive_over_gfp_positive, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = Count_GFP_EdU_positive_nuclei, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

df.nuclei.by.well %>%
    ggplot(mapping = aes(y = gfp_edu_positive_over_gfp_positive, x = Donor, color = Expression)) +
    geom_boxplot() +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    labs(title = "Day 8 Prolif",
         caption = "14 wells/condition")

# Figure ###########################

dir.pdfs <- here("doc/paper/figure6/pdfs/")

# df.wells %>%
#     ggplot(aes(y = Count_Nuclei, x = Donor, color = Expression)) +
#     geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
#     geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
#     plotTheme("figure") +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     scale_color_manual(values = c(paperBlue, paperRed)) +
#     labs(y = "Nuclei Count",
#          x = "Donor") +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# ggsave(paste0(dir.pdfs, "D88_D54_prolif_nuclei_figure.pdf"), height = 2.2, width = 3)
#
# df.wells %>%
#     ggplot(aes(y = Cy5objects_over_nuclei, x = Donor, color = Expression)) +
#     geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
#     geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
#     plotTheme("figure") +
#     theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
#     scale_color_manual(values = c(paperBlue, paperRed)) +
#     labs(y = "EdU Positive / Nuclei",
#          x = "Donor") +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# ggsave(paste0(dir.pdfs, "D88_D54_prolif_edu_figure.pdf"), height = 2.2, width = 3)

df.nuclei.by.well %>%
    ggplot(aes(y = Count_GFP_positive_nuclei, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3, outlier.shape = NA) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c(paperBlue, paperRed)) +
    labs(y = "GFP Positive Nuclei Count",
         x = "Donor") +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D88_D54_prolif_gfp_nuclei_figure.pdf"), height = 2.2, width = 3)

df.nuclei.by.well %>%
    ggplot(aes(y = gfp_edu_positive_over_gfp_positive, x = Donor, color = Expression)) +
    geom_boxplot(position = position_dodge(width = 0.85), size = 0.3) +
    geom_point(position = position_jitterdodge(dodge.width = 0.85, jitter.width = 0.1), size = 0.3) +
    plotTheme("figure") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c(paperBlue, paperRed)) +
    labs(y = "EdU and GFP Positive / GFP Positive Nuclei",
         x = "Donor") +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

ggsave(paste0(dir.pdfs, "D88_D54_prolif_edu_gfp_figure.pdf"), height = 2.2, width = 3)



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
# # GFP objects
# df.w1.GFPobjects <- read_csv(week1.GFPobjects.csv)
#
# df.w1.GFPobjects %<>%
#     select(ImageNumber,
#            ObjectNumber,
#            Well = Metadata_Well,
#            Row = Metadata_Row,
#            Column = Metadata_Column,
#            starts_with("Intensity")) %>%
#     mutate(Column = as.numeric(Column))
#
# df.w1.GFPobjects %<>%
#     left_join(df.w1.plate, by = c("Row", "Column"))
#
# df.w2.GFPobjects <- read_csv(week2.GFPobjects.csv)
#
# df.w2.GFPobjects %<>%
#     select(ImageNumber,
#            ObjectNumber,
#            Well = Metadata_Well,
#            Row = Metadata_Row,
#            Column = Metadata_Column,
#            starts_with("Intensity")) %>%
#     mutate(Column = as.numeric(Column))
#
# df.w2.GFPobjects %<>%
#     left_join(df.w2.plate, by = c("Row", "Column"))
#
# all(colnames(df.w1.GFPobjects) == colnames(df.w2.GFPobjects))
#
# df.w1.GFPobjects %>%
#     bind_rows(df.w2.GFPobjects) -> df.GFPobjects
#
# rm(df.w1.GFPobjects, df.w2.GFPobjects)
#
#
# df.wells %<>%
#     filter(! (Week == 2 & Well == "E11")) %>%
#     filter(! (Week == 2 & Well == "E09")) %>%
#     filter(! (Week == 2 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "E03")) %>%
#     filter(! (Week == 1 & Well == "B02")) %>%
#     filter(! (Week == 1 & Well == "E11")) %>%
#     filter(! (Week == 1 & Well == "B11")) %>%
#     filter(! (Week == 1 & Well == "B09")) %>%
#     filter(! (Week == 1 & Well == "E10")) %>%
#     filter(! (Week == 1 & Well == "B04")) %>%
#     filter(! (Week == 1 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "B08"))
#
# df.GFPobjects %<>%
#     filter(! (Week == 2 & Well == "E11")) %>%
#     filter(! (Week == 2 & Well == "E09")) %>%
#     filter(! (Week == 2 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "E03")) %>%
#     filter(! (Week == 1 & Well == "B02")) %>%
#     filter(! (Week == 1 & Well == "E11")) %>%
#     filter(! (Week == 1 & Well == "B11")) %>%
#     filter(! (Week == 1 & Well == "B09")) %>%
#     filter(! (Week == 1 & Well == "E10")) %>%
#     filter(! (Week == 1 & Well == "B04")) %>%
#     filter(! (Week == 1 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "B08"))
#
# df.Cy5objects %<>%
#     filter(! (Week == 2 & Well == "E11")) %>%
#     filter(! (Week == 2 & Well == "E09")) %>%
#     filter(! (Week == 2 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "E03")) %>%
#     filter(! (Week == 1 & Well == "B02")) %>%
#     filter(! (Week == 1 & Well == "E11")) %>%
#     filter(! (Week == 1 & Well == "B11")) %>%
#     filter(! (Week == 1 & Well == "B09")) %>%
#     filter(! (Week == 1 & Well == "E10")) %>%
#     filter(! (Week == 1 & Well == "B04")) %>%
#     filter(! (Week == 1 & Well == "E02")) %>%
#     filter(! (Week == 1 & Well == "B08"))
#
# df.GFPobjects %<>%
#     filter(! (Expression == "Untransduced"))
#
# df.GFPobjects %>%
#     ggplot(aes(x = Intensity_MaxIntensity_GFP)) +
#     geom_histogram(bins = 50)
#
# df.GFPobjects %>%
#     ggplot(aes(x = Intensity_MaxIntensity_Cy5)) +
#     geom_histogram(bins = 50)
#
# df.GFPobjects %<>%
#     mutate(Cy5positive_GFPobjects = Intensity_MaxIntensity_Cy5 > 0.8)
#
# df.Cy5objects %>%
#     group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
#     summarise(Count_Cy5objects = n()) -> df.Cy5objects.summary
#
# df.GFPobjects %>%
#     group_by(Well, Row, Column, Donor, Expression, DOX, Week, Primary_AB) %>%
#     summarise(Count_GFPobjects = n(),
#               Count_Cy5positive_GFPobjects = sum(Cy5positive_GFPobjects)) -> df.GFPobjects.summary
#
# df.Cy5objects.summary %<>%
#     mutate(Week = factor(Week))
#
# df.GFPobjects.summary %<>%
#     mutate(Week = factor(Week))
#
# df.Cy5objects.summary %<>%
#     left_join(df.wells, by = c("Well", "Row", "Column", "Donor", "Expression", "DOX", "Week", "Primary_AB", "Count_Cy5objects"))
#
# df.GFPobjects.summary %<>%
#     left_join(df.wells, by = c("Well", "Row", "Column", "Donor", "Expression", "DOX", "Week", "Primary_AB", "Count_GFPobjects"))
#
# df.GFPobjects.summary %<>%
#     mutate(Cy5positive_GFPobjects_over_nuclei = Count_Cy5positive_GFPobjects / Count_Nuclei)
#
# df.GFPobjects.summary %<>%
#     mutate(Cy5positive_GFPobjects_over_GFPobjects = Count_Cy5positive_GFPobjects / Count_GFPobjects)
#
#
#
#
# pdf("~/Desktop/d54_d88_diff_ctrl_4707.pdf")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray"))
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Count_Nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Count_GFPobjects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = GFPobjects_over_nuclei, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_GFPobjects_over_GFPobjects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.GFPobjects.summary %>%
#     ggplot(mapping = aes(y = Count_Cy5positive_GFPobjects, x = Week, color = Expression)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange", "gray")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
#
#
#
#
#
#
#
#
#
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Count_Cy5objects, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.5, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.8, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.9, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.5_over_nuclei, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.8_over_nuclei, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# df.Cy5objects.summary %>%
#     ggplot(mapping = aes(y = Cy5positive_.9_over_nuclei, x = Week, color = DOX)) +
#     geom_boxplot() +
#     facet_wrap(~Primary_AB + Donor) +
#     scale_color_manual(values = c("darkblue", "darkorange")) +
#     stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")
#
# dev.off()
