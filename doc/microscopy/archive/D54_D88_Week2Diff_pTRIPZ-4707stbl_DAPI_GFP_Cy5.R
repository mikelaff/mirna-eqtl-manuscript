
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
# D54/D88 Week 2 Diff, pTRIPZ-4707 stable line, +/- DOX, DAPI/GFP/Cy5 (Tuj1 or GFAP)
nuclei.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Nuclei.csv")
Cy5objects.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Cy5objects.csv")
GFPobjects.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/")
images.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/output/D54_D88_Diff2week_pTRIPZ-4707stbl_DAPI_GFP_Cy5_Image.csv")
plate.layout.csv <- here("results/microscopy/D54_D88_Diff2week_pTRIPZ-4707stbl_GFAP_Tuj1/plate_layout.csv")

# GLOBALS ##############################################################################################################


# Import Data ##########################################################################################################

# plate layout metadata
df.plate <- read_csv(plate.layout.csv)


# images
df.images <- read_csv(images.csv)

df.images %<>%
    select(ImageNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           Count_Cy5objects,
           Count_Nuclei) %>%
    mutate(Column = as.numeric(Column))

df.images %<>%
    left_join(df.plate, by = c("Row", "Column"))

df.images %<>%
    mutate(Cy5objects_over_nuclei = Count_Cy5objects / Count_Nuclei)

df.images %>%
    ggplot(mapping = aes(y = Cy5objects_over_nuclei, x = DOX, color = Donor)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB)


# Cy5 objects
df.Cy5objects <- read_csv(Cy5objects.csv)

df.Cy5objects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.Cy5objects %<>%
    left_join(df.plate, by = c("Row", "Column"))

# only cy5 objects with a max intensity > 0.8
df.Cy5objects %<>%
    mutate(Cy5positive = Intensity_MaxIntensity_Cy5 > 0.9)

df.Cy5objects %>%
    group_by(ImageNumber, Donor, Expression, DOX, Primary_AB) %>%
    summarise(Cy5objects_num = n(),
              Cy5positive_num = sum(Cy5positive)) -> df.Cy5objects.summary


df.Cy5objects.summary %<>%
    left_join(df.images, by = c("ImageNumber", "Donor", "Expression", "DOX", "Primary_AB"))



df.Cy5objects.summary %<>%
    mutate(Cy5positive_over_nuclei = Cy5positive_num / Count_Nuclei)

pdf("~/Desktop/week2diff_4707stable.pdf")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(x = Donor, y = Cy5positive_over_nuclei, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB) +
    labs(title = "Week 2 Diff. pTRIPZ-4707stable",
         caption = "48 wells, 192 images") +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(x = Donor, y = Cy5positive_num, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB) +
    labs(title = "Week 2 Diff. pTRIPZ-4707stable",
         caption = "48 wells, 192 images") +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

df.Cy5objects.summary %>%
    ggplot(mapping = aes(x = Donor, y = Count_Nuclei, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB) +
    labs(title = "Week 2 Diff. pTRIPZ-4707stable",
         caption = "48 wells, 192 images") +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format")

dev.off()


stop()



























# nuclei
df.nuclei <- read_csv(nuclei.csv)

df.nuclei %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.nuclei %<>%
    left_join(df.plate, by = c("Row", "Column"))



df.nuclei %>%
    #filter(Expression == "Control") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_GFPenhance)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Donor + DOX, scales = "free")

df.nuclei %>%
    #filter(Expression == "Control") %>%
    ggplot(mapping = aes(x = Intensity_IntegratedIntensity_Cy5enhance)) +
    geom_histogram(bins = 50) +
    facet_wrap(~Primary_AB + Donor + DOX, scales = "fixed", nrow = 2)

df.nuclei %<>%
    mutate(GFP_positive = Intensity_IntegratedIntensity_GFPenhance > 25,
           Cy5_positive = Intensity_IntegratedIntensity_Cy5enhance > 25)

df.nuclei %>%
    group_by(ImageNumber, Donor, Expression, DOX, Primary_AB) %>%
    summarise(num_nuclei = n(),
              GFP_positive_num = sum(GFP_positive),
              Cy5_positive_num = sum(Cy5_positive)) %>%
    mutate(GFP_over_nuclei = GFP_positive_num / num_nuclei,
           Cy5_over_nuclei = Cy5_positive_num / num_nuclei) -> df.nuclei.summary

df.nuclei.summary %>%
    ggplot(mapping = aes(x = DOX, y = GFP_over_nuclei, color = Donor)) +
    geom_boxplot()

library(ggpubr)
df.nuclei.summary %>%
    ggplot(mapping = aes(x = Donor, y = Cy5_over_nuclei, color = DOX)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB) +
    stat_compare_means(label.x=2, label.y=.25, size = 3, method = "t.test", label = "p.format") +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    plotTheme("presentation") +
    labs(y = "Cy5 Positive Cells / Number of Nuclei")


ggsave(file = "~/Desktop/week2diff_cy5_cells.pdf")

write_csv(df.nuclei.summary, file = "desktop.csv")



# GFP objects
df.GFPobjects <- read_csv(GFPobjects.csv)

df.GFPobjects %<>%
    select(ImageNumber,
           ObjectNumber,
           Well = Metadata_Well,
           Row = Metadata_Row,
           Column = Metadata_Column,
           starts_with("Intensity")) %>%
    mutate(Column = as.numeric(Column))

df.GFPobjects %<>%
    left_join(df.plate, by = c("Row", "Column"))





# images
df.images <- read_csv(images.csv)






# Plots ################################################################################################################

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = num_Cy5objects_over_nuclei, x = DOX, color = Donor)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB)

df.Cy5objects.summary %>%
    ggplot(mapping = aes(y = Cy5objects_avg_intensity, x = DOX, color = Donor)) +
    geom_boxplot() +
    facet_wrap(~Primary_AB)





stop()


































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




