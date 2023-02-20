
library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/qpcr/pdfs/")
dir.create(dir.pdfs, recursive = TRUE, showWarnings = FALSE)

# INPUT ################################################################################################################
data.csv <- here("results/qpcr/20210122_HNP_miRNA.xlsx")


# GLOBALS ##############################################################################################################


# Import mRNA Data #####################################################################################################


df.data <- read_xlsx(data.csv, sheet = 3, range = "A47:O95", na = c("", "Undetermined"))

df.data %<>%
    select(Sample = `Sample Name`,
           Target = `Target Name`,
           Well = `Well Position`,
           CT)

df.data %<>%
    mutate(Donor = sapply(strsplit(Sample, "_"), `[`, 1),
           Expression = sapply(strsplit(Sample, "_"), `[`, 2),
           Virus = sapply(strsplit(Sample, "_"), `[`, 3))

df.data %<>%
    select(Sample,
           Donor,
           Expression,
           Virus,
           Well,
           Target,
           CT)

# df.data %<>%
#     filter(! Well == "P1",
#            ! Well == "P2",
#            ! Well == "P3")

df.data$Expression <- factor(df.data$Expression,
                             levels = c("mock", "haus4", "4707-C", "4707-A"),
                             labels = c("Control", "Lenti-HAUS4", "Lenti-miR-4707-3p-REF(C)", "Lenti-miR-4707-3p-ALT(A)"),
                             ordered = TRUE)

df.data$Donor <- factor(df.data$Donor,
                        levels = c("D52", "D111"))

df.data %<>%
    mutate(Name = paste(Donor, Expression, Virus))

# Calculate mean CT values across duplicates, only retain one row per sample/target pair
df.data %<>%
    group_by(Sample, Target) %>%
    mutate(CT_mean = mean(CT),
              CT_sd = sd(CT)) %>%
    select(-CT, -Well) %>%
    distinct()



# loop over samples, calculate delta CT to miR-361
samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL

    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp

    # get miR-361-5p CT value for this sample
    ct.361 <- df.tmp$CT_mean[match("miR-361-5p", df.tmp$Target)]

    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_361 = CT_mean - ct.361)

    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)

}

df.data <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Mock of that donor for each target
donors <- unique(df.data$Donor)
targets <- unique(df.data$Target)

df.new <- tibble()

for (donor in donors) {
    df.tmp <- NULL

    # filter for only this donor
    df.data %>%
        filter(Donor == donor) -> df.tmp

    # loop over targets
    for (target in targets) {
        df.tmp2 <- NULL

        # filter for only this target
        df.tmp %>%
            filter(Target == target) -> df.tmp2

        delta.361.control <- df.tmp2$delta_CT_361[match("Control", df.tmp2$Expression)]

        # calculate delta delta CT within this target and donor
        df.tmp2 %<>%
            mutate(delta_delta_CT_361 = delta_CT_361 - delta.361.control)

        # calculate fold change
        df.tmp2 %<>%
            mutate(fold_change_361 = 2 ^ - (delta_delta_CT_361))

        # combine into all data
        df.new <- bind_rows(df.new, df.tmp2)
    }

}

df.data <- df.new
rm(df.new, df.tmp, df.tmp2)

pdf(paste0(dir.pdfs, "20200122_HNP_expression.pdf"))

df.data %>%
    ggplot(aes(y = fold_change_361, x = Expression, fill = Donor)) +
    geom_col(position = "dodge") +
    facet_wrap(~Target) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 5 days, virus=50uL",
         title = "qPCR")

df.data %>%
    ggplot(aes(y = CT_mean, x = Expression, color = Donor)) +
    geom_point(position = position_dodge(width = 0.85)) +
    facet_wrap(~Target) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_color_manual(values = c("darkblue", "darkorange")) +
    labs(y = "mean(CT)",
         x = "Expression Vector",
         caption = "HNP expression after 5 days, virus=50uL",
         title = "qPCR")

dev.off()
#ggsave("20210122_miRNA_qPCR_miR361.pdf", height = 8, width = 8)

