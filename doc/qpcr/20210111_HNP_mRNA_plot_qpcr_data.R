
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(readxl)


data.xls <- "raw/20210111_HNP_mRNA.xls"

df.data <- read_xls(data.xls, sheet = 3, range = "A48:O192")

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



# loop over samples, calculate delta CT to ACTB
samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL
    
    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp
    
    # get ACTB CT value for this sample
    ct.actb <- df.tmp$CT_mean[match("ACTB", df.tmp$Target)]
    
    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_ACTB = CT_mean - ct.actb)
    
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
        
        delta.ct.control <- df.tmp2$delta_CT_ACTB[match("Control", df.tmp2$Expression)]
        
        # calculate delta delta CT within this target and donor
        df.tmp2 %<>%
            mutate(delta_delta_CT_ACTB = delta_CT_ACTB - delta.ct.control)
        
        # calculate fold change
        df.tmp2 %<>%
            mutate(fold_change_ACTB = 2 ^ - (delta_delta_CT_ACTB))
        
        # combine into all data
        df.new <- bind_rows(df.new, df.tmp2)
    }
    
}

df.data <- df.new
rm(df.new, df.tmp, df.tmp2)


df.data %>%
    ggplot(aes(y = fold_change_ACTB, x = Expression, fill = Donor)) +
    geom_col(position = "dodge") +
    facet_wrap(~Target) +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 5 days, virus=50uL",
         title = "qPCR")

ggsave("20210111_mRNA_qPCR_ACTB.pdf", height = 8, width = 8)

# loop over samples, calculate delta CT to EIF4A2
#samples <- unique(df.data$Sample)

df.new <- tibble()

for (sample in samples) {
    df.tmp <- NULL
    
    # filter for only this sample
    df.data %>%
        filter(Sample == sample) -> df.tmp
    
    # get EIF4A2 CT value for this sample
    ct.EIF4A2 <- df.tmp$CT_mean[match("EIF4A2", df.tmp$Target)]
    
    # calculate delta CT within these samples
    df.tmp %<>%
        mutate(delta_CT_EIF4A2 = CT_mean - ct.EIF4A2)
    
    # combine with all data
    df.new <- bind_rows(df.new, df.tmp)
    
}

df.data <- df.new
rm(df.new, df.tmp)

# loop over donors, calculate delta delta CT to Mock of that donor for each target
#donors <- unique(df.data$Donor)
#targets <- unique(df.data$Target)

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
        
        delta.ct.control <- df.tmp2$delta_CT_EIF4A2[match("Control", df.tmp2$Expression)]
        
        # calculate delta delta CT within this target and donor
        df.tmp2 %<>%
            mutate(delta_delta_CT_EIF4A2 = delta_CT_EIF4A2 - delta.ct.control)
        
        # calculate fold change
        df.tmp2 %<>%
            mutate(fold_change_EIF4A2 = 2 ^ - (delta_delta_CT_EIF4A2))
        
        # combine into all data
        df.new <- bind_rows(df.new, df.tmp2)
    }
    
}

df.data <- df.new
rm(df.new, df.tmp, df.tmp2)

df.data %>%
    ggplot(aes(y = fold_change_EIF4A2, x = Expression, fill = Donor)) +
    geom_col(position = "dodge") +
    facet_wrap(~Target, scales = "free_y") +
    geom_hline(yintercept = 1, linetype = "dashed") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
    scale_fill_manual(values = c("darkblue", "darkorange")) +
    labs(y = "Fold Change",
         x = "Expression Vector",
         caption = "HNP expression after 5 days, virus=50uL",
         title = "qPCR")

ggsave("20210111_mRNA_qPCR_EIF4A2.pdf", height = 8, width = 8)
