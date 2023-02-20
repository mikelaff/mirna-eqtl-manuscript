# create rsid conversion file

# create large list of variant id conversions to rsid from all arrays used for genotyping
# start with most recent conversion list from illumina and work backwards

library(readr)
library(dplyr)
library(magrittr)


set10 <- read_tsv("~/Projects/BACKUP_EXCLUDED/infinium-omni2-5exome-8-v1-4-a1-b150-rsids.zip")
set9 <- read_tsv("~/Projects/BACKUP_EXCLUDED/infinium-omni2-5-exome-8-v1-3-a1-b144-rsids.zip")
set3 <- read_tsv("~/Projects/BACKUP_EXCLUDED/humanomni2-5exome-8v1-1-loci-name-to-rsid-conversion-file.zip")
set1 <- read_tsv("~/Projects/BACKUP_EXCLUDED/HumanOmni2-5-8-v1-2-A-b138-rsIDs.zip")




combo <- bind_rows(set10, filter(set9, ! Name %in% set10$Name))

combo <- bind_rows(combo, filter(set3, ! Name %in% combo$Name))

combo <- bind_rows(combo, filter(set1, ! Name %in% combo$Name))

# remove rows with no rsid
combo %<>%
    filter(RsID != ".")


sum(duplicated(combo$Name))
sum(duplicated(combo$RsID))

write_tsv(combo, "~/Projects/BACKUP_EXCLUDED/rsid_conversion_file.tsv.gz")





