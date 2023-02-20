
# create file to update sample IID and FID among genotype samples

library(readr)
library(dplyr)
library(magrittr)


fam.old <- read_delim("~/Projects/BACKUP_EXCLUDED/updated.fam.ids", delim = " ", col_names = FALSE)

fam.dan <- read_delim("~/Projects/BACKUP_EXCLUDED/allmgergednew.fam", delim = "\t", col_names = FALSE)

fam.dan %<>%
    select(fid = X1, iid = X2)

combo <- left_join(fam.old, fam.dan, by = c("X3" = "iid"))

combo %<>%
    select(X1, X2, fid, X4)

write_delim(combo, "~/Projects/BACKUP_EXCLUDED/updated.fam.ids.new", delim = " ", col_names = FALSE)
