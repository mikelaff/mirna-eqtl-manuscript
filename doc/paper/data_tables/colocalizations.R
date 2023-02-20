
library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
#library(DESeq2)

# OUTPUT FILES #########################################################################################################
# supplementary table 4: colocalizations
output.blood.csv <- here("doc/paper/data_tables/csv/supplementaryTable4_blood.csv")
output.mrna.eqtl.csv <- here("doc/paper/data_tables/csv/supplementaryTable4_mrna-eqtl.csv")
output.mrna.sqtl.csv <- here("doc/paper/data_tables/csv/supplementaryTable4_mrna-sqtl.csv")

# INPUT FILES ##########################################################################################################
# mrna-eqtl colocalizations
eqtl.rds <- here("results/co-localization/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTL_mirQTL_overlaps_r2at0.8.rds")
# mrna-sqtl colocalizations
sqtl.rds <- here("results/co-localization/fetal_brain_sQTL/fetal_brain_sQTL_mirQTL_overlaps_r2at0.8.rds")
# blood colocalizations
blood.rds <- here("results/co-localization/blood_miRNA-eQTL/blood_miRNA-eQTL_mirQTL_overlaps_r2at0.8.rds")

# GLOBALS ##############################################################################################################

# Import #####

df.mrna.eqtl <- read_rds(eqtl.rds)

df.mrna.sqtl <- read_rds(sqtl.rds)

df.blood <- read_rds(blood.rds)

# format columns



write_csv(df.mrna.eqtl, output.mrna.eqtl.csv)
write_csv(df.mrna.sqtl, output.mrna.sqtl.csv)
write_csv(df.blood, output.blood.csv)




