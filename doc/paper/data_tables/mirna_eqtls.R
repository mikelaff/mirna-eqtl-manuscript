
library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
#library(DESeq2)

# OUTPUT FILES #########################################################################################################
# supplementary table 2: mirna-eqtls
output.csv <- here("doc/paper/data_tables/csv/supplementaryTable2_mirna-eqtls.csv")


# INPUT FILES ##########################################################################################################
# mirna-eqtls
eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")

# GLOBALS ##############################################################################################################

# Import #####

df.eqtls <- read_rds(eqtls.rds)

# format columns
df.eqtls %>%
    rename(SNP_CHR = SNP.CHR,
           SNP_BP_hg38 = SNP.BP.hg38,
           EFFECT_ALLELE = EFFECT.ALLELE,
           A1_HOM_COUNT = A1.HOM.COUNT,
           HET_COUNT = HET.COUNT,
           A2_HOM_COUNT = A2.HOM.COUNT,
           miR_CHR = miR.CHR,
           miR_START_hg38 = miR.START.hg38,
           miR_END_hg38 = miR.END.hg38,
           miR_WIDTH = miR.WIDTH,
           miR_STRAND = miR.STRAND,
           SEQUENCE = SEQ,
           NOM_P_VALUE_THRESHOLD = NOM.P.VALUE.THRESH) -> df.export

# Export ######

write_csv(df.export, output.csv)




