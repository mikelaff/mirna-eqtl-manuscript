# compile mirge2.0 miRNA novel miRNAs csv files

library(here)
library(readr)
library(dplyr)
library(magrittr)

# INPUT FILES #####################################################################################
# counts directory root
root_dir <- here("results/mirge2.0/novel/miRBase_v22/")
# OUTPUT FILES ####################################################################################
# output file
outputFile <- paste(here("results/mirge2.0/"), format(Sys.time(), "%Y%m%d"), "_mirge2.0_novel_miRNAs.tsv", sep="")

# Compile miRBase v22 Results #####################################################################
# files
files <- list.files(root_dir)

# load first file
novels <- read_csv(paste(root_dir, files[1], sep=""))
novels$rnaid <- strsplit(files[1], "_")[[1]][1]

# loop through each file
for(i in 2:length(files)){
  # read in counts csv
  suppressMessages(
    df <- read_csv(paste(root_dir, files[i], sep=""))
    )
  if(nrow(df) == 0) {
    next
  }
  df$rnaid <- strsplit(files[i], "_")[[1]][1]
  
  # join to novels df
  novels <- bind_rows(novels, df)
}
rm(df, i)

novels$Position <- paste(novels$Chr, novels$`Start Pos`, novels$`End Pos`, novels$Strand, sep = "_")
novels$Unique_by_position <- !duplicated(novels$Position)

novels$provisional_id <- seq(from=1, to=length(novels$`Novel miRNA name`))

#~63 unique novel mirnas
sum(novels$Unique_by_position)

# write novels data table
write_tsv(novels, outputFile, col_names = TRUE)

