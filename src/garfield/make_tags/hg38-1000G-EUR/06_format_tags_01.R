# convert a plink .tags.list to files suitable to garfield tags

library(here)
library(readr)
library(dplyr)
library(magrittr)

# OUTPUT FILES #########################################################################################################
output.dir <- here("data/garfield-data/tags/hg38-1000G-EUR/r01/")

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
input.dir <- "/pine/scr/m/i/mikelaff/garfield/"

#input.dir <- "~/Downloads/"

# GLOBALS ##############################################################################################################
CHROMS <- c(seq(1,22), "X", "Y")

# Import/Export ########################################################################################################
# Loop over each chr
for (chrom in CHROMS) {

    df.tags <- NULL

    print(paste("Working on chrom:", chrom))

    # import .tags.list
    print("Importing...")
    df.tags <- read_table2(paste0(input.dir, "chr", chrom, "_01.tags.list"), col_types = "cciiiidc")

    # filter for at least 1 tag, remove excess columns
    print("Filtering...")
    df.tags %>%
        dplyr::filter(NTAG > 0) %>%
        dplyr::select(BP, TAGS) %>%
        dplyr::filter(!duplicated(BP)) -> df.tags

    # remove chr and :
    print("Formatting...")
    df.tags$TAGS <- gsub(paste0(chrom, ":"), "", df.tags$TAGS)

    # replace | with ,
    df.tags$TAGS <- gsub("\\|", ",", df.tags$TAGS)

    # export formated tags file
    print("Exporting...")
    write_delim(df.tags, paste0(output.dir, "chr", chrom), delim = " ", col_names = FALSE)

}
