# compile mirna count matrices from featureCounts against small-rna-seq bams
# all_known_and_novel
# mirdeep_mirge
# mirdeep_mirge_friedlander_nowakowski

library(here)
library(readr)
library(dplyr)
library(magrittr)

# INPUT FILES #########################################################################################################
# counts directory roots
root.dir <- "~/scr/small_rna_mapping/featureCounts/"

# samples
samples.txt <- here("src/mapping_bowtie_to_mature_miRNA/RNAIDnumbers.txt")

# OUTPUT FILES ########################################################################################################
dir.output <- here("results/counts/small_rna_seq/bowtie_mir4707/")

output.tsv <- paste(dir.output, format(Sys.time(), "%Y%m%d"),
                    "_mir4707_counts.tsv", sep="")

# Compile Counts #################################################################################
# samples
samples <- read_lines(samples.txt)

# load first count file
countData <- read_tsv(paste0(root.dir, samples[1], "/", samples[1], ".mir4707.counts"),
                      skip = 2,
                      col_names = c("name", "chr", "start", "end", "strand", "length", samples[1]))

#loop through all samples and add count data to countData object
for (i in 2:length(samples)) {
    #load count file
    tempData <- read_tsv(paste0(root.dir, samples[i], "/", samples[i], ".mir4707.counts"),
                         skip = 2,
                         col_names = c("name", "chr", "start", "end", "strand", "length", samples[i]))

    # select only ID and count column
    tempData %<>% dplyr::select(1,7)
    #join tempData with countData
    countData <- left_join(countData, tempData, by="name")
}

#remove objects
rm(tempData, i)

countData %>%
    select(name, starts_with("RNAID")) -> tmp


tmp <- data.frame(tmp[,-1], row.names = tmp$name)
tmp <- as.data.frame(t(tmp))

# Export Count Matrices ###############################################################################################

write_tsv(countData, output.tsv)
