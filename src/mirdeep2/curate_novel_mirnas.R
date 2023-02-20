# curate list of novel miRNAs from miRDeep2 novel miRNA identification run
# export list as fasta file to use as input to miRDeep2 quantification run

library(here)
library(rtracklayer)
library(Biostrings)
library(stringr)
library(dplyr)
library(readr)
library(magrittr)

# mirdeep2 novel mirnas
novel_mirna_file <- here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.csv")

# import novel mirnas
novel_mirnas <- read_csv(novel_mirna_file)

# mature sequences will have 5p
novel_mirnas %>%
  dplyr::select(provisional_id, consensus_mature_sequence) %>%
  dplyr::rename(id = provisional_id, seq = consensus_mature_sequence) %>%
  dplyr::mutate(id = paste("hsa-miR-", id, "-5p", sep="")) -> novel_mature

# star sequences will have 3p
novel_mirnas %>%
  dplyr::select(provisional_id, consensus_star_sequence) %>%
  dplyr::rename(id = provisional_id, seq = consensus_star_sequence) %>%
  dplyr::mutate(id = paste("hsa-miR-", id, "-3p", sep="")) -> novel_star

novel_mirnas %>%
  dplyr::select(provisional_id, consensus_precursor_sequence) %>%
  dplyr::rename(id = provisional_id, seq = consensus_precursor_sequence) %>%
  dplyr::mutate(id = paste("hsa-miR-", id, sep="")) -> novel_precursor

# convert to RNAStringSet
seq <- novel_mature$seq
names(seq) <- novel_mature$id
mature_set <- RNAStringSet(seq)

seq <- novel_star$seq
names(seq) <- novel_star$id
star_set <- RNAStringSet(seq)

seq <- novel_mature$seq
names(seq) <- novel_mature$id
mature_set <- RNAStringSet(seq)

seq <- novel_precursor$seq
names(seq) <- novel_precursor$id
precursor_set <- RNAStringSet(seq)

# output to fast file
writeXStringSet(mature_set, here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mature.fa"))
writeXStringSet(star_set, here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_star.fa"))
writeXStringSet(precursor_set, here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_precursor.fa"))

# combine mature and star to be used as input file to mirdeep2
writeXStringSet(mature_set, here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.fa"))
writeXStringSet(star_set, here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.fa"), append = TRUE)

# import novel mirnas
novel_mirnas <- read_csv(novel_mirna_file)

# mature sequences will have 5p
novel_mirnas %>%
  dplyr::select(provisional_id, consensus_mature_sequence, precursor_coordinate) %>%
  dplyr::rename(id = provisional_id, seq = consensus_mature_sequence) %>%
  dplyr::mutate(id = paste("hsa-miR-", id, "-5p", sep="")) -> novel_mature

novel_mature$seqnames <- sapply(strsplit(novel_mature$precursor_coordinate, ":"), `[`, 1)
novel_mature$strand <- sapply(strsplit(novel_mature$precursor_coordinate, ":"), `[`, 3)
novel_mature %<>% dplyr::select(-precursor_coordinate)

# star sequences will have 3p
novel_mirnas %>%
  dplyr::select(provisional_id, consensus_star_sequence) %>%
  dplyr::rename(id = provisional_id, seq = consensus_star_sequence) %>%
  dplyr::mutate(id = paste("hsa-miR-", id, "-3p", sep="")) -> novel_star

novel_mirnas %>%
  dplyr::select(provisional_id, consensus_precursor_sequence, precursor_coordinate) %>%
  dplyr::rename(id = provisional_id, seq = consensus_precursor_sequence) %>%
  dplyr::mutate(id = paste("hsa-mir-", id, sep="")) -> novel_precursor

novel_precursor$seqnames <- sapply(strsplit(novel_precursor$precursor_coordinate, ":"), `[`, 1)
novel_precursor$strand <- sapply(strsplit(novel_precursor$precursor_coordinate, ":"), `[`, 3)
novel_precursor$start <- sapply(strsplit(novel_precursor$precursor_coordinate, ":"), `[`, 2)
novel_precursor$start <- sapply(strsplit(novel_precursor$start, "\\."), `[`, 1)
novel_precursor$end <- sapply(strsplit(novel_precursor$precursor_coordinate, ":"), `[`, 2)
novel_precursor$end <- sapply(strsplit(novel_precursor$end, "\\."), `[`, 3)

novel_precursor %<>% select(-precursor_coordinate)

novel_precursor$source <- NA
novel_precursor$type <- "miRNA_primary_transcript_novel"
novel_precursor$score <- NA
novel_precursor$phase <- NA

novelgr <- makeGRangesFromDataFrame(novel_precursor, keep.extra.columns = TRUE)
