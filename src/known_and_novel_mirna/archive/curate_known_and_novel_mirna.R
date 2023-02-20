# curate list of known and novel mirnas from
# mirbase v 22
# friedlander 2014 paper
# mirdeep2 novel (dec 2018)
# mirge2.0 novel (dec 2018)

library(here)
library(readr)
library(readxl)
library(dplyr)
library(magrittr)
library(rtracklayer)
library(Biostrings)
library(stringr)

# OUTPUT ###########################
#known_and_novel_gff <- here("results/known_and_novel_mirna/20190130_known_and_novel.gff3")
known_and_novel_hairpin_fa <- here("results/known_and_novel_mirna/20190208_known_and_novel_hairpin.fa")
known_and_novel_mature_fa <- here("results/known_and_novel_mirna/20190208_known_and_novel_mature.fa")
# novel mirnas compiled for all sources in table with more information
novel_mirnas_tsv <- here("results/known_and_novel_mirna/20190214_novel_mirnas.tsv")

# INPUT ############################
# mirbase v22 files
mirbase_hsa_hairpin_fa <- here("data/mirbase22/hsa_hairpin.fa")
mirbase_hsa_mature_fa <- here("data/mirbase22/hsa_mature.fa")
mirbase_hsa_gff <- here("data/mirbase22/hsa.gff3")
# friedlander 2014 paper
friedlander_xls <- here("data/friedlander_novel_2014/13059_2013_3254_MOESM3_ESM.xls")
# mirdeep2 novel
mirdeep2_novel_csv <- here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.csv")
# mirge2.0 novel
mirge_novel_tsv <- here("results/mirge2.0/20190130_mirge2.0_novel_miRNAs.tsv")


# Import mirdeep2 novel ######################

df_mirdeep <- read_csv(mirdeep2_novel_csv)
df_mirdeep$id <- 1:length(df_mirdeep$provisional_id)

# mature sequences will have 5p
df_mirdeep %>%
  dplyr::select(id, seq = consensus_mature_sequence) %>%
  dplyr::mutate(id = paste("hsa-mirD-", id, "-5p", sep="")) -> mirdeep_5p

# star sequences will have 3p
df_mirdeep %>%
  dplyr::select(id, seq = consensus_mature_sequence) %>%
  dplyr::mutate(id = paste("hsa-mirD-", id, "-3p", sep="")) -> mirdeep_3p

df_mirdeep %>%
  dplyr::select(id, seq = consensus_mature_sequence) %>%
  dplyr::mutate(id = paste("hsa-miD-", id, sep="")) -> mirdeep_novel_precursor

# mirdeep novel, combination of 5p and 3p identified mirnas
mirdeep_novel_mature <- bind_rows(mirdeep_5p, mirdeep_3p)

rm(mirdeep_3p, mirdeep_5p)


# Import mirge2.0 novel #####################

df_mirge <- read_tsv(mirge_novel_tsv)

df_mirge %>%
  dplyr::filter(!duplicated(`Mature miRNA sequence`)) %>%
  dplyr::mutate(arm = ifelse(`Arm type` == "arm5", "5p", "3p")) %>%
  dplyr::select(id = provisional_id, seq = `Mature miRNA sequence`, arm) %>%
  dplyr::mutate(id = paste("hsa-mirG-", id, "-", arm, sep="")) %>%
  dplyr::select(-arm) -> mirge_novel_mature

df_mirge %>%
  dplyr::filter(!duplicated(`Mature miRNA sequence`)) %>%
  dplyr::select(id = provisional_id, seq = `Precursor miRNA sequence`) %>%
  dplyr::mutate(id = paste("hsa-mirG-", id, sep="")) -> mirge_novel_precursor

rm()

# Import friedlander novel ###################

df_fried <- read_excel(friedlander_xls)

df_fried %>%
  dplyr::filter(!is.na(`Confidence leve`)) %>%
  dplyr::select(id = Identifier, seq = `Mature sequence`) %>%
  dplyr::mutate(id = paste("hsa-can-", sapply(strsplit(id, "_"), `[`, 3), sep="")) -> fried_novel_mature

df_fried %>%
  dplyr::filter(!is.na(`Confidence leve`)) %>%
  dplyr::select(id = Identifier, seq = `Hairpin sequence`) %>%
  dplyr::mutate(id = paste("hsa-can-", sapply(strsplit(id, "_"), `[`, 3), sep="")) -> fried_novel_precursor

rm()

# Import mirbase files ###################
# import gff3 sequence information
df_mirbase <- as.data.frame(readGFF(mirbase_hsa_gff))

mirbase_mature <- readRNAStringSet(mirbase_hsa_mature_fa)
mirbase_mature <- as.data.frame(mirbase_mature)

mirbase_mature$id <- rownames(mirbase_mature)

mirbase_mature <- tibble(id = mirbase_mature$id, seq = mirbase_mature$x)

mirbase_precursor <- readRNAStringSet(mirbase_hsa_hairpin_fa)
mirbase_precursor <- as.data.frame(mirbase_precursor)

mirbase_precursor$id <- rownames(mirbase_precursor)

mirbase_precursor <- tibble(id = mirbase_precursor$id, seq = mirbase_precursor$x)

# Combine mirna data frames and compare #######################
# dplyr::bind_rows(dplyr::mutate(mirbase_mature, source = "mirbase"),
#                  dplyr::mutate(fried_novel_mature, source = "friedlander"),
#                  dplyr::mutate(mirdeep_novel_mature, source = "mirdeep"),
#                  dplyr::mutate(mirge_novel_mature, source = "mirge")) -> all_mature
# 
# dplyr::bind_rows(dplyr::mutate(mirbase_precursor, source = "mirbase"),
#                  dplyr::mutate(fried_novel_precursor, source = "friedlander"),
#                  dplyr::mutate(mirdeep_novel_precursor, source = "mirdeep"),
#                  dplyr::mutate(mirge_novel_precursor, source = "mirge")) -> all_precursor
# 
# all_mature$seq <- toupper(all_mature$seq)
# all_precursor$seq <- toupper(all_precursor$seq)
# 
# all_mature$dup <- duplicated(all_mature$seq) #| duplicated(all_mature$seq, fromLast=TRUE)
# all_precursor$dup <- duplicated(all_precursor$seq) #| duplicated(all_precursor$seq, fromLast=TRUE)

# Combine data frames and write fasta files #########################
# dplyr::bind_rows(mirbase_mature,
#                  fried_novel_mature,
#                  mirdeep_novel_mature,
#                  mirge_novel_mature) -> all_mature
# 
# dplyr::bind_rows(mirbase_precursor,
#                  fried_novel_precursor,
#                  mirdeep_novel_precursor,
#                  mirge_novel_precursor) -> all_precursor
# 
# all_mature$seq <- toupper(all_mature$seq)
# all_precursor$seq <- toupper(all_precursor$seq)
# 
# #remove duplicates
# all_mature <- all_mature[!duplicated(all_mature$seq),]
# all_precursor <- all_precursor[!duplicated(all_precursor$seq),]
# 
# seq_mature <- all_mature$seq
# names(seq_mature) <- all_mature$id
# mature_set <- RNAStringSet(seq_mature)
# 
# seq_precursor <- all_precursor$seq
# names(seq_precursor) <- all_precursor$id
# precursor_set <- RNAStringSet(seq_precursor)
# 
# writeXStringSet(mature_set, known_and_novel_mature_fa)
# writeXStringSet(precursor_set, known_and_novel_hairpin_fa)

# Combine all novel into data frame for export ################

df_mirge %<>% 
  dplyr::filter(!duplicated(`Mature miRNA sequence`)) %>%
  dplyr::select(chromosome = Chr, start = `Start Pos`, end = `End Pos`, strand = Strand, id = provisional_id) %>%
  dplyr::mutate(width = end - start + 1) %>%
  dplyr::mutate(id = paste("mirge_", id, sep=""), source = "mirge")

df_mirdeep$seqid <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 1)
df_mirdeep$strand <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 3)
df_mirdeep$pos <- sapply(strsplit(df_mirdeep$precursor_coordinate, ":"), `[`, 2)
df_mirdeep$start <- as.integer(sapply(strsplit(df_mirdeep$pos, "[.]"), `[`, 1))
df_mirdeep$end <- as.integer(sapply(strsplit(df_mirdeep$pos, "[.]"), `[`, 3))

df_mirdeep %<>%
  dplyr::select(chromosome = seqid, start, end, strand, id = provisional_id) %>%
  dplyr::mutate(width = end - start + 1) %>%
  dplyr::mutate(id = paste("mirdp_", id, sep=""), source = "mirdeep")

df_all_novel <- bind_rows(df_mirge, df_mirdeep)
df_all_novel$width <- as.integer(df_all_novel$width)

write_tsv(df_all_novel, novel_mirnas_tsv)



































