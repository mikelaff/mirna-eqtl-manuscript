# to import known and novel miRNA data and create GRanges object

# miRBase release 21 homo sapiens data

library(here)
library(rtracklayer)
library(Biostrings)
library(stringr)
library(dplyr)
library(magrittr)

# mirbase data files
hsa_gff_file <- here("data/mirbase/hsa.gff3")
hsa_mature_file <- here("data/mirbase/hsa_mature.fa")
hsa_precursor_file <- here("data/mirbase/hsa_hairpin.fa")

# mirdeep2 novel mirnas
novel_mirna_file <- here("results/mirdeep2/novel_mirna/novel_mirnas.csv")

# import gff3 sequence information
gr <- import(hsa_gff_file, format = "gff3")
df <- data.frame(gr)

# import and format mature mirna sequences
mirna_seqs <- readRNAStringSet(hsa_mature_file)
mirna_seqs <- data.frame(mirna_seqs)
mirna_seqs$Name <- rownames(mirna_seqs)
mirna_seqs %<>% rename(seq = mirna_seqs)

# import and format precursor sequences
precursor_seqs <- readRNAStringSet(hsa_precursor_file)
precursor_seqs <- data.frame(precursor_seqs)
precursor_seqs$Name <- rownames(precursor_seqs)
precursor_seqs %<>% rename(seq = precursor_seqs)

# bind seqs dfs
seqs <- bind_rows(mirna_seqs, precursor_seqs)
rm(mirna_seqs, precursor_seqs)

# join with ranges df
df %<>% left_join(seqs, by = "Name")
rm(seqs)

# fill in derives from so that primary miRNAs also have derives from info
df$Derives_from <- ifelse(is.na(df$Derives_from), df$ID, df$Derives_from)
df$Derives_from <- factor(df$Derives_from)

# GRanges for precursor and mature miRNA
gr2 <- makeGRangesFromDataFrame(df, keep.extra.columns = TRUE, seqinfo = Seqinfo(genome = "hg38"))

# GRangesList divided by precursor miRNAs
grl <- makeGRangesListFromDataFrame(df,
                                    split.field = c("Derives_from"),
                                    names.field = c("Name"),
                                    keep.extra.columns = TRUE,
                                    seqinfo = Seqinfo(genome = "hg38"))

