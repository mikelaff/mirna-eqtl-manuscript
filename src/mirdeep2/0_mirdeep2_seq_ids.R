# script to create 3 character IDs for fastq file names used in miRDeep2

library(tidyverse)

# read in fastq.gz file names and label column
files <- read_csv("fastq_files.txt", col_names = FALSE)
files <- rename(files, gz_file = X1)

# make column for RNAIDs
files$rnaid <- sapply(strsplit(files$gz_file, "_"), `[`,1)

# extract RNAIDs into separate df
rnaids <- distinct(select(files, rnaid))

# new column for seqletters
rnaids$seqletters <- NA

# index for letters
let1 <- 1
let2 <- 1

for(i in 1:length(rnaids$rnaid)) {
  # assing letters by index into letters
  letter1 <- letters[let1]
  letter2 <- letters[let2]
  
  # update let2
  let2 <- let2 + 1
  # check for let2 rollover
  if(let2 > 26) {
    # reset let2
    let2 <- 1
    # advance let1
    let1 <- let1 + 1
  }
  # check for let1 rollover
  if(let1 > 26) {
    # reset let1
    let1 <- 1
  }
  
  # assign letters to each RNAID
  rnaids$seqletters[i] <- paste(letter1, letter2, sep="")
}

# join rnaids by rnaid with files
files <- left_join(files, rnaids, by = "rnaid")

# new column for seqnumber
files$seqnumber <- NA

# seqnumber iterator
num <- 1

# list of duplicate RNAIDs
dupIDs <- duplicated(files$rnaid)

# loop over all rnaids
for(i in 1:length(files$rnaid)) {
  # if duplicated rnaid is false, reset num
  if(dupIDs[i] == FALSE) {
    num <- 1
  }
  # assign number of seqnumber
  files$seqnumber[i] <- num
  # advance seqnumber
  num <- num + 1
}

# combine seqletters and seqnumber into seqid
files$seqid <- paste(files$seqletters, files$seqnumber, sep="")

# fastq file names
files$fastq_file <- paste(sapply(strsplit(files$gz_file, "\\."), "[", 1), ".fastq", sep="")

# fastq file paths
files$fastq_path <- paste("fastq/", files$fastq_file, sep="")

# output fastq paths and seqids
write_tsv(select(files, fastq_path, seqid), "config.txt", col_names = FALSE)
