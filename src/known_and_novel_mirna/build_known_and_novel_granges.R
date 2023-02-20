# import and format novel mirnas from mirge and mirdeep2 for export as GRanges and GTF/GFF files
# import and format mirbase and friedlander mirnas for export as GRanges objects

# export gtf and granges of non-overlapping mirnas with priority of annotation:
# mirbase_v22 > friedlander > nowakowski > mirge > mirdeep2

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(readxl)
library(rtracklayer)
library(liftOver)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(mikelaffr)

datestamp <- format(Sys.time(), "%Y%m%d")

# OUTPUT ##############################################################################################################
# rds files for GRanges objects
all.known.novel.non.overlapping.rds.output.file <- paste0(here("data/gtf_and_granges/"), datestamp, "_all_known_and_novel_mirna_non_overlapping_granges.rds")

# gtf/gff files
# mirdeep.mirge.gtf.output.file <- paste(here("data/gtf_and_granges/"),
#                                        format(Sys.time(), "%Y%m%d"), "_mirdeep_mirge_novel_mirna.gtf", sep="")
# mirdeep.mirge.friedlander.nowakowski.gtf.output.file <- paste(here("data/gtf_and_granges/"),
#                                                               format(Sys.time(), "%Y%m%d"),
#                                                               "_mirdeep_mirge_friedlander_nowakowski_mirna.gtf",
#                                                               sep="")
# all.known.novel.gtf.output.file <- paste(here("data/gtf_and_granges/"),
#                                          format(Sys.time(), "%Y%m%d"), "_all_known_and_novel_mirna.gtf", sep="")
# all.known.novel.non.overlapping.gtf.output.file <- paste(here("data/gtf_and_granges/"),
#                                                          format(Sys.time(), "%Y%m%d"), "_all_known_and_novel_mirna_non_overlapping.gtf", sep="")

# INPUT ################################################################################################################
# miRBase v22 gff file
mirbase.gff <- here("data/mirbase22/hsa.gff3")

# mirge2.0 fasta sequences (with -5p/-3p annotations that are not included in miRBase)
mirge.pseudo.fa <- here("data/mirge2.0/miRge.Libs/human/fasta.Libs/human_mirna_SNP_pseudo_miRBase.fa")

# friedlander 2014 paper
friedlander.xls <- here("data/friedlander_novel_2014/13059_2013_3254_MOESM3_ESM.xls")
# nowakowski 2018 paper
nowakowski.xlsx <- here("data/nowakowski_2018/41593_2018_265_MOESM3_ESM.xlsx")
# mirdeep2 novel
mirdeep2.novel.csv <- here("results/mirdeep2/mirbase_v22/20181212_novel_mirna_prediction/novel_mirnas.csv")
# mirge2.0 novel
mirge.novel.tsv <- here("results/mirge2.0/20190130_mirge2.0_novel_miRNAs.tsv")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import miRBase miRNAs ################################################################################################
print("Importing miRBase v22 miRNAs.")
print("Converting to GRanges object.")
# import mirbase data as genome ranges object
gr.mirbase <- rtracklayer::import(mirbase.gff)
# set seqinfo of gr
seqinfo(gr.mirbase) <- Seqinfo(genome = "hg38")
# set source as miRBase_v22
gr.mirbase$source <- factor("miRBase_v22")

print("Getting RNA sequence")
# get sequence from Hsapiens: convert View to character, convert to DNAStringSet, convert to RNAStringSet
gr.mirbase$sequence <- RNAStringSet(DNAStringSet(as.character(Views(Hsapiens, gr.mirbase))))

# convert alias to character
gr.mirbase$Alias <- as.character(gr.mirbase$Alias)

gr.mirbase

print("Finished importing miRBase v22 miRNAs.")

# Metadata columns
# source (factor): miRBase_v22
# type (factor): miRNA or miRNA_primary_transcript
# score (numeric): NA
# phase (int): NA
# ID (character): unique ID
# Alias (character): often same as ID, not unique
# Name (character): common name
# Derives_from (character): NA for primary_transcript or ID of primary_transcript
# sequence (RNAStringSet)

# Import miRge2.0 Pseudo Sequences #####################################################################################
# miRge2.0 ammended the miRBase annotations to include -5p and -3p annotations when miRBase only specified one

# import fasta sequences
pseudo.fa <- readDNAStringSet(mirge.pseudo.fa)

# expand fasta sequences to include combined miRNAs
combo.indexes <- grep("/", names(pseudo.fa))

for (i in combo.indexes) {

    print(names(pseudo.fa)[i])

    # suffixes to create new names
    suffixes <- strsplit(names(pseudo.fa)[i], "/")[[1]][-1]
    # first in list name to modify the original sequence name
    keep <- strsplit(names(pseudo.fa)[i], "/")[[1]][1]

    # used to create new names
    base <- paste(strsplit(keep, "-")[[1]][1:2], collapse = "-")
    # append suffixes to base to create new names
    newnames <- sapply(suffixes, function(x) paste(base, x, sep = "-"))

    # the sequences to create a new string set
    newSeqs <- rep(as.character(pseudo.fa[[i]]), length(newnames))
    # name sequences with newnames
    names(newSeqs) <- newnames

    # create now DNAStringSet
    newSet <- DNAStringSet(newSeqs)

    # modify name of original string
    names(pseudo.fa)[i] <- keep

    # append newSet to pseudo.fa
    pseudo.fa <- c(pseudo.fa, newSet)

    rm(suffixes, keep, base, newnames, newSeqs, newSet)
}

# sequences not in miRBase
pseudo.fa.new <- pseudo.fa[!names(pseudo.fa) %in% gr.mirbase$Name]

# convert to RNAStringSet
pseudo.fa.new.rna <- RNAStringSet(pseudo.fa.new)


grepl("-5p", names(pseudo.fa.new)) | grepl("-3p", names(pseudo.fa.new))

# loop over sequences not in miRBase and ammend miRBase annotations
# add mirge SNP annotation
# if mirge 5/3p annotation overlaps with a miRBase miRNA annotation, update name only
# if mirge 5/3p annotation overlaps with a miRBase primary but not miRNA annotation, add new annotation derived from miRBase primary
# if not, something is wrong
for (i in 1:length(pseudo.fa.new.rna)) {

    new.RNAstring <- pseudo.fa.new.rna[[i]]
    new.name <- names(pseudo.fa.new.rna)[i]
    new.base <- paste(strsplit(new.name, "-")[[1]][1:3], collapse = "-")
    # there are now let miRNAs in new set, so primary mir should be hsa-mir
    new.pri <- paste0("hsa-mir-", strsplit(new.name, "-")[[1]][3])

    printMessage(paste(i, "of", length(pseudo.fa.new.rna), ":", new.name, "of", new.pri))

    # if base is in miRBase annotations
    if (new.base %in% gr.mirbase$Name) {

        print("base match")

        base.index <- which(gr.mirbase$Name == new.base)
        if (length(base.index) != 1) {
            stop("base.index != 1")
        }

        pri.index <- which(gr.mirbase$Name == new.pri)
        if (length(pri.index) != 1) {
            stop("pri.index != 1")
        }

    } else if (new.pri %in% gr.mirbase$Name) {

        print("pri match")

        pri.index <- which(gr.mirbase$Name == new.pri)
        if (length(pri.index) != 1) {
            stop("pri.index != 1")
        }
    } else {
        printMessage("no mirbase match", width = 50, fillChar = "!")
    }



}

# Import miRDeep2 Novel miRNAs #########################################################################################
print("Importing miRDeep2 putatively novel miRNAs.")
df.mirdeep <- read_csv(mirdeep2.novel.csv)
# determine genome position of precursor sequence
df.mirdeep$seqnames <- sapply(strsplit(df.mirdeep$precursor_coordinate, ":"), `[`, 1)
df.mirdeep$strand <- sapply(strsplit(df.mirdeep$precursor_coordinate, ":"), `[`, 3)
df.mirdeep$pos <- sapply(strsplit(df.mirdeep$precursor_coordinate, ":"), `[`, 2)
# start and end positions are not correct to reference genome?
df.mirdeep$start <- as.integer(sapply(strsplit(df.mirdeep$pos, "[.]"), `[`, 1)) + 1
df.mirdeep$end <- as.integer(sapply(strsplit(df.mirdeep$pos, "[.]"), `[`, 3))

# construct modified mirdeep2 df with relavent informaiton
df.mirdeep %<>%
    dplyr::select(name = provisional_id,
                  score = miRDeep2_score,
                  mature_seq = consensus_mature_sequence,
                  star_seq = consensus_star_sequence,
                  precursor_seq = consensus_precursor_sequence,
                  chr = seqnames,
                  strand,
                  precursor_start = start,
                  precursor_end = end)

df.mirdeep$mature_seq <- toupper(df.mirdeep$mature_seq)
df.mirdeep$star_seq <- toupper(df.mirdeep$star_seq)
df.mirdeep$precursor_seq <- toupper(df.mirdeep$precursor_seq)

df.mirdeep

print("Converting to GRanges object.")

# rna string set
precursor.set <- RNAStringSet(df.mirdeep$precursor_seq)
mature.set <- RNAStringSet(df.mirdeep$mature_seq)
star.set <- RNAStringSet(df.mirdeep$star_seq)

# check precursor length is same as width
# width of sequence should be 1 more than coord. start - coord. end (start and end coords are included in sequence)
if ( ! all(df.mirdeep$precursor_end - df.mirdeep$precursor_start == width(precursor.set) - 1) ) {
    stop("Precursor coordinate length not same as sequence width.")
}

# define precursor width
df.mirdeep$precursor_width <- width(precursor.set)

# variables to find with string matching
df.mirdeep$mature_start <- NA
df.mirdeep$mature_end <- NA
df.mirdeep$mature_width <- NA
df.mirdeep$star_start <- NA
df.mirdeep$star_end <- NA
df.mirdeep$star_width <- NA

# loop over all precursors and define mature and star coordinates
# check coordinates correspond to correct genomic sequence
for (i in 1:length(df.mirdeep$name)) {
    # status message
    print(paste("Working on novel ", i, " of ", length(precursor.set), ": ",
                df.mirdeep$name[i], ". With precursor sequence: ", as.character(precursor.set[[i]]), sep = ""))

    # check RNAStringSet position is same as df position
    if ( ! (df.mirdeep$precursor_seq[i] == as.character(precursor.set[[i]]) &
            df.mirdeep$mature_seq[i] == as.character(mature.set[[i]]) &
            df.mirdeep$star_seq[i] == as.character(star.set[[i]])) ){
        print(paste("Error at sequence:", i))
        print(df.mirdeep[i,])
        print(precursor.set[[i]])
        print(mature.set[[i]])
        print(star.set[[i]])
        stop("Construction of miRDeep2 GRanges object failed.
             Missmatch between data frame index and RNAStringSet index.")
    }

    # Flags to locate mature and star at start or end of precursor sequence
    matStartFlag <- FALSE
    matEndFlag <- FALSE
    starStartFlag <- FALSE
    starEndFlag <- FALSE

    #print("Finding mature sequence in precursor sequence.")

    # View of pattern matched mature into precursor
    matView <- matchPattern(mature.set[[i]], precursor.set[[i]])
    matView

    # check for no view
    if ( length(matView) == 0 ) {
        print(paste("Error at mature sequence:", i))
        print(df.mirdeep[i,])
        print(mature.set[[i]])
        stop("Construction of miRDeep2 GRanges object failed. Mature sequence not found in precursor sequence.")
    } else if ( length(matView) > 1 ) { # check for more than one view
        print("More than 1 view found")
        # check for begining of precursor
        if ( 1 %in% start(matView) ) {
            print("Mature seq found at start of precursor seq")
            matStartFlag <- TRUE
        }
        # check for end of precursor
        if ( width(precursor.set[i]) %in% end(matView) ) {
            print("Mature seq found at end of precursor seq")
            matEndFlag <- TRUE
        }
        # if both at beginning and end, error
        if ( matStartFlag & matEndFlag ) {
            print("Mature seq found at start and end of precursor seq!")
            print(paste("Error at mature sequence:", i))
            print(df.mirdeep[i,])
            print(mature.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
        # if neither beginning or end, error
        if ( ! ( matStartFlag | matEndFlag ) ) {
            print("Mature seq not found at start or end of precursor seq!")
            print(paste("Error at mature sequence:", i))
            print(df.mirdeep[i,])
            print(mature.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    } else { # only one view, find at start or end of precursor
        # check for begining of precursor
        if ( 1 == start(matView) ) {
            print("Mature seq found at start of precursor seq")
            matStartFlag <- TRUE
        }
        # check for end of precursor
        if ( width(precursor.set[i]) == end(matView) ) {
            print("Mature seq found at end of precursor seq")
            matEndFlag <- TRUE
        }
        # if both at beginning and end, error
        if ( matStartFlag & matEndFlag ) {
            print("Mature seq found at start and end of precursor seq!")
            print(paste("Error at mature sequence:", i))
            print(df.mirdeep[i,])
            print(mature.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
        # if neither beginning or end, error
        if ( ! ( matStartFlag | matEndFlag ) ) {
            print("Mature seq not found at start or end of precursor seq!")
            print(paste("Error at mature sequence:", i))
            print(df.mirdeep[i,])
            print(mature.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    }

    #print("Finding star sequence in precursor sequence.")

    # View of pattern matched star into precursor
    starView <- matchPattern(star.set[[i]], precursor.set[[i]])
    starView

    # check for no view
    if ( length(starView) == 0 ) {
        print(paste("Error at star sequence:", i))
        print(df.mirdeep[i,])
        print(star.set[[i]])
        stop("Construction of miRDeep2 GRanges object failed. Star sequence not found in precursor sequence.")
    } else if ( length(starView) > 1 ) { # check for more than one view
        print("More than 1 view found")
        # check for begining of precursor
        if ( 1 %in% start(starView) ) {
            print("Star seq found at start of precursor seq")
            starStartFlag <- TRUE
        }
        # check for end of precursor
        if ( width(precursor.set[i]) %in% end(starView) ) {
            print("Star seq found at end of precursor seq")
            starEndFlag <- TRUE
        }
        # if both at beginning and end, error
        if ( starStartFlag & starEndFlag ) {
            print("Star seq found at start and end of precursor seq!")
            print(paste("Error at star sequence:", i))
            print(df.mirdeep[i,])
            print(star.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
        # if neither beginning or end, error
        if ( ! ( starStartFlag | starEndFlag ) ) {
            print("Star seq not found at start or end of precursor seq!")
            print(paste("Error at star sequence:", i))
            print(df.mirdeep[i,])
            print(star.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    } else { # only one view, find at start or end of precursor
        # check for begining of precursor
        if ( 1 == start(starView) ) {
            print("Star seq found at start of precursor seq")
            starStartFlag <- TRUE
        }
        # check for end of precursor
        if ( width(precursor.set[i]) == end(starView) ) {
            print("Star seq found at end of precursor seq")
            starEndFlag <- TRUE
        }
        # if both at beginning and end, error
        if ( starStartFlag & starEndFlag ) {
            print("Star seq found at start and end of precursor seq!")
            print(paste("Error at star sequence:", i))
            print(df.mirdeep[i,])
            print(star.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
        # if neither beginning or end, error
        if ( ! ( starStartFlag | starEndFlag ) ) {
            print("Star seq not found at start or end of precursor seq!")
            print(paste("Error at star sequence:", i))
            print(df.mirdeep[i,])
            print(star.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    }

    # positive strand
    if (df.mirdeep$strand[i] == "+") {
        # mature sequence at beginning star at end
        if (matStartFlag & starEndFlag) {
            print("Mature Beginning, Star End")

            df.mirdeep$mature_start[i] <- df.mirdeep$precursor_start[i]
            df.mirdeep$mature_end[i] <- df.mirdeep$precursor_start[i] + width(mature.set[i]) - 1
            df.mirdeep$mature_width[i] <- width(mature.set[i])
            df.mirdeep$star_start[i] <- df.mirdeep$precursor_end[i] - width(star.set[i]) + 1
            df.mirdeep$star_end[i] <- df.mirdeep$precursor_end[i]
            df.mirdeep$star_width[i] <- width(star.set[i])

        } else if (starStartFlag & matEndFlag) { # star sequence at beginning mature at end
            print("Star Beginning, Mature End")

            df.mirdeep$mature_start[i] <- df.mirdeep$precursor_end[i] - width(mature.set[i]) + 1
            df.mirdeep$mature_end[i] <- df.mirdeep$precursor_end[i]
            df.mirdeep$mature_width[i] <- width(mature.set[i])
            df.mirdeep$star_start[i] <- df.mirdeep$precursor_start[i]
            df.mirdeep$star_end[i] <- df.mirdeep$precursor_start[i] + width(star.set[i]) - 1
            df.mirdeep$star_width[i] <- width(star.set[i])

        } else {
            print("Mature and star sequence not found on opposite ends of precursor sequence.")
            print(paste("Error at precursor sequence:", i))
            print(df.mirdeep[i,])
            print(precursor.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    }
    # negative strand
    if (df.mirdeep$strand[i] == "-") {
        # mature sequence at beginning star at end
        if (matStartFlag & starEndFlag) {
            print("Mature Beginning, Star End")

            df.mirdeep$mature_start[i] <- df.mirdeep$precursor_end[i] - width(mature.set[i]) + 1
            df.mirdeep$mature_end[i] <- df.mirdeep$precursor_end[i]
            df.mirdeep$mature_width[i] <- width(mature.set[i])
            df.mirdeep$star_start[i] <- df.mirdeep$precursor_start[i]
            df.mirdeep$star_end[i] <- df.mirdeep$precursor_start[i] + width(star.set[i]) - 1
            df.mirdeep$star_width[i] <- width(star.set[i])

        } else if (starStartFlag & matEndFlag) { # star sequence at beginning mature at end
            print("Star Beginning, Mature End")

            df.mirdeep$mature_start[i] <- df.mirdeep$precursor_start[i]
            df.mirdeep$mature_end[i] <- df.mirdeep$precursor_start[i] + width(mature.set[i]) - 1
            df.mirdeep$mature_width[i] <- width(mature.set[i])
            df.mirdeep$star_start[i] <- df.mirdeep$precursor_end[i] - width(star.set[i]) + 1
            df.mirdeep$star_end[i] <- df.mirdeep$precursor_end[i]
            df.mirdeep$star_width[i] <- width(star.set[i])

        } else {
            print("Mature and star sequence not found on opposite ends of precursor sequence.")
            print(paste("Error at precursor sequence:", i))
            print(df.mirdeep[i,])
            print(precursor.set[[i]])
            stop("Construction of miRDeep2 GRanges object failed.")
        }
    }

    # confirm coordinates and sequences match to genome sequence
    precursorSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                     names = df.mirdeep$chr[i],
                                                     start = df.mirdeep$precursor_start[i],
                                                     end = df.mirdeep$precursor_end[i],
                                                     strand = df.mirdeep$strand[i])))
    if ( precursorSeq != df.mirdeep$precursor_seq[i] ) {
        stop("Precursor coordinates do not match precursor sequence.")
    }
    matSeq <- as.character((RNAStringSet(getSeq(x = Hsapiens,
                                                names = df.mirdeep$chr[i],
                                                start = df.mirdeep$mature_start[i],
                                                end = df.mirdeep$mature_end[i],
                                                strand = df.mirdeep$strand[i]))))
    if ( matSeq != df.mirdeep$mature_seq[i] ) {
        stop("Mature coordinates do not match mature sequence.")
    }
    starSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                names = df.mirdeep$chr[i],
                                                start = df.mirdeep$star_start[i],
                                                end = df.mirdeep$star_end[i],
                                                strand = df.mirdeep$strand[i])))
    if ( starSeq != df.mirdeep$star_seq[i] ) {
        stop("Star coordinates do not match star sequence.")
    }

}

# separate into 3 dfs, mature, star, and precursor
df.mirdeep.precursor <- dplyr::select(df.mirdeep,
                                      seqnames = chr,
                                      start = precursor_start,
                                      end = precursor_end,
                                      strand,
                                      score,
                                      ID = name,
                                      sequence = precursor_seq)

df.mirdeep.precursor$Alias <- df.mirdeep.precursor$ID
df.mirdeep.precursor$source <- "miRDeep2"
df.mirdeep.precursor$type <- "miRNA_putative_precursor"
df.mirdeep.precursor$phase <- NA
df.mirdeep.precursor$Derives_from <- NA
df.mirdeep.precursor$Name <- df.mirdeep.precursor$ID

df.mirdeep.mature <- dplyr::select(df.mirdeep,
                                   seqnames = chr,
                                   start = mature_start,
                                   end = mature_end,
                                   strand,
                                   score,
                                   Derives_from = name,
                                   sequence = mature_seq)

df.mirdeep.mature$source <- "miRDeep2"
df.mirdeep.mature$type <- "miRNA_putative_mature"
df.mirdeep.mature$phase <- NA
df.mirdeep.mature$Alias <- paste(df.mirdeep.mature$Derives_from, "mat", sep="_")
df.mirdeep.mature$ID <- df.mirdeep.mature$Alias
df.mirdeep.mature$Name <- df.mirdeep.mature$Alias

df.mirdeep.star <- dplyr::select(df.mirdeep,
                                 seqnames = chr,
                                 start = star_start,
                                 end = star_end,
                                 strand,
                                 score,
                                 Derives_from = name,
                                 sequence = star_seq)

df.mirdeep.star$source <- "miRDeep2"
df.mirdeep.star$type <- "miRNA_putative_star"
df.mirdeep.star$phase <- NA
df.mirdeep.star$Alias <- paste(df.mirdeep.star$Derives_from, "star", sep="_")
df.mirdeep.star$ID <- df.mirdeep.star$Alias
df.mirdeep.star$Name <- df.mirdeep.star$Alias

print("Building miRDeep2 GRanges object.")
# combine into single data frame for creation of GRanges object
df.mirdeep.combination <- bind_rows(df.mirdeep.precursor, df.mirdeep.mature, df.mirdeep.star)
# format and reorder
df.mirdeep.combination$sequence <- RNAStringSet(df.mirdeep.combination$sequence)
df.mirdeep.combination$source <- factor(df.mirdeep.combination$source)
df.mirdeep.combination$type <- factor(df.mirdeep.combination$type)
df.mirdeep.combination <- dplyr::select(df.mirdeep.combination,
                                        seqnames,
                                        start,
                                        end,
                                        strand,
                                        source,
                                        type,
                                        score,
                                        phase,
                                        ID,
                                        Alias,
                                        Name,
                                        Derives_from,
                                        sequence)

gr.mirdeep <- makeGRangesFromDataFrame(df = df.mirdeep.combination,
                                       keep.extra.columns = TRUE,
                                       ignore.strand = FALSE)

# set seqinfo of gr
seqinfo(gr.mirdeep, new2old = match(seqlevels(Seqinfo(genome = "hg38")), seqlevels(gr.mirdeep))) <- Seqinfo(genome = "hg38")

# Metadata columns
# source (factor): miRDeep2
# type (factor): miRNA_putative_precursor or miRNA_putative_mature or miRNA_putative_star
# score (numeric): miRDeep2 score
# phase (int): NA
# ID (character): unique ID
# Alias (character): same as ID
# Name (character): same as ID
# Derives_from (character): NA for primary_transcript or ID of primary_transcript
# sequence (RNAStringSet)

print("miRDeep2 GRanges object complete.")

rm(df.mirdeep, df.mirdeep.combination, df.mirdeep.mature, df.mirdeep.precursor, df.mirdeep.star,
   mature.set, precursor.set, star.set, matView, starView, i)

gr.mirdeep

# Import miRge Novel miRNAs ###########################################################################################
print("Importing miRge putatively novel miRNAs.")
df.mirge <- read_tsv(mirge.novel.tsv)

# filter unique by position
df.mirge %<>% dplyr::filter(Unique_by_position)

df.mirge$Name <- paste("miRge", df.mirge$provisional_id, sep="_")
df.mirge$ID <- paste(df.mirge$Name, df.mirge$`Arm type`, sep="_")

df.mirge %<>%
    dplyr::select(seqnames = Chr,
                  mature_start = `Start Pos`,
                  mature_end = `End Pos`,
                  strand = Strand,
                  Name,
                  score = Probability,
                  ID,
                  precursor_seq = `Precursor miRNA sequence`)

df.mirge$mature_seq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                        names = df.mirge$seqnames,
                                                        start = df.mirge$mature_start,
                                                        end = df.mirge$mature_end,
                                                        strand = df.mirge$strand)))

df.mirge

print("Converting to GRanges object.")

# rna string set
precursor.set <- RNAStringSet(df.mirge$precursor_seq)
mature.set <- RNAStringSet(df.mirge$mature_seq)

# define precursor width
df.mirge$precursor_width <- width(precursor.set)
df.mirge$mature_width <- width(mature.set)

# variables to find with string matching
df.mirge$precursor_start <- NA
df.mirge$precursor_end <- NA

# loop over all precursors and define coordinates
# check coordinates correspond to correct genomic sequence
for (i in 1:length(df.mirge$Name)) {
    # Flags to locate mature at start or end of precursor sequence
    matStartFlag <- FALSE
    matEndFlag <- FALSE

    # View of pattern matched mature into precursor
    matView <- matchPattern(mature.set[[i]], precursor.set[[i]])
    #print(matView)
    if (length(matView) != 1) {
        stop("Incorrect number of views.")
    }
    distanceToStart <- start(matView) - 1

    distanceToEnd <- width(precursor.set[i]) - end(matView)

    if ( ! (width(precursor.set[i]) == distanceToStart + distanceToEnd + width(mature.set[i]) ) ) {
        stop("Distances don't match.")
    }

    # set start and end coords.
    if (df.mirge$strand[i] == "+") {
        df.mirge$precursor_start[i] <- df.mirge$mature_start[i] - distanceToStart
        df.mirge$precursor_end[i] <- df.mirge$mature_end[i] + distanceToEnd
    }
    if (df.mirge$strand[i] == "-") {
        df.mirge$precursor_start[i] <- df.mirge$mature_start[i] - distanceToEnd
        df.mirge$precursor_end[i] <- df.mirge$mature_end[i] + distanceToStart
    }

    # confirm coordinates and sequences match to genome sequence
    precursorSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                     names = df.mirge$seqnames[i],
                                                     start = df.mirge$precursor_start[i],
                                                     end = df.mirge$precursor_end[i],
                                                     strand = df.mirge$strand[i])))
    if ( precursorSeq != df.mirge$precursor_seq[i] ) {
        print(paste("index:", i))
        print(matView)
        stop("Precursor coordinates do not match precursor sequence.")
    }
    # confirm coordinates and sequences match to genome sequence
    matureSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                  names = df.mirge$seqnames[i],
                                                  start = df.mirge$mature_start[i],
                                                  end = df.mirge$mature_end[i],
                                                  strand = df.mirge$strand[i])))
    if ( matureSeq != df.mirge$mature_seq[i] ) {
        print(paste("index:", i))
        print(matView)
        stop("Mature coordinates do not match mature sequence.")
    }
}

# separate into  dfs, mature and precursor
df.mirge.precursor <- dplyr::select(df.mirge,
                                    seqnames,
                                    start = precursor_start,
                                    end = precursor_end,
                                    strand,
                                    score,
                                    Name,
                                    sequence = precursor_seq)

df.mirge.precursor$Alias <- df.mirge.precursor$Name
df.mirge.precursor$ID <- df.mirge.precursor$Name
df.mirge.precursor$source <- "miRge"
df.mirge.precursor$type <- "miRNA_putative_precursor"
df.mirge.precursor$phase <- NA
df.mirge.precursor$Derives_from <- NA

df.mirge.mature <- dplyr::select(df.mirge,
                                 seqnames,
                                 start = mature_start,
                                 end = mature_end,
                                 strand,
                                 score,
                                 Derives_from = Name,
                                 ID,
                                 sequence = mature_seq)

df.mirge.mature$source <- "miRge"
df.mirge.mature$type <- "miRNA_putative_mature"
df.mirge.mature$phase <- NA
df.mirge.mature$Alias <- df.mirge.mature$ID
df.mirge.mature$Name <- df.mirge.mature$ID

# combine into single data frame for creation of GRanges object
df.mirge.combination <- bind_rows(df.mirge.precursor, df.mirge.mature)
# format and reorder
df.mirge.combination$sequence <- RNAStringSet(df.mirge.combination$sequence)
df.mirge.combination$type <- factor(df.mirge.combination$type)
df.mirge.combination$source <- factor(df.mirge.combination$source)
df.mirge.combination <- dplyr::select(df.mirge.combination,
                                      seqnames,
                                      start,
                                      end,
                                      strand,
                                      source,
                                      type,
                                      score,
                                      phase,
                                      ID,
                                      Alias,
                                      Name,
                                      Derives_from,
                                      sequence)

gr.mirge <- makeGRangesFromDataFrame(df = df.mirge.combination,
                                     keep.extra.columns = TRUE,
                                     ignore.strand = FALSE)

# set seqinfo of gr
seqinfo(gr.mirge, new2old = match(seqlevels(Seqinfo(genome = "hg38")), seqlevels(gr.mirge))) <- Seqinfo(genome = "hg38")

# Metadata columns
# source (factor): miRge
# type (factor): miRNA_putative_precursor or miRNA_putative_mature
# score (numeric): miRge probability
# phase (int): NA
# ID (character): unique ID
# Alias (character): same as ID
# Name (character): same as ID
# Derives_from (character): NA for primary_transcript or ID of primary_transcript
# sequence (RNAStringSet)

print("miRge GRanges object complete.")

rm(df.mirge, df.mirge.combination, df.mirge.mature, df.mirge.precursor, mature.set, precursor.set, matView, i)

gr.mirge

# Import Friedlander miRNAs ###########################################################################################
print("Importing Friedlander 2014 novel miRNAs.")
df.fried <- read_excel(friedlander.xls)

# remove novels with no experimental evidence
df.fried %<>%
    dplyr::filter(!is.na(`Confidence leve`))

# select and rename columns
df.fried %<>%
    dplyr::select(ID = Identifier,
                  seqnames = Chromosome,
                  strand = Strand,
                  precursor_start = `Begin pos (1-indexed)`,
                  precursor_end = `End position (pos included)`,
                  mature_seq = `Mature sequence`,
                  precursor_seq = `Hairpin sequence`,
                  score = `Confidence leve`)

# convert strand word to strand symbol
df.fried$strand <- ifelse(df.fried$strand == "plus", "+", "-")

# convert to uppercase
df.fried$mature_seq <- toupper(df.fried$mature_seq)
df.fried$precursor_seq <- toupper(df.fried$precursor_seq)

# convert hg19 to hg38 coordinates
print("Converting hg19 to hg38 coordinates.")
# chain for hg19 to hg38 conversion
path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch <- import.chain(path)
# GRanges of fried precursor coordinates
gr <- GRanges(seqnames = df.fried$seqnames,
              ranges = IRanges(start = df.fried$precursor_start,
                               end = df.fried$precursor_end),
              strand = df.fried$strand)
# GRangesList of GRanges conversion
lo <- liftOver(gr, ch)
# Loop over all converted coordinates and replace them in data frame
for (i in 1:length(lo)) {
    #print(i)
    if (length(lo[[i]]) == 1) {
        df.fried$seqnames[i] <- as.character(seqnames(lo[[i]]))
        df.fried$strand[i] <- as.character(strand(lo[[i]]))
        df.fried$precursor_start[i] <- start(lo[[i]])
        df.fried$precursor_end[i] <- end(lo[[i]])
    }
    if (length(lo[[i]]) > 1) { # if more than one conversion, choose the first
        df.fried$seqnames[i] <- as.character(seqnames(lo[[i]]))[1]
        df.fried$strand[i] <- as.character(strand(lo[[i]]))[1]
        df.fried$precursor_start[i] <- start(lo[[i]])[1]
        df.fried$precursor_end[i] <- end(lo[[i]])[1]
    }
}

print("Converting to GRanges object.")

# rna string set
precursor.set <- RNAStringSet(df.fried$precursor_seq)
mature.set <- RNAStringSet(df.fried$mature_seq)

# define precursor width
df.fried$precursor_width <- width(precursor.set)
df.fried$mature_width <- width(mature.set)

# variables to find with string matching
df.fried$mature_start <- NA
df.fried$mature_end <- NA

# loop over all precursors and define coordinates
# check coordinates correspond to correct genomic sequence
for (i in 1:length(df.fried$ID)) {

    # View of pattern matched mature into precursor
    matView <- matchPattern(mature.set[[i]], precursor.set[[i]])
    #print(matView)
    if (length(matView) != 1) {
        print("Incorrect number of views.")
        print(paste("index:", i))

        distanceToStart <- start(matView) - 1

        distanceToEnd <- width(precursor.set[i]) - end(matView)

        whichView <- 0

        for (k in 1:length(matView)) {
            if (distanceToStart[k] == 10 | distanceToEnd[k] == 10) {
                whichView <- k
            }
        }

        if (whichView == 0) {
            stop("No good view")
        } else {
            matView <- matView[whichView]
        }
    }
    distanceToStart <- start(matView) - 1

    distanceToEnd <- width(precursor.set[i]) - end(matView)

    if ( ! (width(precursor.set[i]) == distanceToStart + distanceToEnd + width(mature.set[i]) ) ) {
        stop("Distances don't match.")
    }

    # set start and end coords.
    if (df.fried$strand[i] == "+") {
        df.fried$mature_start[i] <- df.fried$precursor_start[i] + distanceToStart
        df.fried$mature_end[i] <- df.fried$precursor_end[i] - distanceToEnd
    }
    if (df.fried$strand[i] == "-") {
        df.fried$mature_start[i] <- df.fried$precursor_start[i] + distanceToEnd
        df.fried$mature_end[i] <- df.fried$precursor_end[i] - distanceToStart
    }

    # confirm coordinates and sequences match to genome sequence
    precursorSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                     names = df.fried$seqnames[i],
                                                     start = df.fried$precursor_start[i],
                                                     end = df.fried$precursor_end[i],
                                                     strand = df.fried$strand[i])))
    if ( precursorSeq != df.fried$precursor_seq[i] ) {
        print(paste("index:", i))
        print(df.fried$precursor_seq[i])
        print(precursorSeq)
        print(matView)
        print("Precursor coordinates do not match precursor sequence.")
    }
    # confirm coordinates and sequences match to genome sequence
    matureSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                  names = df.fried$seqnames[i],
                                                  start = df.fried$mature_start[i],
                                                  end = df.fried$mature_end[i],
                                                  strand = df.fried$strand[i])))
    if ( matureSeq != df.fried$mature_seq[i] ) {
        print(paste("index:", i))
        print(df.fried$mature_seq[i])
        print(matureSeq)
        print(matView)
        print("Mature coordinates do not match mature sequence.")
    }
}

# separate into  dfs, mature and precursor
df.fried.precursor <- dplyr::select(df.fried,
                                    seqnames,
                                    start = precursor_start,
                                    end = precursor_end,
                                    strand,
                                    score,
                                    Name = ID,
                                    sequence = precursor_seq)

df.fried.precursor$Alias <- df.fried.precursor$Name
df.fried.precursor$ID <- df.fried.precursor$Name
df.fried.precursor$source <- "Friedlander2014"
df.fried.precursor$type <- "miRNA_putative_precursor"
df.fried.precursor$phase <- NA
df.fried.precursor$Derives_from <- NA

df.fried.mature <- dplyr::select(df.fried,
                                 seqnames,
                                 start = mature_start,
                                 end = mature_end,
                                 strand,
                                 score,
                                 Derives_from = ID,
                                 sequence = mature_seq)

df.fried.mature$source <- "Friedlander2014"
df.fried.mature$type <- "miRNA_putative_mature"
df.fried.mature$phase <- NA
df.fried.mature$Alias <- paste(df.fried.mature$Derives_from, "mat", sep="_")
df.fried.mature$ID <- df.fried.mature$Alias
df.fried.mature$Name <- df.fried.mature$Alias

# combine into single data frame for creation of GRanges object
df.fried.combination <- bind_rows(df.fried.precursor, df.fried.mature)
# format and reorder
df.fried.combination$sequence <- RNAStringSet(df.fried.combination$sequence)
df.fried.combination$source <- factor(df.fried.combination$source)
df.fried.combination$type <- factor(df.fried.combination$type)
df.fried.combination <- dplyr::select(df.fried.combination,
                                      seqnames,
                                      start,
                                      end,
                                      strand,
                                      source,
                                      type,
                                      score,
                                      phase,
                                      ID,
                                      Alias,
                                      Name,
                                      Derives_from,
                                      sequence)

gr.fried <- makeGRangesFromDataFrame(df = df.fried.combination,
                                     keep.extra.columns = TRUE,
                                     ignore.strand = FALSE)

# set seqinfo of gr
seqinfo(gr.fried, new2old = match(seqlevels(Seqinfo(genome = "hg38")), seqlevels(gr.fried))) <- Seqinfo(genome = "hg38")

print("Friedlander GRanges object complete.")

# Metadata columns
# source (factor): Friedlander2014
# type (factor): miRNA_putative_precursor or miRNA_putative_mature
# score (numeric): Friedlander evidence score
# phase (int): NA
# ID (character): unique ID
# Alias (character): same as ID
# Name (character): same as ID
# Derives_from (character): NA for primary_transcript or ID of primary_transcript
# sequence (RNAStringSet)

rm(df.fried, df.fried.combination, df.fried.mature, df.fried.precursor, mature.set, precursor.set, matView, i,
   k, lo, gr, ch)

gr.fried

# Import Nowakowski miRNAs ############################################################################################
print("Importing Nowakowski 2018 novel miRNAs.")
df.nowa <- read_xlsx(nowakowski.xlsx, sheet = "Table 12", skip = 4)

# select and rename columns
df.nowa %<>%
    dplyr::select(ID = `#miRNA`,
                  mature_pos = `chr_mature (hg19)`,
                  precursor_pos = `chr_hairpin (hg19)`,
                  mature_seq = seq_mature,
                  precursor_seq = seq_hairpin)

df.nowa$seqnames <- sapply(strsplit(df.nowa$precursor_pos, "_"), `[`, 1)
df.nowa$strand <- sapply(strsplit(df.nowa$precursor_pos, "_"), `[`, 2)
df.nowa$precursor_start <- as.integer(sapply(strsplit(df.nowa$precursor_pos, "_"), `[`, 3))
df.nowa$precursor_end <- as.integer(sapply(strsplit(df.nowa$precursor_pos, "_"), `[`, 4))
df.nowa$mature_start <- as.integer(sapply(strsplit(df.nowa$mature_pos, "_"), `[`, 3))
df.nowa$mature_end <- as.integer(sapply(strsplit(df.nowa$mature_pos, "_"), `[`, 4))

# convert to RNA
df.nowa$mature_seq <- as.character(RNAStringSet(DNAStringSet(df.nowa$mature_seq)))
df.nowa$precursor_seq <- as.character(RNAStringSet(DNAStringSet(df.nowa$precursor_seq)))

# convert hg19 to hg38 coordinates
print("Converting hg19 to hg38 coordinates.")
# chain for hg19 to hg38 conversion
path <- system.file(package="liftOver", "extdata", "hg19ToHg38.over.chain")
ch <- import.chain(path)
# GRanges of precursor coordinates
gr.precursor <- GRanges(seqnames = df.nowa$seqnames,
                        ranges = IRanges(start = df.nowa$precursor_start,
                                         end = df.nowa$precursor_end),
                        strand = df.nowa$strand)
# GRangesList of GRanges conversion
lo.precursor <- liftOver(gr.precursor, ch)
# Loop over all converted coordinates and replace them in data frame
for (i in 1:length(lo.precursor)) {
    #print(i)
    if (length(lo.precursor[[i]]) == 1) {
        df.nowa$seqnames[i] <- as.character(seqnames(lo.precursor[[i]]))
        df.nowa$strand[i] <- as.character(strand(lo.precursor[[i]]))
        df.nowa$precursor_start[i] <- start(lo.precursor[[i]])
        df.nowa$precursor_end[i] <- end(lo.precursor[[i]])
    }
    if (length(lo.precursor[[i]]) > 1) { # if more than one conversion, choose the first
        df.nowa$seqnames[i] <- as.character(seqnames(lo.precursor[[i]]))[1]
        df.nowa$strand[i] <- as.character(strand(lo.precursor[[i]]))[1]
        df.nowa$precursor_start[i] <- start(lo.precursor[[i]])[1]
        df.nowa$precursor_end[i] <- end(lo.precursor[[i]])[1]
    }
}

# GRanges of mature coordinates
gr.mature <- GRanges(seqnames = df.nowa$seqnames,
                     ranges = IRanges(start = df.nowa$mature_start,
                                      end = df.nowa$mature_end),
                     strand = df.nowa$strand)
# GRangesList of GRanges conversion
lo.mature <- liftOver(gr.mature, ch)
# Loop over all converted coordinates and replace them in data frame
for (i in 1:length(lo.mature)) {
    #print(i)
    if (length(lo.mature[[i]]) == 1) {
        df.nowa$seqnames[i] <- as.character(seqnames(lo.mature[[i]]))
        df.nowa$strand[i] <- as.character(strand(lo.mature[[i]]))
        df.nowa$mature_start[i] <- start(lo.mature[[i]])
        df.nowa$mature_end[i] <- end(lo.mature[[i]])
    }
    if (length(lo.precursor[[i]]) > 1) { # if more than one conversion, choose the first
        df.nowa$seqnames[i] <- as.character(seqnames(lo.mature[[i]]))[1]
        df.nowa$strand[i] <- as.character(strand(lo.mature[[i]]))[1]
        df.nowa$mature_start[i] <- start(lo.mature[[i]])[1]
        df.nowa$mature_end[i] <- end(lo.mature[[i]])[1]
    }
}

# check coordinates and sequences match
for (i in 1:length(df.nowa$ID)) {
    # confirm coordinates and sequences match to genome sequence
    precursorSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                     names = df.nowa$seqnames[i],
                                                     start = df.nowa$precursor_start[i],
                                                     end = df.nowa$precursor_end[i],
                                                     strand = df.nowa$strand[i])))
    if ( precursorSeq != df.nowa$precursor_seq[i] ) {
        print(paste("index:", i))
        print(df.nowa$precursor_seq[i])
        print(precursorSeq)
        print("Precursor coordinates do not match precursor sequence.")
    }
    # confirm coordinates and sequences match to genome sequence
    matureSeq <- as.character(RNAStringSet(getSeq(x = Hsapiens,
                                                  names = df.nowa$seqnames[i],
                                                  start = df.nowa$mature_start[i],
                                                  end = df.nowa$mature_end[i],
                                                  strand = df.nowa$strand[i])))
    if ( matureSeq != df.nowa$mature_seq[i] ) {
        print(paste("index:", i))
        print(df.nowa$mature_seq[i])
        print(matureSeq)
        print("Mature coordinates do not match mature sequence.")
    }
}

# separate into  dfs, mature and precursor
df.nowa.precursor <- dplyr::select(df.nowa,
                                   seqnames,
                                   start = precursor_start,
                                   end = precursor_end,
                                   strand,
                                   ID,
                                   sequence = precursor_seq)

df.nowa.precursor$ID <- sapply(strsplit(df.nowa.precursor$ID, "_"), `[`, 1)
df.nowa.precursor$Alias <- df.nowa.precursor$ID
df.nowa.precursor$Name <- df.nowa.precursor$ID
df.nowa.precursor$source <- "Nowakowski2018"
df.nowa.precursor$type <- "miRNA_putative_precursor"
df.nowa.precursor$phase <- NA
df.nowa.precursor$Derives_from <- NA
df.nowa.precursor$score <- NA

# one precursor is duplicated, remove it
df.nowa.precursor <- distinct(df.nowa.precursor)

df.nowa.mature <- dplyr::select(df.nowa,
                                seqnames,
                                start = mature_start,
                                end = mature_end,
                                strand,
                                ID,
                                sequence = mature_seq)

df.nowa.mature$source <- "Nowakowski2018"
df.nowa.mature$type <- "miRNA_putative_mature"
df.nowa.mature$phase <- NA
df.nowa.mature$score <- NA
df.nowa.mature$Derives_from <- sapply(strsplit(df.nowa.mature$ID, "_"), `[`, 1)
df.nowa.mature$Alias <- df.nowa.mature$ID
df.nowa.mature$Name <- df.nowa.mature$ID

# combine into single data frame for creation of GRanges object
df.nowa.combination <- bind_rows(df.nowa.precursor, df.nowa.mature)
# format and reorder
df.nowa.combination$sequence <- RNAStringSet(df.nowa.combination$sequence)
df.nowa.combination$source <- factor(df.nowa.combination$source)
df.nowa.combination$type <- factor(df.nowa.combination$type)
df.nowa.combination <- dplyr::select(df.nowa.combination,
                                     seqnames,
                                     start,
                                     end,
                                     strand,
                                     source,
                                     type,
                                     score,
                                     phase,
                                     ID,
                                     Alias,
                                     Name,
                                     Derives_from,
                                     sequence)

print("Converting to GRanges object.")

gr.nowa <- makeGRangesFromDataFrame(df = df.nowa.combination,
                                    keep.extra.columns = TRUE,
                                    ignore.strand = FALSE)

# set seqinfo of gr
seqinfo(gr.nowa, new2old = match(seqlevels(Seqinfo(genome = "hg38")), seqlevels(gr.nowa))) <- Seqinfo(genome = "hg38")

rm(ch, df.nowa, df.nowa.combination, df.nowa.mature, df.nowa.precursor, gr.mature, gr.precursor, lo.mature,
   lo.precursor, i, path, matureSeq, precursorSeq)

print("Nowakowski GRanges object complete.")

gr.nowa

# Create and Export GTF/GFF FIles #####################################################################################
# # mirdeep and mirge novel mirnas
# # combine novel granges objects into one
# gr.novel <- c(gr.mirge, gr.mirdeep)
#
# # remove score terms
# gr.novel$score <- NA
# # move type to seq_type, change type to miRNA
# gr.novel$seq_type <- gr.novel$type
# gr.novel$type <- "miRNA"
# # sequence as charaacter
# gr.novel$sequence <- as.character(gr.novel$sequence)
#
# # export as GTF/GFF
# rtracklayer::export(gr.novel, mirdeep.mirge.gtf.output.file)
#
# # mirdeep and mirge novel and friedlander mirnas
# # combine novel granges objects into one
# gr.novel <- c(gr.mirge, gr.mirdeep, gr.fried, gr.nowa)
#
# # remove score terms
# gr.novel$score <- NA
# # move type to seq_type, change type to miRNA
# gr.novel$seq_type <- gr.novel$type
# gr.novel$type <- "miRNA"
# # sequence as charaacter
# gr.novel$sequence <- as.character(gr.novel$sequence)
#
# # export as GTF/GFF
# rtracklayer::export(gr.novel, mirdeep.mirge.friedlander.nowakowski.gtf.output.file)
#
# # all known and novel mirnas
# gr.all <- c(gr.mirbase, gr.mirdeep, gr.mirge, gr.fried, gr.nowa)
#
# # remove score terms
# gr.all$score <- NA
# # move type to seq_type, change type to miRNA
# gr.all$seq_type <- gr.all$type
# gr.all$type <- "miRNA"
# # sequence as charaacter
# gr.all$sequence <- as.character(gr.all$sequence)
#
# # export as GTF/GFF
# rtracklayer::export(gr.all, all.known.novel.gtf.output.file)
#
# rm(gr.novel, gr.all)

# Export GRanges Objects ##############################################################################################

# saveRDS(gr.mirbase, mirbase.granges.rds.output.file)
# saveRDS(gr.mirdeep, mirdeep.granges.rds.output.file)
# saveRDS(gr.mirge, mirge.granges.rds.output.file)
# saveRDS(gr.fried, friedlander.granges.rds.output.file)
# saveRDS(gr.nowa, nowakowski.granges.rds.output.file)

# Non-Overlapping Set #################################################################################################
# mirbase_v22 > friedlander > nowakowski > mirge > mirdeep2
# ignore strand with novel annotation because we didn't use stranded RNA-seq

#findOverlaps(gr.mirge, gr.mirbase)

# find overlaps only between mature sequences
gr.mirbase.mirna <- gr.mirbase[gr.mirbase$type == "miRNA"]
gr.mirge.mirna <- gr.mirge[gr.mirge$type != "miRNA_putative_precursor"]
gr.mirdeep.mirna <- gr.mirdeep[gr.mirdeep$type != "miRNA_putative_precursor"]
gr.fried.mirna <- gr.fried[gr.fried$type != "miRNA_putative_precursor"]
gr.nowa.mirna <- gr.nowa[gr.nowa$type != "miRNA_putative_precursor"]

# mirbase overlaps to self
findOverlaps(gr.mirbase.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = FALSE)

# friedlander overlaps to self
overlaps.fried.self <- findOverlaps(gr.fried.mirna,
                                    drop.self = TRUE,
                                    drop.redundant = TRUE,
                                    ignore.strand = TRUE)
# remove overlaps to self
gr.fried.mirna <- gr.fried.mirna[! 1:length(gr.fried.mirna) %in% queryHits(overlaps.fried.self)]


# nowakowski overlaps to self
findOverlaps(gr.nowa.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = TRUE)

# mirge overlaps to self
overlaps.mirge.self <- findOverlaps(gr.mirge.mirna,
                                    drop.self = TRUE,
                                    drop.redundant = TRUE,
                                    ignore.strand = TRUE)
# remove overlaps to self
gr.mirge.mirna <- gr.mirge.mirna[! 1:length(gr.mirge.mirna) %in% queryHits(overlaps.mirge.self)]

# mirdeep overlaps to self
overlaps.mirdeep.self <- findOverlaps(gr.mirdeep.mirna,
                                      drop.self = TRUE,
                                      drop.redundant = TRUE,
                                      ignore.strand = TRUE)
# remove overlaps to self
gr.mirdeep.mirna <- gr.mirdeep.mirna[! 1:length(gr.mirdeep.mirna) %in% queryHits(overlaps.mirdeep.self)]


# merges
# merge friedlander into mirbase
overlaps <- findOverlaps(gr.fried.mirna, gr.mirbase.mirna,
                         ignore.strand = TRUE)
gr.combo.mirna <- c(gr.mirbase.mirna, gr.fried.mirna[! 1:length(gr.fried.mirna) %in% queryHits(overlaps)])

findOverlaps(gr.combo.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = FALSE)

# merge nowakowski into combo
overlaps <- findOverlaps(gr.nowa.mirna, gr.combo.mirna,
                         ignore.strand = TRUE)
gr.combo.mirna <- c(gr.combo.mirna, gr.nowa.mirna[! 1:length(gr.nowa.mirna) %in% queryHits(overlaps)])

findOverlaps(gr.combo.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = FALSE)

# merge mirge into combo
overlaps <- findOverlaps(gr.mirge.mirna, gr.combo.mirna,
                         ignore.strand = TRUE)
gr.combo.mirna <- c(gr.combo.mirna, gr.mirge.mirna[! 1:length(gr.mirge.mirna) %in% queryHits(overlaps)])

findOverlaps(gr.combo.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = FALSE)

# merge mirdeep into combo
overlaps <- findOverlaps(gr.mirdeep.mirna, gr.combo.mirna)
gr.combo.mirna <- c(gr.combo.mirna, gr.mirdeep.mirna[! 1:length(gr.mirdeep.mirna) %in% queryHits(overlaps)])

findOverlaps(gr.combo.mirna,
             drop.self = TRUE,
             drop.redundant = TRUE,
             ignore.strand = FALSE)



# all known and novel mirnas, including precursor/primary mirnas
gr.all <- c(gr.mirbase, gr.mirdeep, gr.mirge, gr.fried, gr.nowa)

# all precursor/primary sequences within the combo mature list
derives.from <- unique(gr.combo.mirna$Derives_from)

# select from full list of mirnas just the precursor/primary
gr.precursor <- gr.all[gr.all$ID %in% derives.from]

# add precursors to combo list
gr.combo <- c(gr.combo.mirna, gr.precursor)

sum(duplicated(gr.combo$ID))

# create unique name of pri-mir + mir
# find derives from name
gr.combo$Derives_from_name <- gr.combo$Name[match(gr.combo$Derives_from, gr.combo$ID)]
# create UniqueName
gr.combo$UniqueName <- ifelse(is.na(gr.combo$Derives_from_name), gr.combo$Name, paste(gr.combo$Derives_from_name, gr.combo$Name, sep = "_"))
# some unique names are still duplicated because the mirbase database has errors
sum(duplicated(gr.combo$UniqueName))
# modify uniqueNames of only the duplicated ones
dup.index <- which(duplicated(gr.combo$UniqueName) | duplicated(gr.combo$UniqueName, fromLast = TRUE))
gr.combo[dup.index]$UniqueName <- make.unique(gr.combo$UniqueName[dup.index])
# check for duplicates
sum(duplicated(gr.combo$uniqueName))


# # output gr
# gr.output <- gr.combo
# # remove score terms
# gr.output$score <- NA
# # move type to seq_type, change type to miRNA
# gr.output$seq_type <- gr.output$type
# gr.output$type <- "miRNA"
# # sequence as charaacter
# gr.output$sequence <- as.character(gr.output$sequence)
#
# # export as GTF/GFF
# rtracklayer::export(gr.output, all.known.novel.non.overlapping.gtf.output.file)
# export as granges
saveRDS(gr.combo, all.known.novel.non.overlapping.rds.output.file)



