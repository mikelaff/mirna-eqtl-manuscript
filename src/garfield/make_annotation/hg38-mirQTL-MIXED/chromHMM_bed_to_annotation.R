# process chromHMM files for use as garfield annotations

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(rtracklayer)

# OUTPUT FILES #########################################################################################################
output.dir <- here("data/garfield-data/annotation/hg38-mirQTL-MIXED/chromHMM_fetal_brain/")

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# INPUT FILES ##########################################################################################################
# directory with maf and tss distance variants. only these variants will be used to create annotation files
maftssd.dir <- here("data/garfield-data/maftssd/hg38-mirQTL-MIXED/")

# chromHMM bed files
fetal.brain.male.bed <- here("data/chromHMM/E081_15_coreMarks_hg38lift_mnemonics.bed.gz")
fetal.brain.female.bed <- here("data/chromHMM/E082_15_coreMarks_hg38lift_mnemonics.bed.gz")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Bed Files #####################################################################################################

gr.male <- rtracklayer::import(fetal.brain.male.bed)
gr.female <- rtracklayer::import(fetal.brain.female.bed)

# modify annotation names
gr.male$name <- paste(gr.male$name, "male", sep = "_")
gr.female$name <- paste(gr.female$name, "female", sep = "_")

# Format Annotations ###################################################################################################

# unique annotation names from granges
df.anno <- tibble(Annotation = c(unique(gr.male$name), unique(gr.female$name)))

# deconstruct names
df.anno$number <- as.integer(sapply(strsplit(df.anno$Annotation, "_"), `[`, 1))
df.anno$sex <- sapply(strsplit(df.anno$Annotation, "_"), `[`, 3)

# arrange into meaningful order
df.anno %<>%
    arrange(sex, number)

# add additional columns
df.anno %<>%
    mutate(Index = seq(0, nrow(df.anno)-1),
           Celltype = "PrimaryTissue",
           Tissue = "FetalBrain",
           Type = "ChromHMM",
           Category = "Core_15-state_Mnemonics")

df.anno$Tissue <- paste(df.anno$Tissue, df.anno$sex, sep = "_")

# select columns for link_file.txt
df.anno %<>%
    dplyr::select(Index,
                  Annotation,
                  Celltype,
                  Tissue,
                  Type,
                  Category)

# Process Variants #####################################################################################################

# loop over chroms
for (chrom in CHROMS) {
    print(paste("Processing variants from chromosome:", chrom))

    # load mafftssd file
    print("Loading variants...")
    df.vars <- read_table2(paste0(maftssd.dir, chrom), col_names = c("BP", "MAF", "TSSD"), col_types = "idi")
    df.vars$seqnames <- chrom

    df.vars %<>%
        dplyr::select(seqnames,
                      BP)

    # convert to GRanges
    gr.vars <- makeGRangesFromDataFrame(df.vars,
                                        ignore.strand = TRUE,
                                        seqnames.field = "seqnames",
                                        start.field = "BP",
                                        end.field = "BP")

    print("Finding overlaps...")
    # overlap with male annotations
    hits.male <- GenomicRanges::findOverlaps(gr.vars, gr.male, select = "first", ignore.strand = TRUE)
    df.vars$male_annotation <- gr.male$name[hits.male]

    # overlap with female annotations
    hits.female <- GenomicRanges::findOverlaps(gr.vars, gr.female, select = "first", ignore.strand = TRUE)
    df.vars$female_annotation <- gr.female$name[hits.female]

    # get annotation index
    df.vars$male_index <- df.anno$Index[match(df.vars$male_annotation, df.anno$Annotation)]
    df.vars$female_index <- df.anno$Index[match(df.vars$female_annotation, df.anno$Annotation)]

    # make -1 index for NA values
    df.vars$male_index[is.na(df.vars$male_index)] <- -1
    df.vars$female_index[is.na(df.vars$female_index)] <- -1

    print("Formatting annotation strings...")
    # create annotation string
    df.vars$anno_string <- paste(rep("0", max(df.anno$Index) + 1), collapse = "")
    # replace 0 with 1 at anno index, -1's will replace at position 0, meaning no replacement
    substr(df.vars$anno_string, df.vars$male_index + 1, df.vars$male_index + 1) <- "1"
    substr(df.vars$anno_string, df.vars$female_index + 1, df.vars$female_index + 1) <- "1"

    print("Writing annotation file...")
    # export position and annotation string
    df.vars %<>%
        dplyr::select(BP, anno_string)
    write_delim(df.vars, path = paste0(output.dir, chrom), col_names = FALSE, delim = " ")

}

# Export Link File #####################################################################################################

print("Writing link_file.txt")
write_delim(df.anno, path = paste0(output.dir, "link_file.txt"), delim = " ")
