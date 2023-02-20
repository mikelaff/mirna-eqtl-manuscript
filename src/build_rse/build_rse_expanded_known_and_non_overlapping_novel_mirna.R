# build ranged summarized experiment with mirna count data for both miRBase and novel mirnas
# mirge quantified counts are expanded to multiple genomic loci

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(SummarizedExperiment)
library(Biostrings)
library(RMySQL)

# OUTPUT FILES #########################################################################################################
output.se.rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"),
                        "_rse_expanded_known_and_non_overlapping_novel_mirna_counts.rds", sep="")

# INPUT FILES ##########################################################################################################
# mirge quantified expression counts for miRBase miRNAs
mirbase.counts.tsv <- here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv")
# featureCounts quantified counts from bowtie mapped bams for novel miRNAs
novel.counts.tsv <- here("results/counts/small_rna_seq/20190501_mirdeep_mirge_friedlander_nowakowski_counts.tsv")

# sample metadata
samples.metadata.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")

# mirge mirna merges
mirge.merges.csv <- here("data/mirge2.0/miRge.Libs/human/annotation.Libs/human_merges_miRBase.csv")

# Non-overlapping miRNA annotations for known and novel miRNAs
known.novel.granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Import Non-overlapping miRNA #########################################################################################
gr.mirna <- readRDS(known.novel.granges.rds)
names(gr.mirna) <- gr.mirna$UniqueName

# Import Counts and Format #############################################################################################
# read in count data frame
df.counts <- read_tsv(mirbase.counts.tsv, col_types = cols(.default = "i", miRNA = "c"))

# some rows have SNP countifications, not sure why, they don't have any counts, so i'm removing them
df.counts %<>% filter(!grepl("SNP", miRNA))
# miRge merged some miRNAs because their sequences are the same, this includes different SNP versions of a miRNA
# I want to expand the merged miRNAs so that both genomic loci attributed to a similar sequence will both have counts
# I don't want to expand SNP version, those should stay collapsed.
mirge.merges <- read_lines(mirge.merges.csv)
# loop over merges and expand to both miRNA names
for (mir in mirge.merges) {
    # look for slash, otherwise merge is of SNPs and dont want to expand
    if (grepl("/", mir)) {
        # the mir which will be expanded
        expand.mir <- unlist(strsplit(mir, ","))[1]
        # a list of mirs to expand into
        merge.list <- unlist(strsplit(mir, ","))[-1]
        # check expand mir is in list only once
        if (sum(grepl(expand.mir, df.counts$miRNA)) == 1) {
            # row to be expanded
            copyRow <- filter(df.counts, miRNA == expand.mir)
            # new df with expanded rows
            df.new <- data.frame(miRNA = merge.list, copyRow[2:241], stringsAsFactors = FALSE)
            # bind rows to df.counts
            df.counts <- bind_rows(list(df.counts, df.new))
            # remove old mir that was expanded
            df.counts %<>% filter(miRNA != expand.mir)

        } else {
            stop("error")
        }
    } else {
        next
    }
}
rm(df.new, copyRow, mir, merge.list, mirge.merges, expand.mir)

# romove SNP suffix, get rid of repeats
df.counts$miRNA <- sapply(strsplit(df.counts$miRNA, "\\."), `[`, 1)
df.counts %<>% distinct()

# create count matrix
mat.counts <- as.matrix(df.counts[,2:241])
rownames(mat.counts) <- df.counts$miRNA

# Import novel counts
df.counts.novel <- read_tsv(novel.counts.tsv)
# remove precursor counts
df.counts.novel %<>%
    filter(seq_type != "miRNA_putative_precursor")

# keep only novel counts in the non-overlapping list
df.counts.novel %<>%
    filter(name %in% gr.mirna$Name)

# create count matrix
mat.counts.novel <- as.matrix(df.counts.novel[,11:250])
rownames(mat.counts.novel) <- df.counts.novel$name

stopifnot(all(colnames(mat.counts) == colnames(mat.counts.novel)))

# bind count matrices
mat.counts <- rbind(mat.counts, mat.counts.novel)
rm(df.counts, df.counts.novel, mat.counts.novel)

# SteinLab Database ####################################################################################################

# database connection
con <- dbConnect(MySQL(),
                 user = user_name,
                 password = user_password,
                 dbname = db_name,
                 host = host_name)

# get donors table for gestation week
donors <- dbGetQuery(con, "SELECT * FROM `Donors`")
rm(con, user_name, user_password, db_name, host_name)

# Load Metadata ########################################################################################################
df.metadata <- read_tsv(samples.metadata.tsv, col_types = cols(Lanes = "c"))

# combine with donor table to get gestation week and tissue acquisition date
df.metadata %<>%
    left_join(dplyr::select(donors, DonorID, TissueAcquisitionDate = AcquisitionDate, GestationWeek), by = "DonorID")

# rearrange columns
df.metadata %<>%
    dplyr::select(RNAID, DonorID, GestationWeek, TissueSection, TissueAcquisitionDate, everything())

rm(donors)

# format colData
# assign row names
col.data <- as.data.frame(df.metadata)
row.names(col.data) <- col.data$RNAID
rm(df.metadata)

# Import miRge Merged miRNAs ###########################################################################################
# Import mirge merged miRNAs to be added to row data
mirge.merges <- read_lines(mirge.merges.csv)

df.merges <- data.frame(miRNA = sapply(strsplit(mirge.merges, ","), `[`, 1),
                        merge.list = sapply(strsplit(mirge.merges, ","), function(x) paste(x[-1], collapse = ",")),
                        stringsAsFactors = FALSE)

stopifnot(all(df.merges$mirna %in% rownames(mat.counts)))

# Match Rows to GRanges ################################################################################################
# rows of count matrix not in granges object
# most likely 5p and 3p version of mirbase mirnas that are only annotated as
# miR w/o 5p and 3p destinctions in mirbase, and mirge merges
# build data frame to search for name of GRange
df.mirs <- data.frame(row_name = rownames(mat.counts),
                      stringsAsFactors = FALSE)

# name w/o 5p or 3p for secondary match
df.mirs$base <- paste("hsa-miR", sapply(strsplit(df.mirs$row_name, "-"), `[`, 3), sep="-")
# name of hairpin for tertiary match
df.mirs$hairpin <- paste("hsa-mir", sapply(strsplit(df.mirs$row_name, "-"), `[`, 3), sep="-")

# primary, secondary, or tertiary index into gr.mirna
df.mirs$gr_index <- match(df.mirs$row_name, gr.mirna$Name)
df.mirs$gr_index2 <- match(df.mirs$base, gr.mirna$Name)
df.mirs$gr_index3 <- match(df.mirs$hairpin, gr.mirna$Name)
# select index to use
df.mirs$gr_ind <- ifelse(!is.na(df.mirs$gr_index),
                         df.mirs$gr_index,
                         ifelse(!is.na(df.mirs$gr_index2),
                                df.mirs$gr_index2,
                                NA))

# mirs to remove from count matrix (only 4)
mat.counts <- mat.counts[df.mirs$row_name[!is.na(df.mirs$gr_ind)],]
# remove last 4 NAs
df.mirs <- df.mirs[!is.na(df.mirs$gr_ind),]
# final rowranges for construction of ranged summarized experiment
row.ranges <- gr.all[df.mirs$gr_ind]

rm(gr.all, df.mirs)





# format rowData
row.data <- data.frame(miRNA = rownames(mat.counts),
                       stringsAsFactors = FALSE)
row.data %<>%
    left_join(df.merges, by = "miRNA")
rownames(row.data) <- row.data$miRNA

# Build SummarizedExperiment ######################################################################

mode(mat.counts) <- "integer"

# format and check col.data
col.data <- col.data[colnames(mat.counts),]
stopifnot(all(rownames(col.data) == colnames(mat.counts)))


stopifnot(all(rownames(row.data) == rownames(mat.counts)))

# build se
se <- SummarizedExperiment(assays = SimpleList(counts = mat.counts),
                           rowData = row.data,
                           colData = col.data)

# save se
saveRDS(se, output.se.rds)


