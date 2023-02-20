# find co-localization miRNA-eQTL to blood miRNA-eQTLs

# for each miRNA-eQTL:
#       look for overlap with miRNA-eQTL at r2 >= 0.8

library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# working directory for co-localization analysis
dir.working <- paste0(here("results/co-localization/"), "blood_miRNA-eQTL", "/")
dir.create(dir.working, recursive = TRUE, showWarnings = FALSE)

# table of overlaps
overlap.output.rds <- paste0(dir.working, "blood_miRNA-eQTL", "_mirQTL_overlaps_r2at0.8.rds")

# list of miR-eQTL index snps that have overlaps
index.snps.txt <- paste0(dir.working, "blood_miRNA-eQTL", "_overlap.index.snps.txt")

# table of emirs and their sets
emir.sets.rds <- paste0(dir.working, "blood_miRNA-eQTL", "mirQTL_blood-miRNA-eQTL_emiR_sets.rds")

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
mirQTL.df.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")
# blood miRNA-eQTLs
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# Directory for LD at each mirQTL index SNP
mirQTL.ld.dir <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/ld/")

# hg19 to hg38 chain file
chain.file <- here("data/hg19ToHg38.over.chain")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import eQTL Data #####################################################################################################

df.mirQTL <- readRDS(mirQTL.df.rds)

df.mirQTL %>%
    filter(SIGNIFICANCE == "eigenMT_fdr5percent") -> df.mirQTL.high

df.blood <- read_xlsx(huan.eqtls.xlsx, skip = 1)

blood.emirs <- unique(df.blood$hsa_miR_name)

brain.emirs <- unique(df.mirQTL.high$NAME)

joint.emirs <- brain.emirs[brain.emirs %in% blood.emirs]

brain.unique.emirs <- brain.emirs[! brain.emirs %in% blood.emirs]

blood.unique.emirs <- blood.emirs[! blood.emirs %in% brain.emirs]

all.emirs <- unique(c(blood.emirs, brain.emirs))

df.emirs <- tibble(emir = all.emirs)

df.emirs$set <- NA
df.emirs$set[df.emirs$emir %in% brain.unique.emirs] <- "brain_only"
df.emirs$set[df.emirs$emir %in% blood.unique.emirs] <- "blood_only"
df.emirs$set[df.emirs$emir %in% joint.emirs] <- "brain_and_blood"

df.emirs$coloc <- FALSE


# chain for hg19 to hg38 conversion
ch <- import.chain(chain.file)

# GRanges of huan eqtls
gr <- makeGRangesFromDataFrame(df.blood,
                               seqnames.field = "chr.SNP",
                               start.field = "SNP.pos",
                               end.field = "SNP.pos",
                               strand.field = "SNP.strand",
                               ignore.strand = TRUE,
                               keep.extra.columns = TRUE)
# GRangesList of GRanges conversion
lo <- liftOver(gr, ch)

gr.huan.hg38 <- unlist(lo)

rm(lo, gr, ch)

df.blood <- as_tibble(as.data.frame(gr.huan.hg38))







df.mirQTL.high %>%
    filter(NAME %in% df.blood$hsa_miR_name) -> df.mirQTL.possible



# Find Overlaps ########################################################################################################

# build table of overlaps
df.overlaps <- tibble()

# loop over each miRNA-eQTL
for (i in 1:nrow(df.mirQTL.possible)) {

    printMessage(df.mirQTL.possible$eQTL[i])

    df.blood.this.chrom <- NULL
    df.blood.possible <- NULL

    mirQTL.ld.file <- NULL
    df.mirQTL.ldbuddies <- NULL

    # filter for blood QTLs on this chromosome
    df.blood %>%
        filter(seqnames == df.mirQTL.possible$SNP.CHR[i]) -> df.blood.this.chrom

    # filter for blood QTLs within 1MB of this mirQTL
    bp.start <- df.mirQTL.possible$SNP.BP.hg38[i] - 1e6
    bp.end <- df.mirQTL.possible$SNP.BP.hg38[i] + 1e6
    df.blood.this.chrom %>%
        filter(start >= bp.start & start <= bp.end) -> df.blood.possible

    # if no possible overlaps, go to next mirQTL
    if (nrow(df.blood.possible) == 0) {
        print(paste("No overlaps possible for:", df.mirQTL$eQTL[i]))
        next()
    }

    # get LD buddies for mirQTL index SNP
    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    mirQTL.ld.file <- paste0(mirQTL.ld.dir, "conditional.mirQTLor.index.", df.mirQTL.possible$eSNP[i], ".ld")

    # ld for this index SNP
    suppressWarnings(
        df.mirQTL.ldbuddies <- read_table(mirQTL.ld.file, col_types = cols())
    )
    df.mirQTL.ldbuddies %<>%
        filter(R2 >= 0.8)

    if (any(df.mirQTL.ldbuddies$BP_B %in% df.blood.possible$start)) {

        print("Overlap found!")

        tmp.overlaps <- NULL
        tmp.mirQTL <- NULL
        tmp.mQTL <- NULL

        # all overlaps
        overlap.blood.indexes <- which(df.blood.possible$start %in% df.mirQTL.ldbuddies$BP_B)
        overlap.snps <- df.mirQTL.ldbuddies$SNP_B[df.mirQTL.ldbuddies$BP_B %in% df.blood.possible$start]
        print(overlap.snps)

        df.emirs$coloc[match(df.mirQTL.possible$NAME[i], df.emirs$emir)] <- TRUE


        tmp.overlaps <- tibble(overlap.snps = paste(overlap.snps, collapse = ","))

        # build overlap dataframe
        tmp.mirQTL <- df.mirQTL.possible[i, ]
        colnames(tmp.mirQTL) <- paste0(colnames(tmp.mirQTL), ".mirQTL")
        tmp.blood <- df.blood.possible[overlap.blood.indexes, ]
        colnames(tmp.blood) <- paste0(colnames(tmp.blood), ".bloodQTL")
        df.overlaps %<>%
            bind_rows(bind_cols(tmp.mirQTL, tmp.blood, tmp.overlaps))


    }

}

# Export Overlaps ######################################################################################################

saveRDS(df.overlaps, overlap.output.rds)

write_lines(unique(df.overlaps$eSNP.mirQTL), index.snps.txt)



# plot ######

df.emirs %<>%
    mutate(set = factor(set, levels = c("blood_only", "brain_only", "brain_and_blood")),
           coloc = factor(coloc, levels = c(TRUE, FALSE), ordered = TRUE))



df.emirs %>%
    ggplot(aes(x = set, fill = coloc)) +
    geom_bar() +
    scale_fill_manual(values = c("darkorange", "navy")) +
    labs(y = "emiR Count",
         fill = "Co-localized eQTL") +
    plotTheme("figure")

ggsave(here("doc/paper/figure4/pdf/emir_sets.pdf"), height = 2, width = 2)


write_rds(df.emirs, emir.sets.rds)

# Target Prediction, GO Analysis ############################











