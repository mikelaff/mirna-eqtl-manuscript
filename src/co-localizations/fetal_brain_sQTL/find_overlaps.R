# find co-localization miRNA-eQTL to fetal brain sQTLs

# for each miRNA-eQTL:
#       look for overlap with sQTL at r2 >= 0.8

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
# working directory for co-localization analysis
dir.working <- paste0(here("results/co-localization/"), "fetal_brain_sQTL", "/")
dir.create(dir.working, recursive = TRUE, showWarnings = FALSE)

# table of overlaps
overlap.output.rds <- paste0(dir.working, "fetal_brain_sQTL", "_mirQTL_overlaps_r2at0.8.rds")

# list of miR-eQTL index snps that have overlaps
index.snps.txt <- paste0(dir.working, "fetal_brain_sQTL", "_overlap.index.snps.txt")

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
mirQTL.df.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")
# mQTL sQTLs
mQTL.df.rds <- here("results/external_data/fetal_brain_sQTL/fetal_brain_sQTLs.rds")

# Directory for LD at each mirQTL index SNP
mirQTL.ld.dir <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/ld/")
# Directory for LD at each mQTL index SNP
mQTL.ld.dir <- here("results/external_data/fetal_brain_sQTL/ld/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import eQTL Data #####################################################################################################

df.mirQTL <- readRDS(mirQTL.df.rds)

df.mQTL <- readRDS(mQTL.df.rds)
df.mQTL %<>%
    mutate(sQTL = paste(snp, intron, sep = "_"))

# Find Overlaps ########################################################################################################

# build table of overlaps
df.overlaps <- tibble()

# loop over each sQTL
for (i in 1:nrow(df.mirQTL)) {

    printMessage(df.mirQTL$eQTL[i])

    df.mQTL.this.chrom <- NULL
    df.mQTL.possible <- NULL

    mirQTL.ld.file <- NULL
    df.mirQTL.ldbuddies <- NULL

    # filter for mQTLs on this chromosome
    df.mQTL %>%
        filter(chr == df.mirQTL$SNP.CHR[i]) -> df.mQTL.this.chrom

    # filter for mQTLs within 1MB of this mirQTL
    bp.start <- df.mirQTL$SNP.BP.hg38[i] - 1e6
    bp.end <- df.mirQTL$SNP.BP.hg38[i] + 1e6
    df.mQTL.this.chrom %>%
        filter(BP >= bp.start & BP <= bp.end) -> df.mQTL.possible

    # if no possible overlaps, go to next mirQTL
    if (nrow(df.mQTL.possible) == 0) {
        print(paste("No overlaps possible for:", df.mirQTL$eQTL[i]))
        next()
    }

    # get LD buddies for mirQTL index SNP
    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    mirQTL.ld.file <- paste0(mirQTL.ld.dir, "conditional.mirQTLor.index.", df.mirQTL$eSNP[i], ".ld")

    # ld for this index SNP
    suppressWarnings(
    df.mirQTL.ldbuddies <- read_table(mirQTL.ld.file, col_types = cols())
    )
    df.mirQTL.ldbuddies %<>%
        filter(R2 >= 0.8)

    # for each possible overlap
    for (j in 1:nrow(df.mQTL.possible)) {

        mQTL.ld.file <- NULL
        df.mQTL.ldbuddies <- NULL

        # get LD buddies for mQTL index SNP
        # ld file: mQTL.index.chrX:19203648:C:T.ld
        mQTL.ld.file <- paste0(mQTL.ld.dir, "sQTL.index.", df.mQTL.possible$snp[j], ".ld")

        # ld for this index SNP
        suppressWarnings(
            df.mQTL.ldbuddies <- read_table(mQTL.ld.file, col_types = cols())
        )
        df.mQTL.ldbuddies %<>%
            filter(R2 >= 0.8)

        if (any(df.mirQTL.ldbuddies$BP_B %in% df.mQTL.ldbuddies$BP_B)) {

            print("Overlap found!")

            tmp.overlaps <- NULL
            tmp.mirQTL <- NULL
            tmp.mQTL <- NULL

            # all overlaps
            overlap.snps <- df.mirQTL.ldbuddies$SNP_B[df.mirQTL.ldbuddies$BP_B %in% df.mQTL.ldbuddies$BP_B]
            tmp.overlaps <- tibble(overlap.snps = paste(overlap.snps, collapse = ","))

            # build overlap dataframe
            tmp.mirQTL <- df.mirQTL[i, ]
            colnames(tmp.mirQTL) <- paste0(colnames(tmp.mirQTL), ".mirQTL")
            tmp.mQTL <- df.mQTL.possible[j, ]
            colnames(tmp.mQTL) <- paste0(colnames(tmp.mQTL), ".mQTL")
            df.overlaps %<>%
                bind_rows(bind_cols(tmp.mirQTL, tmp.mQTL, tmp.overlaps))

        }

    }

}

# Export Overlaps ######################################################################################################

saveRDS(df.overlaps, overlap.output.rds)

write_lines(unique(df.overlaps$eSNP.mirQTL), index.snps.txt)





