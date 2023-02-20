# munge blood miRNA-eQTL data from Huan 2015

# only significant snps reported on hg19

library(here)
library(readr)
library(dplyr)
library(readxl)
library(magrittr)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)

# OUTPUT FILES #########################################################################################################
# directory for results
output.dir <- here("results/external_data/huan_2015_blood_miRNA-eQTL/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# compiled eQTLs rds file
compiled.eqtls.output.rds <- here("results/external_data/huan_2015_blood_miRNA-eQTL/huan_2015_blood_miRNA-eQTLs.rds")

# row eQTLs rds file
raw.eqtls.output.rds <- here("results/external_data/huan_2015_blood_miRNA-eQTL/huan_2015_blood_miRNA-eQTL_hg38_variants_dataFrame.rds")

# directory for LD for each mRNA-eQTL index snp
ld.output.dir <- here("results/external_data/huan_2015_blood_miRNA-eQTL/ld/")
dir.create(ld.output.dir, recursive = TRUE, showWarnings = FALSE)

# index snp list
index.snps.txt <- here("results/external_data/huan_2015_blood_miRNA-eQTL/index.snps.txt")

# INPUT FILES ##########################################################################################################
# blood miRNA-eQTLs
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# hg19 to hg38 chain file
chain.file <- here("data/hg19ToHg38.over.chain")



# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import eQTLs #########################################################################################################

df.blood <- read_xlsx(huan.eqtls.xlsx, skip = 1)

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

df.blood.hg38 <- as_tibble(as.data.frame(gr.huan.hg38))

rm(df.blood, gr.huan.hg38)

# export raw data
write_rds(df.blood.hg38, raw.eqtls.output.rds)

# Find Index SNPs ######################################################################################################

# list of emirs
emirs <- unique(df.blood.hg38$hsa_miR_name)

df.eqtls <- tibble()

# loop over emirs and get eqtl
for (i in 1:length(emirs)) {

    df.blood.hg38 %>%
        filter(hsa_miR_name == emirs[i]) -> df.this.emir

    df.this.emir %>%
        top_n(wt = Pval, n = -1) -> df.tmp

    df.tmp <- df.tmp[1,]

    df.eqtls %<>%
        bind_rows(df.tmp)

}

# export eqtl list
write_rds(df.eqtls, compiled.eqtls.output.rds)

snp.list <- paste(df.eqtls$seqnames, df.eqtls$snpID, sep = ":")

# export index snp list
write_lines(unique(snp.list), index.snps.txt)











# system(command = "/bin/bash -c \"source ~/.bashrc && module list && module add plink && module list\"")





