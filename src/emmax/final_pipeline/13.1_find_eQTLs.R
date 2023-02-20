# get emiRs and eSNPs to get eQTLs.

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(DESeq2)
library(GenomicRanges)
library(mikelaffr)

date.prefix <- format(Sys.time(), "%Y%m%d")
date.prefix <- "20200120"

# association results directory name
results.name <- "20200120_mirQTLor"

# OUTPUT FILES ########################################################################################################
# eQTL output files
eqtl.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eQTLs_dataFrame.rds")
eqtl.gr.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eQTLs_GRanges.rds")

# emiR output files
emir.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_emiRs_dataFrame.rds")
emir.gr.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_emiRs_GRanges.rds")

# eSNP output files
esnp.df.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eSNPs_dataFrame.rds")
esnp.gr.rds <- paste0(here("results/emmax/association_results/"), results.name, "/compiled/", results.name, "_eSNPs_GRanges.rds")

# INPUT FILES ##########################################################################################################
# Clumped results directory
results.dir <- paste0(here("results/emmax/association_results/"), results.name, "/clumped/")

# Summarized association results for every variant within 1MB of each expressed miRNA
summarized.results.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# non-overlapping granges of known and novel mirnas
#granges.rds <- here("data/gtf_and_granges/20191203_all_known_and_novel_mirna_non_overlapping_granges.rds")

# vsd
vsd.rds <- here("results/emmax/phenotype_files/20200120_mirQTLor_VST_miRNA_expression_residual/20200120_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Seqinfo object
SEQINFO <- Seqinfo(genome = "hg38")[CHROMS]

# Import Clumped Data ##################################################################################################

clump.files <- list.files(results.dir, pattern = "*.clumped")

# loop over clump files and build summary data frame
df.eqtls <- tibble()

for (clump.file in clump.files) {
    print(clump.file)

    emir <- NULL
    df.clump <- NULL
    df.tmp <- NULL

    emir <- strsplit(clump.file, "\\.")[[1]][5]

    df.clump <- read_table(paste0(results.dir, clump.file),
                           col_names = c("CHR", "F", "SNP", "BP", "P", "TOTAL", "NSIG", "S05", "S01", "S001", "S0001", "SP2"),
                           skip = 1)

    df.tmp <- tibble(eSNP = df.clump$SNP, emiR = emir)

    df.eqtls <- bind_rows(df.eqtls, df.tmp)

}

rm(df.tmp, df.clump)

# add eqtl column
df.eqtls %<>%
    mutate(eQTL = paste(emiR, eSNP, sep = "+"))


print(paste("Number of eQTLs:", sum(!duplicated(df.eqtls$eQTL))))
print(paste("Number of eSNPs:", sum(!duplicated(df.eqtls$eSNP))))
print(paste("Number of emiRs:", sum(!duplicated(df.eqtls$emiR))))

# Import Variant Summarized Results ####################################################################################
df.variants <- as_tibble(readRDS(summarized.results.dataframe.rds))

df.variants %<>%
    mutate(UniName_SNP = paste(UniName, SNP, sep = "+"))

# Import miRNA GRanges #################################################################################################
vsd <- readRDS(vsd.rds)

gr.mirs <- rowRanges(vsd)

rm(vsd)

df.mirs <- as_tibble(gr.mirs)

df.mirs %<>%
    select(UniName = uniqueName,
           miR.CHR = seqnames,
           miR.START.hg38 = start,
           miR.END.hg38 = end,
           miR.WIDTH = width,
           miR.STRAND = strand,
           SOURCE = source,
           TYPE = type,
           ID,
           ALIAS = Alias,
           NAME = Name,
           DERIVES_FROM = Derives_from,
           DERIVES_FROM_NAME = Derives_from_name,
           SEQ = sequence)

#all(df.eqtls$emiR %in% df.mirs$UniName)

# Format eQTL Table ####################################################################################################

# join with SNP information
df.eqtls %<>%
    left_join(df.variants, by = c("eQTL" = "UniName_SNP")) %>%
    select(eQTL,
           emiR,
           eSNP,
           BETA,
           P,
           SNP.CHR = CHR,
           SNP.BP.hg38 = BP.hg38,
           EFFECT.ALLELE,
           REF,
           ALT,
           ALT_CTS,
           OBS_CT,
           A1,
           A2,
           A1.HOM.COUNT = A1.HOM.count,
           HET.COUNT = HET.count,
           A2.HOM.COUNT = A2.HOM.count)

# join with miR information
df.eqtls %<>%
    left_join(df.mirs, by = c("emiR" = "UniName"))

# create GRanges
gr.eqtls <- makeGRangesFromDataFrame(df.eqtls,
                                     keep.extra.columns = TRUE,
                                     ignore.strand = TRUE,
                                     seqinfo = SEQINFO,
                                     seqnames.field = "SNP.CHR",
                                     start.field = "SNP.BP.hg38",
                                     end.field = "SNP.BP.hg38")
names(gr.eqtls) <- gr.eqtls$eQTL

# Format eSNP Table ####################################################################################################

# select only SNP information
df.eqtls %>%
    select(eSNP,
           CHR = SNP.CHR,
           BP.hg38 = SNP.BP.hg38,
           EFFECT.ALLELE,
           REF,
           ALT,
           ALT_CTS,
           OBS_CT,
           A1,
           A2,
           A1.HOM.COUNT,
           HET.COUNT,
           A2.HOM.COUNT) -> df.esnps

# remove duplicates
df.esnps %<>%
    distinct()

# create GRanges
gr.esnps <- makeGRangesFromDataFrame(df.esnps,
                                     keep.extra.columns = TRUE,
                                     ignore.strand = TRUE,
                                     seqinfo = SEQINFO,
                                     seqnames.field = "CHR",
                                     start.field = "BP.hg38",
                                     end.field = "BP.hg38")
names(gr.esnps) <- gr.esnps$eSNP

# Format emiR Table ####################################################################################################

# select only miR information
df.eqtls %>%
    select(emiR,
           CHR = miR.CHR,
           START.hg38 = miR.START.hg38,
           END.hg38 = miR.END.hg38,
           WIDTH = miR.WIDTH,
           STRAND = miR.STRAND,
           SOURCE,
           TYPE,
           ID,
           ALIAS,
           NAME,
           DERIVES_FROM,
           DERIVES_FROM_NAME,
           SEQ) -> df.emirs

# remove duplicates
df.emirs %<>%
    distinct()

# create GRanges
gr.emirs <- makeGRangesFromDataFrame(df.emirs,
                                     keep.extra.columns = TRUE,
                                     ignore.strand = FALSE,
                                     seqinfo = SEQINFO,
                                     seqnames.field = "CHR",
                                     start.field = "START.hg38",
                                     end.field = "END.hg38",
                                     strand.field = "STRAND")
names(gr.emirs) <- gr.emirs$emiR

# Export Tables and GRanges ############################################################################################

print("Saving dataFrames and GRanges .rds files...")

saveRDS(df.eqtls, eqtl.df.rds)
saveRDS(gr.eqtls, eqtl.gr.rds)

saveRDS(df.esnps, esnp.df.rds)
saveRDS(gr.esnps, esnp.gr.rds)

saveRDS(df.emirs, emir.df.rds)
saveRDS(gr.emirs, emir.gr.rds)
