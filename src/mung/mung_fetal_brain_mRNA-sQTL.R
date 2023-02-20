# munge fetal brain sQTL data from Nil

# Raw association data provided by Nil from:
# /proj/steinlab/projects/R00/eQTLanalysis/Raw_data/splicing/bulk/

# Significant conditionally independent sQTLs (per intron junction):
# /proj/steinlab/projects/R00/eQTLanalysis/sumstats/sQTL/bulk_sqtl_sumstat.csv

# Genotype directory to get major/minor alleles
# /proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/
# chromosome files labeled as chr*_F.uniq.bim
# columns: chromosome_number SNP 0 position minor_allele major_allele
# major_allele is the effect allele

# Per sample splice percent
# Quantile normalized
# /proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/bulk_unique_perind.counts.gz.qqnorm_chr$i (fetal bulk)
# Raw ratios
# /proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/bulk_unique_perind.counts.gz.phen_chr$i (fetal bulk)

library(here)
library(readr)
library(dplyr)
library(magrittr)

# OUTPUT FILES #########################################################################################################
# directory for raw compiled results by chr
raw.output.dir <- here("results/external_data/fetal_brain_sQTL/raw/")
dir.create(raw.output.dir, recursive = TRUE, showWarnings = FALSE)

# compiled sQTLs rds file
compiled.sqtls.output.rds <- here("results/external_data/fetal_brain_sQTL/fetal_brain_sQTLs.rds")

# directory for LD for each sQTL index snp
ld.output.dir <- here("results/external_data/fetal_brain_sQTL/ld/")
dir.create(ld.output.dir, recursive = TRUE, showWarnings = FALSE)

# index snp list
index.snps.txt <- here("results/external_data/fetal_brain_sQTL/index.snps.txt")

# matrix of splice ratio across donors
splice.percent.raw.output.rds <- here("results/external_data/fetal_brain_sQTL/fetal_brain_sQTL_raw_ratios.rds")
splice.percent.norm.output.rds <- here("results/external_data/fetal_brain_sQTL/fetal_brain_sQTL_norm_ratios.rds")

# INPUT FILES ##########################################################################################################
# raw association data by chr
raw.data.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/Raw_data/splicing/bulk/"

# Significant SNP-junction pairs (conditional sQTLs):
sqtls.csv <- "/proj/steinlab/projects/R00/eQTLanalysis/sumstats/sQTL/bulk_sqtl_sumstat.csv"

# Genotype directory
genotypes.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/"

# Raw and normalized ratios dir
ratios.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/Splicing_QTL/"

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import sQTLs #########################################################################################################

df.sqtls <- read_csv(sqtls.csv)

# df.cond %<>%
#     select(SNP = snp,
#            ENSG = gene,
#            BETA = beta,
#            P = pvalue,
#            CHR = chr,
#            RANK = rank,
#            BETA.CONDITIONAL = cond.beta,
#            P.CONDITIONAL = cond.pval)

# Raw Data #############################################################################################################


# # new dataframe to build conditional eqtl table with alleles and positions
# df.new <- tibble()
#
# for (chr in CHROMS) {
#
#     print(chr)
#
#     df.genos <- NULL
#     df.raw <- NULL
#     df.tmp <- NULL
#     raw.output.file <- NULL
#
#     # read in variant information for this chromosome
#     print("Reading in variants...")
#     if (chr == "chrX") {
#         df.par <- read_table(paste0(genotypes.dir, chr, ".PAR_F.uniq.bim"),
#                               col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
#         df.nonpar <- read_table(paste0(genotypes.dir, chr, ".Non.PAR_F.uniq.bim"),
#                                  col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
#         df.genos <- bind_rows(df.par, df.nonpar)
#         rm(df.par, df.nonpar)
#     } else {
#         df.genos <- read_table(paste0(genotypes.dir, chr, "_F.uniq.bim"),
#                                 col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
#     }
#     print("Finished.")
#
#     # read in raw association data for this chromosome
#     print("Reading in raw association data...")
#     df.raw <- read_table(paste0(raw.data.dir, "F-", chr, ".txt.gz"), col_names = c("SNP", "BETA", "P", "ENSG"))
#     print("Finished.")
#
#     print("Modifying columns...")
#     # isolate ENSG string
#     df.raw %<>%
#         mutate(ENSG = sapply(strsplit(sapply(strsplit(ENSG, "-"), `[`, 3), "\\."), `[`, 1))
#
#     # Join with genotype data to get effect allele and base pair positions
#     df.raw %<>%
#         left_join(df.genos, by = "SNP")
#
#     # modify chrom column
#     df.raw %<>%
#         mutate(CHR = chr)
#
#     # select and order columns
#     df.raw %<>%
#         select(CHR,
#                BP,
#                SNP,
#                BETA,
#                P,
#                ENSG,
#                ALLELE_MINOR,
#                ALLELE_MAJOR_EFFECT)
#
#     # check for NAs in dataframe
#     if (any(is.na(df.raw))) {
#         stop("NAs present in df.raw after merge!")
#     }
#
#     print("Finished.")
#
#     # filter conditional eqtls for this chrom
#     print("Getting conditional eqtls...")
#     df.cond %>%
#         filter(CHR == chr) -> df.tmp
#
#     # join with raw information to get position and allele information
#     df.tmp %<>%
#         left_join(df.raw, by = c("SNP", "ENSG"))
#
#     if (any(is.na(df.tmp))) {
#         stop("NAs present in df.tmp after merge!")
#     }
#     if (! all(df.tmp$BETA.x == df.tmp$BETA.y)) {
#         stop("BETAs not the same in df.tmp!")
#     }
#     if (! all(df.tmp$P.x == df.tmp$P.y)) {
#         stop("Ps not the same in df.tmp!")
#     }
#     if (! all(df.tmp$CHR.x == df.tmp$CHR.y)) {
#         stop("CHRs not the same in df.tmp!")
#     }
#
#     # select needed columns
#     df.tmp %<>%
#         select(SNP,
#                ENSG,
#                BETA = BETA.x,
#                P = P.x,
#                CHR = CHR.x,
#                BP,
#                RANK,
#                BETA.CONDITIONAL,
#                P.CONDITIONAL,
#                ALLELE_MINOR,
#                ALLELE_MAJOR_EFFECT)
#
#     # grow table of conditional eqtls
#     df.new %<>%
#         bind_rows(df.tmp)
#
#     print("Finished.")
#
#     # save raw data as RDS
#     raw.output.file <- paste0(raw.output.dir, chr, "_fetal_brain_mRNA-eQTL_raw.rds")
#
#     print(paste0("Saving raw data to: ", raw.output.file))
#     saveRDS(df.raw, raw.output.file)
#     print("Finished")
#
# }

# save sqtl table
print(paste0("Saving conditional sqtls to: ", compiled.sqtls.output.rds))
saveRDS(df.sqtls, paste0(compiled.sqtls.output.rds))
print("Finished.")

# Get LD for Index SNPs ################################################################################################
df.sqtls <- readRDS(compiled.sqtls.output.rds)

# export index snp list
write_lines(df.sqtls$snp, index.snps.txt)

# run the shell script file to get ld using plink

# Ratios Files #################################################################################################

# for each "intron" loop over, pull raw and normalized data, build table

df.raw.ratios <- tibble()
df.norm.ratios <- tibble()


for (chr in CHROMS) {

    print(paste("Working on", chr))

    # introns on this chrom
    df.sqtls %>%
        filter(chr == !!chr) %>%
        pull(intron) -> introns

    introns <- unique(introns)

    if (chr == "chrX") {
        chr <- "chr23"
    }



    print("Importing...")

    df.this.chrom.raw <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.phen_", chr))

    df.this.chrom.norm <- read_tsv(paste0(ratios.dir, "bulk_unique_perind.counts.gz.qqnorm_", chr))

    print("Checking all IDs present...")
    all(introns %in% df.this.chrom.raw$ID)
    all(introns %in% df.this.chrom.norm$ID)


    df.this.chrom.raw %<>%
        filter(ID %in% introns) %>%
        select(ID, ends_with(".unique"))

    df.this.chrom.norm %<>%
        filter(ID %in% introns) %>%
        select(ID, ends_with(".unique"))

    df.raw.ratios %<>%
        bind_rows(df.this.chrom.raw)

    df.norm.ratios %<>%
        bind_rows(df.this.chrom.norm)




}


all(df.raw.ratios$ID %in% df.sqtls$intron)
all(df.norm.ratios$ID %in% df.sqtls$intron)

colnames(df.raw.ratios) <- sapply(strsplit(colnames(df.raw.ratios), "\\."), `[`, 1)
colnames(df.norm.ratios) <- sapply(strsplit(colnames(df.norm.ratios), "\\."), `[`, 1)

write_rds(df.raw.ratios, splice.percent.raw.output.rds)
write_rds(df.norm.ratios, splice.percent.norm.output.rds)


