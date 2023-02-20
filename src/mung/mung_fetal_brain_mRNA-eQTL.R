# munge fetal brain eQTL data from Nil

# Raw association data provided by Nil from:
# /proj/steinlab/projects/R00/eQTLanalysis/Raw_data/bulk/
# chromosome files labeled as F-chr*.txt.gz
# columns: SNP BETA P ENSG_string

# Significant eSNP-eGene pairs (primary eQTLs):
# /proj/steinlab/projects/R00/eQTLanalysis/multiple_test_correction/eigenMT/input/bulk/sig.genes
# after checking the files, all the primary eQTLs in this file are also in the conditional eQTLs file

# Significant eSNP-eGene pairs (conditional eQTLs):
# /proj/steinlab/projects/R00/eQTLanalysis/conditional_analysis/eqtl_bulk_real.parameters.csv

# Genotype directory to get major/minor alleles
# /proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/
# chromosome files labeled as chr*_F.uniq.bim
# columns: chromosome_number SNP 0 position minor_allele major_allele
# major_allele is the effect allele

library(here)
library(readr)
library(dplyr)
library(magrittr)

# OUTPUT FILES #########################################################################################################
# directory for raw compiled results by chr
raw.output.dir <- here("results/external_data/fetal_brain_mRNA-eQTL/raw/")
dir.create(raw.output.dir, recursive = TRUE, showWarnings = FALSE)

# compiled eQTLs rds file
compiled.eqtls.output.rds <- here("results/external_data/fetal_brain_mRNA-eQTL/fetal_brain_mRNA-eQTLs.rds")

# directory for LD for each mRNA-eQTL index snp
ld.output.dir <- here("results/external_data/fetal_brain_mRNA-eQTL/ld/")
dir.create(ld.output.dir, recursive = TRUE, showWarnings = FALSE)

# index snp list
index.snps.txt <- here("results/external_data/fetal_brain_mRNA-eQTL/index.snps.txt")

# INPUT FILES ##########################################################################################################
# raw association data by chr
raw.data.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/Raw_data/bulk/"

# Significant eSNP-eGene pairs (primary eQTLs):
#eqtls.txt <- "/proj/steinlab/projects/R00/eQTLanalysis/multiple_test_correction/eigenMT/input/bulk/sig.genes"

# Significant eSNP-eGene pairs (conditional eQTLs):
conditional.eqtls.csv <- "/proj/steinlab/projects/R00/eQTLanalysis/conditional_analysis/eqtl_bulk_real.parameters.csv"

# Genotype directory
genotypes.dir <- "/proj/steinlab/projects/R00/eQTLanalysis/genofiles/allele_order_corrected_genofiles/bulk/"

# mRNA Expression Data (used to check effect alleles)
rse.mrna.rds <- here("results/rdata_files/20200803_rse_gene_counts.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import eQTLs #########################################################################################################

#df.eqtls <- read_table2(eqtls.txt)

df.cond <- read_csv(conditional.eqtls.csv)

df.cond %<>%
    select(SNP = snp,
           ENSG = gene,
           BETA = beta,
           P = pvalue,
           CHR = chr,
           RANK = rank,
           BETA.CONDITIONAL = cond.beta,
           P.CONDITIONAL = cond.pval)

# Raw Data #############################################################################################################


# new dataframe to build conditional eqtl table with alleles and positions
df.new <- tibble()

for (chr in CHROMS) {

    print(chr)

    df.genos <- NULL
    df.raw <- NULL
    df.tmp <- NULL
    raw.output.file <- NULL

    # read in variant information for this chromosome
    print("Reading in variants...")
    if (chr == "chrX") {
        df.par <- read_table(paste0(genotypes.dir, chr, ".PAR_F.uniq.bim"),
                              col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
        df.nonpar <- read_table(paste0(genotypes.dir, chr, ".Non.PAR_F.uniq.bim"),
                                 col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
        df.genos <- bind_rows(df.par, df.nonpar)
        rm(df.par, df.nonpar)
    } else {
        df.genos <- read_table(paste0(genotypes.dir, chr, "_F.uniq.bim"),
                                col_names = c("CHR", "SNP", "CM", "BP", "ALLELE_MINOR" , "ALLELE_MAJOR_EFFECT"))
    }
    print("Finished.")

    # read in raw association data for this chromosome
    print("Reading in raw association data...")
    df.raw <- read_table(paste0(raw.data.dir, "F-", chr, ".txt.gz"), col_names = c("SNP", "BETA", "P", "ENSG"))
    print("Finished.")

    print("Modifying columns...")
    # isolate ENSG string
    df.raw %<>%
        mutate(ENSG = sapply(strsplit(sapply(strsplit(ENSG, "-"), `[`, 3), "\\."), `[`, 1))

    # Join with genotype data to get effect allele and base pair positions
    df.raw %<>%
        left_join(df.genos, by = "SNP")

    # modify chrom column
    df.raw %<>%
        mutate(CHR = chr)

    # select and order columns
    df.raw %<>%
        select(CHR,
               BP,
               SNP,
               BETA,
               P,
               ENSG,
               ALLELE_MINOR,
               ALLELE_MAJOR_EFFECT)

    # check for NAs in dataframe
    if (any(is.na(df.raw))) {
        stop("NAs present in df.raw after merge!")
    }

    print("Finished.")

    # filter conditional eqtls for this chrom
    print("Getting conditional eqtls...")
    df.cond %>%
        filter(CHR == chr) -> df.tmp

    # join with raw information to get position and allele information
    df.tmp %<>%
        left_join(df.raw, by = c("SNP", "ENSG"))

    if (any(is.na(df.tmp))) {
        stop("NAs present in df.tmp after merge!")
    }
    if (! all(df.tmp$BETA.x == df.tmp$BETA.y)) {
        stop("BETAs not the same in df.tmp!")
    }
    if (! all(df.tmp$P.x == df.tmp$P.y)) {
        stop("Ps not the same in df.tmp!")
    }
    if (! all(df.tmp$CHR.x == df.tmp$CHR.y)) {
        stop("CHRs not the same in df.tmp!")
    }

    # select needed columns
    df.tmp %<>%
        select(SNP,
               ENSG,
               BETA = BETA.x,
               P = P.x,
               CHR = CHR.x,
               BP,
               RANK,
               BETA.CONDITIONAL,
               P.CONDITIONAL,
               ALLELE_MINOR,
               ALLELE_MAJOR_EFFECT)

    # grow table of conditional eqtls
    df.new %<>%
        bind_rows(df.tmp)

    print("Finished.")

    # save raw data as RDS
    raw.output.file <- paste0(raw.output.dir, chr, "_fetal_brain_mRNA-eQTL_raw.rds")

    print(paste0("Saving raw data to: ", raw.output.file))
    saveRDS(df.raw, raw.output.file)
    print("Finished")

}

# save conditional eqtl table
print(paste0("Saving conditional eqtls to: ", compiled.eqtls.output.rds))
saveRDS(df.new, paste0(compiled.eqtls.output.rds))
print("Finished.")

# Get LD for Index SNPs ################################################################################################
df.eqtls <- readRDS(compiled.eqtls.output.rds)

# export index snop list
write_lines(df.eqtls$SNP, index.snps.txt)

# run the shell script file to get ld using plink

# Check Effect Alleles #################################################################################################
stop()
# check effect alleles against the donor genotypes and expression data

library(DESeq2)

df.eqtls <- readRDS(compiled.eqtls.output.rds)

# import expression data
rse.mrna <- readRDS(rse.mrna.rds)
dds <- DESeqDataSet(rse.mrna, design = ~1)
vsd <- vst(dds)

df.effect.results <- tibble()

# loop over each eqtl, pull samples, genotypes, and expression data. calculate beta and compare
for (i in 1:nrow(df.eqtls)) {

    chr <- df.eqtls$CHR[i]








}



# system(command = "/bin/bash -c \"source ~/.bashrc && module list && module add plink && module list\"")





