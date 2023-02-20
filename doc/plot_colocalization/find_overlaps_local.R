# gviz colocalization scratch

library(here)
library(dplyr)
library(readxl)
library(readr)
library(magrittr)
#library(Gviz)
library(GenomicRanges)
#library(GenomicAlignments)
#library(DESeq2)
#library(biomaRt)
library(rtracklayer)
library(liftOver)

# OUTPUT ##############################################################################################################
#dir.pdf <- here("doc/plot_colocalization/pdfs/")

#dir.create(file.path(dir.pdf))

# INPUT FILES ##########################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# mirQTL eSNPs
esnps.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_dataFrame.rds")
esnps.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eSNPs_GRanges.rds")

# mirQTL emiRs
emirs.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_dataFrame.rds")
emirs.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_emiRs_GRanges.rds")

# Directory for LD at each mirQTL index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# gwas dataset root dir
dir.gwas <- "/Users/mikelaff/Projects/BACKUP_EXCLUDED/gwas_datasets/"

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

GWAS.PVAL.THRESHOLD <- 5e-8

HG19.SEQINFO <- Seqinfo(genome = "hg19")

# Import ###############################################################################################################
# mirQTL eqtls
df.mirQTL.eqtls <- readRDS(eqtls.dataframe.rds)

# loop over eqtl, get ld buddies, combine into table
df.mirQTL.eqtls.ldbuddies <- tibble()

for (i in 1:nrow(df.mirQTL.eqtls)) {

    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, df.mirQTL.eqtls$SNP.CHR[i], ".hardcall.prefiltered.mirQTLor.", df.mirQTL.eqtls$eSNP[i], ".ld")

    # ld for this index SNP
    df.mirQTL.ldbuddies <- read_table(ld.file, col_types = cols())
    df.mirQTL.ldbuddies %<>%
        mutate(CHR_A = paste0("chr", CHR_A),
               CHR_B = paste0("chr", CHR_B))
    df.mirQTL.ldbuddies %<>%
        dplyr::rename(BP.hg38_A = BP_A,
                      BP.hg38_B = BP_B)
    df.mirQTL.ldbuddies %<>%
        dplyr::mutate(eQTL = df.mirQTL.eqtls$eQTL[i])

    # combine with eqtl table
    df.mirQTL.ldbuddies %<>%
        left_join(df.mirQTL.eqtls, by = "eQTL")



    # r2 > 0.8
    df.mirQTL.ldbuddies %<>%
        dplyr::filter(R2 > 0.8)

    df.mirQTL.eqtls.ldbuddies <- bind_rows(df.mirQTL.eqtls.ldbuddies, df.mirQTL.ldbuddies)

    rm(df.mirQTL.ldbuddies)

}

rm(i, ld.file)

# change chr23 to chrX
df.mirQTL.eqtls.ldbuddies %<>%
    mutate(CHR_A = replace(CHR_A, CHR_A == "chr23", "chrX"),
           CHR_B = replace(CHR_B, CHR_B == "chr23", "chrX"))

# create granges from ld snps
gr.eqtl.snps.in.ld <- makeGRangesFromDataFrame(df.mirQTL.eqtls.ldbuddies,
                                               keep.extra.columns = TRUE,
                                               ignore.strand = TRUE,
                                               seqinfo = Seqinfo(genome = "hg38"),
                                               start.field = "BP.hg38_B",
                                               end.field = "BP.hg38_B",
                                               seqnames.field = "CHR_B")

# liftover to hg19
# chain for hg38 to hg19 conversion
path <- system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch <- import.chain(path)

# GRangesList of GRanges conversion
lo <- liftOver(gr.eqtl.snps.in.ld, ch)

gr.eqtl.snps.in.ld.hg19 <- unlist(lo)
seqlevels(gr.eqtl.snps.in.ld.hg19) <- CHROMS

seqinfo(gr.eqtl.snps.in.ld.hg19) <- HG19.SEQINFO
gr.eqtl.snps.in.ld.hg19 <- keepSeqlevels(gr.eqtl.snps.in.ld.hg19,
                                         CHROMS,
                                         pruning.mode = "coarse")

rm(ch, lo, path)

# Colocalizations ######################################################################################################

# load GWAS snps (hg19), look for overlaps with eqtl ld buddies

# AD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        df.gwas1 <- read_delim(paste0(dir.gwas, "ad/Kunkle_etal_Stage1_results.txt?file=1"), delim = " ")
        df.gwas2 <- read_delim(paste0(dir.gwas, "ad/Kunkle_etal_Stage2_results.txt?file=1"), delim = " ")
        df.gwas <- bind_rows(df.gwas1, df.gwas2)

        rm(df.gwas1, df.gwas2)

        df.gwas.sig <- filter(df.gwas, Pvalue <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$Chromosome <- paste0("chr", df.gwas.sig$Chromosome)



        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$Chromosome, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$Position,
                                                end = df.gwas.sig$Position,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to AD:\n")
    print(length(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig)))

    rm(gr.gwas.sig)
}

# ADHD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing ADHD data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "adhd/adhd_eur_jun2017.gz"), delim = "\t")

        df.gwas.sig <- filter(df.gwas, P <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to ADHD:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# BD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing BD data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "bd/daner_PGC_BIP32b_mds7a_0416a.gz"), delim = "\t")

        df.gwas.sig <- filter(df.gwas, P <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to BD:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Dep_symptoms #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Depressive Symptoms data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "depressive_symptoms/DS_Full.txt.gz"), delim = "\t")


        df.gwas.sig <- filter(df.gwas, Pval <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to Depressive Symptoms:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Edu_attain #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Educational Attainment data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "educational_attainment/GWAS_EA_excl23andMe.txt.gz"), delim = "\t")


        df.gwas.sig <- filter(df.gwas, Pval <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to Educational Attainment:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Epilepsy #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Epilepsy data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "epilepsy/all_epilepsy_METAL.gz"), delim = "\t")


        df.gwas.sig <- filter(df.gwas, `P-value` <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to Epilepsy:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Insomnia #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Insomnia data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "insomnia/Insomnia_sumstats_Jansenetal.txt.gz"), delim = " ")

        df.gwas.sig <- filter(df.gwas, P <= GWAS.PVAL.THRESHOLD)

        rm(df.gwas)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- HG19.SEQINFO

        rm(df.gwas.sig)
    }
    cat("Number of overlaps to Insomnia:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# IQ #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing IQ data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "iq/SavageJansen_2018_intelligence_metaanalysis.txt.gz"), delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to IQ:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Neuroticism #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Neuroticism data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "neuroticism/sumstats_neuroticism_ctg_format.txt.gz"), delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$RSID),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Neuroticism:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# PD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing PD data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "pd/PD_rsid_added.txt.gz"), delim = "\t")
        df.gwas.sig <- filter(df.gwas, p <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$chr)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$bp,
                                                end = df.gwas.sig$bp,
                                                names = df.gwas.sig$rs),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to PD:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Sub_well_being #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Subjective Well Being data...\n")
        df.gwas <- read_delim(paste0(dir.gwas, "subjective_well_being/SWB_Full.txt.gz"), delim = "\t")
        df.gwas.sig <- filter(df.gwas, Pval <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Subjective Well Being:\n")
    print(findOverlaps(gr.eqtl.snps.in.ld.hg19, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Anorexia_nervosa #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Anorexia Nervosa data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/anorexia_nervosa/pgc.ed.freeze1.July2017.zip", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Anorexia Nervosa:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Antisocial_behavior #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Antisocial Behavior data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/antisocial_behavior/BroadABC_METALoutput_CombinedwrsID.tbl.gz", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P.value <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$chr)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$bp,
                                                end = df.gwas.sig$bp,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Antisocial Behavior:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Anxiety_case_control #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Anxiety case/control data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Anxiety/anxiety.meta.full.cc.tbl.gz", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P.value <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNPID),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Anxiety case/control:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Anxiety_factor_score #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Anxiety factor score data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Anxiety/anxiety.meta.full.fs.tbl.gz", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P.value <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNPID),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Anxiety factor score:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# ASD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing ASD data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/ASD/iPSYCH-PGC_ASD_Nov2017.gz", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to ASD:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Internalizing #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Internalizing data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Internalizing/meta3.INTsnplist_F.txt", delim = " ")
        df.gwas.sig <- filter(df.gwas, meta.pval <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$chr)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$pos,
                                                end = df.gwas.sig$pos,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Internalizing:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Loneliness #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Loneliness data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Loneliness/Linear4PGC.txt", delim = " ")
        df.gwas.sig <- filter(df.gwas, PVALUE <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNPID),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Loneliness:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# MDD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing MDD data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/ENIGMA3/diseasetraitmunged/origdownloads/mdd/MDD2018_ex23andMe.gz", delim = "\t")
        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to MDD:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# SCZ #################
if (TRUE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing SCZ data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/ENIGMA3/diseasetraitmunged/origdownloads/clozuk/clozuk_pgc2.meta.sumstats.txt.gz", delim = " ")
        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        df.gwas.sig$CHR <- paste0("chr", df.gwas.sig$CHR)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        #seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to SCZ:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}
