# gviz colocalization scratch

library(here)
library(dplyr)
library(readxl)
library(readr)
library(magrittr)
library(Gviz)
library(GenomicRanges)
library(GenomicAlignments)
library(DESeq2)
library(biomaRt)
library(rtracklayer)

# OUTPUT ##############################################################################################################
#dir.pdf <- here("doc/plot_colocalization/pdfs/")

#dir.create(file.path(dir.pdf))

# INPUT ###############################################################################################################
# mirna ranged summarized experiment for rowranges
#rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# GRangesList containing SNP-miR associations within 1MB of each miR (only loaded if needed)
grangeslist.rds <- here("results/rdata_files/20190808_mirQTL.GRangesList.rds")
# Dataframe containing the GRangesList data (computation is faster on dataframe than the GRangesList)
dataframe.rds <- here("results/rdata_files/20190808_mirQTL.dataframe.rds")

# list of emiRs (eGenes)
#emiRs.tsv <- here("results/emmax/association_results/20190808_mirQTL/emiRs.tsv")
# nominal p-value threshold
nom.pval.thresh.txt <- here("results/emmax/association_results/20190808_mirQTL/nominal.pval.threshold")

# directory of clumped results
#clump.dir <- here("results/emmax/association_results/20190808_mirQTL/clumped/")

# GLOBALS ###########################3
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import SNP Associations #############################################################################################
# if the convertered dataframe is available, load it
# else import grangeslist and convert to dataframe
if (! exists("gr.snp.mir.sig") ) {
    if (! exists("df.snp.mir.assoc") ) {
        if (file.exists(dataframe.rds)) {
            df.snp.mir.assoc <- readRDS(dataframe.rds)
        } else {
            grl <- readRDS(grangeslist.rds)

            # unlist GRangesList into GRanges
            gr.snps <- unlist(grl)
            # split name of each range into snpid and mir
            gr.snps$snpid <- sapply(strsplit(names(gr.snps), "\\."), `[`, 2)
            gr.snps$mir <- sapply(strsplit(names(gr.snps), "\\."), `[`, 1)
            # remove names from GRanges
            names(gr.snps) <- NULL

            # convert GRanges into dataframe
            df.snp.mir.assoc <- as.data.frame(gr.snps)

            # select and rename columns
            df.snp.mir.assoc %<>%
                select(chr = seqnames,
                       coord.hg19 = start,
                       snpid,
                       mir, A1, A2.effect, A1.HOM.count, HET.count, A2.HOM.count, beta, pval)

            # save dataframe
            saveRDS(df.snp.mir.assoc, dataframe.rds)

            # clean up
            rm(grl, gr.snps)
        }
    }

    # filter snp-mir associations to not be dependent on only 1 homozygous minor sample
    df.snp.mir.assoc %<>% filter(A1.HOM.count != 1, HET.count >= 2)

    nom.pval.thresh <- as.numeric(read_lines(nom.pval.thresh.txt))

    # filter for only SNPs below nominal p-value threshold
    df.snp.mir.sig <- filter(df.snp.mir.assoc, pval <= nom.pval.thresh)

    # GRanges of significant snps
    gr.snp.mir.sig <- GRanges(seqnames = factor(df.snp.mir.sig$chr, levels = CHROMS),
                              ranges = IRanges(start = df.snp.mir.sig$coord.hg19,
                                               end = df.snp.mir.sig$coord.hg19,
                                               names = df.snp.mir.sig$snpid),
                              strand = "*")

    seqinfo(gr.snp.mir.sig) <- Seqinfo(genome = "hg19")
    rm(df.snp.mir.assoc)
}

# AD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        df.gwas1 <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/AD/Kunkle_etal_Stage1_results.txt?file=1", delim = " ")
        df.gwas2 <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/AD/Kunkle_etal_Stage2_results.txt?file=1", delim = " ")
        df.gwas <- bind_rows(df.gwas1, df.gwas2)
        df.gwas$Chromosome <- paste0("chr", df.gwas$Chromosome)

        df.gwas.sig <- filter(df.gwas, Pvalue <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$Chromosome, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$Position,
                                                end = df.gwas.sig$Position,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas1, df.gwas2, df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to AD:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# ADHD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing ADHD data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/ENIGMA3/diseasetraitmunged/origdownloads/ADHD/adhd_eur_jun2017.gz", delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to ADHD:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# BD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing BD data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/BD/daner_PGC_BIP32b_mds7a_0416a.gz", delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to BD:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Dep_symptoms #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Depressive Symptoms data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Depressive_symptoms/DS_Full.txt.gz", delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, Pval <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Depressive Symptoms:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Edu_attain #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Educational Attainment data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Educational_attainment/GWAS_EA_excl23andMe.txt.gz", delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, Pval <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$POS,
                                                end = df.gwas.sig$POS,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Educational Attainment:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Epilepsy #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Epilepsy data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Epilepsy/all_epilepsy_METAL.gz", delim = "\t")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, `P-value` <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$MarkerName),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Epilepsy:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Insomnia #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Insomnia data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Insomnia/Insomnia_sumstats_Jansenetal.txt.gz", delim = " ")
        df.gwas$CHR <- paste0("chr", df.gwas$CHR)

        df.gwas.sig <- filter(df.gwas, P <= 5e-8)

        gr.gwas.sig <- GRanges(seqnames = factor(df.gwas.sig$CHR, levels = CHROMS),
                               ranges = IRanges(start = df.gwas.sig$BP,
                                                end = df.gwas.sig$BP,
                                                names = df.gwas.sig$SNP),
                               strand = "*")

        seqinfo(gr.gwas.sig) <- Seqinfo(genome = "hg19")

        rm(df.gwas, df.gwas.sig)
    }
    cat("Number of overlaps to Insomnia:\n")
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# IQ #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing IQ data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/IQ/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt.gz", delim = "\t")
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
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Neuroticism #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Neuroticism data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Neuroticism/sumstats_neuroticism_ctg_format.txt.gz", delim = "\t")
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
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# PD #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing PD data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/PD/PD_rsid_added.txt.gz", delim = "\t")
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
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

    rm(gr.gwas.sig)
}

# Sub_well_being #################
if (FALSE) {
    if (! exists("gr.gwas.sig")) {
        cat("Importing Subjective Well Being data...\n")
        df.gwas <- read_delim("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/Subjective_well-being/SWB_Full.txt.gz", delim = "\t")
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
    print(findOverlaps(gr.snp.mir.sig, gr.gwas.sig))

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
