# find co-localization miRNA-eQTL to GWAS data

library(here)
library(dplyr)
library(readr)
library(magrittr)
library(GenomicRanges)
library(mikelaffr)

# OUTPUT FILES #########################################################################################################
gwas.dataset <- "scz"

# working directory for temp plink files
dir.working <- paste0(here("results/external_data/"), gwas.dataset, "/")
dir.create(dir.working, recursive = TRUE, showWarnings = FALSE)

# table of colocalizations
coloc.output.rds <- paste0(dir.working, gwas.dataset, "_mirQTL_colocalizations.rds")

# INPUT FILES ##########################################################################################################
# GWAS data
gr.gwas.rds <- paste0(here("data/gwas_datasets/"), gwas.dataset, "/", gwas.dataset, ".hg38.GRanges.rds")

# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_dataFrame.rds")
eqtls.granges.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_eQTLs_GRanges.rds")

# mirQTL all variants
mirqtl.results.variants.dataframe.rds <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_variants_dataFrame.rds")

# Directory for LD at each mirQTL index SNP
dir.ld <- here("results/emmax/association_results/20200120_mirQTLor/ld/")

# 1000 Genomes EUR plink files
dir.1kg.eur.plink <- here("data/1000genomes_phase3_hg38/EUR.plink/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

GWAS.PVAL.THRESHOLD <- 5e-8

# FUNCTIONS ############################################################################################################
getRSIDfrom1kgEUR <- function(chr, pos) {

    # check for correct chr format
    if (! chr %in% CHROMS) {
        stop(paste("chromosome", chr, "not in correct format!"))
    }

    # load .bim plink file
    df.bim <- read_tsv(paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38.bim"), col_names = c("chr", "rsid", "cm", "pos", "a1", "a2"))

    # return corresponding RSID
    rsid <- df.bim$rsid[match(pos, df.bim$pos)]

    rm(df.bim)

    return(rsid)
}

# Import ###############################################################################################################

# mirQTL eqtls
df.mirQTL.eqtls <- readRDS(eqtls.dataframe.rds)

# import gwas data
gr.gwas <- readRDS(gr.gwas.rds)

# Find Co-Localizations ################################################################################################

df.coloc <- tibble()

# loop over each mirQTL
for (i in 1:nrow(df.mirQTL.eqtls)) {

    printMessage(paste(i, "of", nrow(df.mirQTL.eqtls)))

    eqtl <- df.mirQTL.eqtls$eQTL[i]
    emir <- df.mirQTL.eqtls$emiR[i]
    esnp <- df.mirQTL.eqtls$eSNP[i]

    chr <- df.mirQTL.eqtls$SNP.CHR[i]
    pos <- df.mirQTL.eqtls$SNP.BP.hg38[i]

    # look for GWAS genome wide significant snps within 1MB of index SNP
    gr.gwas.window <- subsetByOverlaps(gr.gwas, ranges = GRanges(seqnames = chr,
                                                                 ranges = IRanges(start = pos - 1e6,
                                                                                  end = pos + 1e6)))
    if (!any(gr.gwas.window$Pval <= GWAS.PVAL.THRESHOLD)) {
        print("No genome wide significant snps within 1MB of index SNP")
        next
    }

    # use PLINK to clump data in 1MB window
    # 1000Genomes EUR bfile
    bfile <- paste0(dir.1kg.eur.plink, "EUR.", chr, "_GRCh38")
    # results file for clumping
    clumpFile <- paste0(dir.working, "tmp.results.window.for.clump.tsv")
    # output file prefix
    outputPrefix <- paste0(dir.working, "tmp.clump")
    # write clumpFile
    write_tsv(tibble(SNP = gr.gwas.window$SNPID, P = gr.gwas.window$Pval), clumpFile)
    # plink command
    command <- sprintf("plink --bfile %s --clump %s --clump-kb 250 --clump-p1 %s --clump-p2 0.01 --clump-r2 0.2 --out %s",
                       bfile, clumpFile, GWAS.PVAL.THRESHOLD, outputPrefix)
    # call PLINK
    system(command = command)

    # list of GWAS index snps in this window
    gwas.index.snps <- NULL
    command <- sprintf("awk '{if (NR!=1) {print $3}}' %s | awk NF", paste0(dir.working, "tmp.clump.clumped"))
    gwas.index.snps <- system(command = command, intern = TRUE)

    # check for at least one index snp
    if (length(gwas.index.snps) == 0) {
        stop("No GWAS index SNPs found?!")
    }

    # calculate LD buddies for GWAS index SNPs
    for (j in 1:length(gwas.index.snps)) {

        outputPrefix <- paste0(dir.working, "tmp.", gwas.index.snps[j])

        command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 2000 --ld-window-kb 1000 --ld-window-r2 0.2 --out %s",
                           bfile, gwas.index.snps[j], outputPrefix)

        # call PLINK
        system(command = command)
    }

    # get LD buddies for mirQTL index SNP
    # ld file: chrX.hardcall.prefiltered.mirQTLor.chrX:49189007:G:A.ld
    ld.file <- paste0(dir.ld, chr, ".hardcall.prefiltered.mirQTLor.", esnp, ".ld")

    # ld for this index SNP
    df.mirQTL.ldbuddies <- read_table(ld.file, col_types = cols())

    # get RSID from 1kg EUR dataset
    esnp.rsid <- getRSIDfrom1kgEUR(chr = chr,
                                   pos = pos)

    # check for NA rsid
    if (is.na(esnp.rsid)) {
        print("No 1KG EUR RSID found for index SNP!")
    } else {
        # get 1kg EUR ld to the mirQTL index SNP
        outputPrefix <- paste0(dir.working, "tmp.", esnp.rsid)
        command <- sprintf("plink --bfile %s --r2 --ld-snp %s --ld-window 2000 --ld-window-kb 1000 --ld-window-r2 0.2 --out %s",
                           bfile, esnp.rsid, outputPrefix)
        # call PLINK
        system(command = command)

        tmp <- read_table(paste0(dir.working, "tmp.", esnp.rsid, ".ld"))

        # combine mirQTL LD buddies with 1kg EUR LD buddies
        df.mirQTL.ldbuddies <- bind_rows(df.mirQTL.ldbuddies, tmp)

        rm(tmp)
    }

    # filter for only high LD buddies (r2 >= 0.8)
    df.mirQTL.ldbuddies.high <- dplyr::filter(df.mirQTL.ldbuddies, R2 >= 0.8)
    # filter for distinct positions
    df.mirQTL.ldbuddies.high <- dplyr::filter(df.mirQTL.ldbuddies.high, !duplicated(BP_B))

    # loop over GWAS index SNPs and look for co-localizations to GWAS index LD buddies
    for (j in 1:length(gwas.index.snps)) {

        # gwas index snp ld file
        gwas.index.ld.file <- paste0(dir.working, "tmp.", gwas.index.snps[j], ".ld")

        # read in ld buddies
        df.gwas.ldbuddies <- read_table(gwas.index.ld.file, col_types = cols())

        # retain high ld buddies
        df.gwas.ldbuddies %<>%
            dplyr::filter(R2 >= 0.8)

        # do any gwas ld buddies overlap with mirQTL ld buddies?
        if (any(df.mirQTL.ldbuddies.high$BP_B %in% df.gwas.ldbuddies$BP_B)) {
            df.overlap <- inner_join(df.mirQTL.ldbuddies.high, df.gwas.ldbuddies, by = "BP_B", suffix = c(".mirQTL", ".gwas"))
            df.overlap$eQTL <- eqtl
            df.overlap$emiR <- emir
            df.overlap$eSNP <- esnp
            df.coloc <- bind_rows(df.coloc, df.overlap)

            printMessage("Colocalization Found!", fillChar = "!")
        }
    }

    # remove all plink tmp files
    system(command = paste0("rm ", dir.working, "tmp.*"))

    print(paste0("Total ", gwas.dataset, " colocalizations found so far: ", nrow(df.coloc)))
}

print("Saving colocalizations...")

if (nrow(df.coloc) > 0) {
    saveRDS(df.coloc, coloc.output.rds)
}

print("Finished!")




