# Create phenotype files of miRNA expression per sample for EMMAX analysis by primary emiR
# pre-regress out the covariate matrix and primary eQTL genotype

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)
library(DESeq2)

# date.prefix <- format(Sys.time(), "%Y%m%d")
# date.prefix <- "20200120"

# association results directory name
results.name <- "20200120_mirQTLor_fdr10percent"

# OUTPUT FILES #########################################################################################################
# output directory
output.dir <- paste0(here("results/conditional_eqtls/"), results.name, "/tertiary/", results.name, "_VST_miRNA_expression_residual_after_tertiary_eqtl/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# gene list and coords output file
genes.hg38.output.txt <- paste0(output.dir, results.name, "_genes_hg38.txt")

# variance stabilized data set used for phenotype files exported as a ranged summarized experiment
vsd.output.rds <- paste0(output.dir, results.name, "_VST_miRNA_expression_and_residual_after_tertiary_eqtl_rse.rds")

# INPUT FILES ##########################################################################################################
# samples file
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# non-overlapping granges of known and novel mirnas
granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Autosomes covariate matrix
cov.mat.auto.file <- here("results/emmax/covariates/20200120_chrAutosomes.mirQTLor.columnLabeled.cov")
# Chromosome X covariate matrix
# chrX now treated the same way as autosomes!

# data frame of primary eqtls
primary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/primary/20200120_mirQTLor_fdr10percent_primary_eQTLs_dataFrame.rds")

# data frame of secondary eqtls
secondary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/secondary/20200120_mirQTLor_fdr10percent_secondary_eQTLs_dataFrame.rds")

# data frame of tertiary eqtls
tertiary.eqtls.rds <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/tertiary/20200120_mirQTLor_fdr10percent_tertiary_eQTLs_dataFrame.rds")

# directory with primary eqtl genotypes
dir.primary.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/primary/primary_eSNP_genotypes/")

# directory with secondary eqtl genotypes
dir.secondary.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/secondary/secondary_eSNP_genotypes/")

# directory with tertiary eqtl genotypes
dir.tertiary.genotypes <- here("results/conditional_eqtls/20200120_mirQTLor_fdr10percent/tertiary/tertiary_eSNP_genotypes/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX", "chrXPAR")

PAR.LOWER.BP <- 2781479
PAR.UPPER.BP <- 155701383
# Import Samples #######################################################################################################
samples <- read_delim(samples.txt, col_names = c("RNAID", "DonorID", "DNAID"), delim = " ")

# Import RSE ###########################################################################################################
rse <- readRDS(rse.rds)

# threshold for expression
rse <- rse[rowSums(assay(rse) >= 10) >= 10, ]

# convert to DESeqDataSet for vst normalization
dds <- DESeqDataSet(rse, design = ~1)

vsd <- varianceStabilizingTransformation(dds)

# subset samples
vsd <- vsd[,samples$RNAID]

rm(rse, dds)

# Import GRanges #######################################################################################################
gr <- readRDS(granges.rds)

# Format Expression Table ##############################################################################################
# remove rows (miRNAs) that are not in granges (these have been removed as overlapping another annotation)
vsd <- vsd[rowData(vsd)$Name %in% gr$Name,]

# expression data
df.exprs <- data.frame(Name = rowData(vsd)$Name, assay(vsd), stringsAsFactors = FALSE)
# duplicate from 3p/5p in mirge but not in mirbase
#df.exprs$Name[duplicated(df.exprs$Name)]
# hsa-miR-4473
# modify Name variable to include 5p/3p
df.exprs[which(df.exprs$Name == "hsa-miR-4473")[1], "Name"] <- rownames(df.exprs)[which(df.exprs$Name == "hsa-miR-4473")[1]]
df.exprs[which(df.exprs$Name == "hsa-miR-4473")[1], "Name"] <- rownames(df.exprs)[which(df.exprs$Name == "hsa-miR-4473")[1]]
# should have zero duplicates
stopifnot(sum(duplicated(df.exprs$Name)) == 0)

# Format miRNA Table ###################################################################################################
# granges data
df.mirs <- as.data.frame(gr)

# add row and modify existing to have 5p/3p for miR-4437
df.tmp <- df.mirs[df.mirs$Name == "hsa-miR-4473",]
df.tmp[1,"Name"] <- paste0(df.tmp[1,"Name"], "-3p")
df.tmp[1,"UniqueName"] <- paste(df.tmp[1,"Name"], df.tmp[1,"Derives_from_name"], sep = "_")

df.mirs[df.mirs$Name == "hsa-miR-4473","Name"] <- paste0(df.mirs[df.mirs$Name == "hsa-miR-4473","Name"], "-5p")
df.mirs[df.mirs$Name == "hsa-miR-4473-5p","UniqueName"] <- paste(df.mirs[df.mirs$Name == "hsa-miR-4473-5p","Name"], df.mirs[df.mirs$Name == "hsa-miR-4473-5p","Derives_from_name"], sep = "_")
df.mirs <- bind_rows(df.mirs, df.tmp)
rm(df.tmp)

# all Names in expression table should be names in mirna table
# there should be no duplicate names in expression table
# there should be no duplicate uniqueNames in the mir table
stopifnot(all(df.exprs$Name %in% df.mirs$Name))
stopifnot(sum(duplicated(df.exprs$Name)) == 0)
stopifnot(sum(duplicated(df.mirs$UniqueName)) == 0)

# join expression table with mirna table, expanding expression table to fill repeated names
# remove mirnas not in expression table
df.mirs <- df.mirs[df.mirs$Name %in% df.exprs$Name,]
df.joined <- left_join(df.mirs, df.exprs, by = "Name")

#sum(duplicated(df.joined$Name))
stopifnot(sum(duplicated(df.joined$UniqueName)) == 0)

# remove non-standard chromosomes
df.joined %<>%
    filter(seqnames %in% CHROMS)

# Import Primary eQTLs #################################################################################################
# primary eqtls
df.primary.eqtls <- read_rds(primary.eqtls.rds)

# secondary eqtls
df.secondary.eqtls <- read_rds(secondary.eqtls.rds)

# tertiary eqtls
df.tertiary.eqtls <- read_rds(tertiary.eqtls.rds)

# filter df.joined for secondary eqtl emiRs
df.joined %<>%
    filter(UniqueName %in% df.tertiary.eqtls$UniName)

# Import Covariate Matrices ############################################################################################
# autosome covariates file
cov.mat.auto <- read_delim(cov.mat.auto.file, delim = " ")
cov.mat.auto %<>%
    dplyr::rename(DonorID = X1, DNAID = X2)

# covariates
batchVars.auto <- c(colnames(cov.mat.auto)[4:length(colnames(cov.mat.auto))], "primaryDosage", "secondaryDosage", "tertiaryDosage")

# formula strings
formula.string.auto <- paste("expression ~", paste(batchVars.auto, collapse = " + "))

# Export Gene Lists ####################################################################################################

# # export gene list and coords hg19
# write_tsv(dplyr::select(df.lo, uniqename = mcols.uniqename, seqnames, start, end), genes.hg19.output.txt, col_names = FALSE)

# export gene list and coords hg38
write_tsv(dplyr::select(df.joined, UniqueName, seqnames, start, end), genes.hg38.output.txt, col_names = FALSE)

# Export Phenotype Files ###############################################################################################

# new expression table
df.exprs <- dplyr::select(df.joined, UniqueName, starts_with("RNAID"))
mat.exprs <- as.matrix(df.exprs[,2:213])
rownames(mat.exprs) <- df.exprs$UniqueName

# empty data frame for residual expression
df.residual <- data.frame(RNAID = samples$RNAID, stringsAsFactors = FALSE)

# loop over each chromosome

numMirs <- nrow(mat.exprs)
numFiles <- 0

print(paste0("Looping over ", numMirs, " miRNAs..."))

for (chr in CHROMS) {
    print(paste("On ", chr, "......", sep=""))

    # create output folder
    path <- file.path(output.dir, chr)
    dir.create(path, recursive = TRUE, showWarnings = FALSE)

    # mirnas on this chrom
    mirs <- NULL
    mirs <- df.joined$UniqueName[which(df.joined$seqnames == chr)]

    # loop over each mirna
    for (mir in mirs) {

        #printMessage(paste(i, mir, sep = ": "))
        print(mir)

        # reset data for safety
        filename <- NULL
        df <- NULL
        output <- NULL
        formula.string <- NULL

        # filename
        filename <- file.path(path, paste0(chr, ".mirQTLor.", mir, ".tertiary.pheno"))

        # expression matrix
        df <- data.frame(expression = mat.exprs[mir,], stringsAsFactors = FALSE)
        df$RNAID <- rownames(df)

        # join with samples to get DNAID
        df %<>%
            left_join(samples, by = "RNAID")

        # import primary sample genotypes
        esnp.primary <- df.primary.eqtls$SNP[which(df.primary.eqtls$UniName == mir)]
        df.primary.genotypes <- read_table2(paste0(dir.primary.genotypes, esnp.primary, ".raw"))
        df.primary.genotypes %<>%
            select(DNAID = IID, primaryDosage = 7)

        # import secondary sample genotypes
        esnp.secondary <- df.secondary.eqtls$SNP[which(df.secondary.eqtls$UniName == mir)]
        df.secondary.genotypes <- read_table2(paste0(dir.secondary.genotypes, esnp.secondary, ".raw"))
        df.secondary.genotypes %<>%
            select(DNAID = IID, secondaryDosage = 7)

        # import tertiary sample genotypes
        esnp.tertiary <- df.tertiary.eqtls$SNP[which(df.tertiary.eqtls$UniName == mir)]
        df.tertiary.genotypes <- read_table2(paste0(dir.tertiary.genotypes, esnp.tertiary, ".raw"))
        df.tertiary.genotypes %<>%
            select(DNAID = IID, tertiaryDosage = 7)

        # join with cov matrix
        df %<>%
            left_join(dplyr::select(cov.mat.auto, -DonorID), by = "DNAID") %>%
            left_join(df.primary.genotypes, by = "DNAID") %>%
            left_join(df.secondary.genotypes, by = "DNAID") %>%
            left_join(df.tertiary.genotypes, by = "DNAID")

        formula.string <- formula.string.auto

        # set rownames
        rownames(df) <- df$DonorID
        # fit model
        fit <- lm(formula = as.formula(formula.string), data = df)

        # expression residuals
        df.resid <- data.frame(expression.residual = resid(fit),
                                   DonorID = names(resid(fit)),
                                   stringsAsFactors = FALSE)

        df %<>%
            left_join(df.resid, by = "DonorID")

        # save residuals in data frame
        df.residual %<>%
            left_join(dplyr::select(df, RNAID, !!mir := expression.residual), by = "RNAID")

        # create output data frame
        output <- dplyr::select(df, DonorID, DNAID, expression.residual)

        stopifnot(all(output$DonorID == samples$DonorID))

        # write file
        write_delim(output, filename, delim = " ", col_names = FALSE)
        numFiles <- numFiles + 1
    }
}

print(paste0("Finished writing ", numFiles, " files."))

# Export VSD ###########################################################################################################
# modified VSD which has correct rowranges and expression for mirQTLor samples and mirnas

# format expression matrix
df.assay <- select(df.joined, starts_with("RNAID"))
rownames(df.assay) <- df.joined$UniqueName

# format residual matrix
mat.residual <- as.matrix(df.residual[,-1])
rownames(mat.residual) <- df.residual$RNAID
colnames(mat.residual) <- colnames(df.residual)[2:ncol(df.residual)]
mat.residual <- t(mat.residual)
mat.residual <- mat.residual[rownames(df.assay), colnames(df.assay)]


# format rowRanges
gr.mirs <- GRanges(select(df.joined, -starts_with("RNAID")))
names(gr.mirs) <- gr.mirs$uniqueName

all(rownames(df.assay) == rownames(mat.residual))
all(colnames(df.assay) == colnames(mat.residual))

all(names(gr.mirs) == rownames(df.assay))

# construct Ranged Summarized Experiment using colData from vsd
rse.vsd.mirQTLor <- SummarizedExperiment(assays = list("VST" = as.matrix(df.assay), "VST.RESIDUAL" = mat.residual),
                                         rowRanges = gr.mirs,
                                         colData = colData(vsd))

saveRDS(rse.vsd.mirQTLor, vsd.output.rds)

