# Create phenotype files of miRNA expression per sample for EMMAX analysis by CHR
# pre-regress out the covariate matrix created in the 20191204 run

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(mikelaffr)
library(DESeq2)

date.prefix <- format(Sys.time(), "%Y%m%d")
#date.prefix <- "20200120"

# OUTPUT FILES #########################################################################################################
# Phenotype file output directory
output.dir <- paste0(here("results/emmax/phenotype_files/"), date.prefix, "_mirQTLor_VST_miRNA_expression_residual_4707_edu_attain/")

dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# gene list and coords output file
#genes.hg19.output.txt <- paste0(output.dir, date.prefix, "_mirQTL_genes_hg19.txt")
genes.hg38.output.txt <- paste0(output.dir, date.prefix, "_mirQTLor_genes_hg38.txt")

# variance stabilized data set used for phenotype files exported as a ranged summarized experiment
vsd.output.rds <- paste0(output.dir, date.prefix, "_mirQTLor_VST_miRNA_expression_and_residual_rse.rds")

# INPUT FILES ##########################################################################################################
# samples file
samples.txt <- here("results/emmax/samples/20200120_mirQTLor_RNAID_DonorID_DNAID.txt")

# ranged summarized experiment for expression values
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# non-overlapping granges of known and novel mirnas
granges.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# Autosomes covariate matrix
cov.mat.auto.file <- here("results/emmax/covariates/20201012_chrAutosomes.mirQTLor.columnLabeled.edu_index.cov")
# Chromosome X covariate matrix
#cov.mat.x.file <- here("results/emmax/covariates/20200120_chrX.mirQTLor.columnLabeled.cov")

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

# Modify chrX for PAR and nonPAR #######################################################################################
# indecies of mirs on chrX
ind.mirs.chrX <- which(df.joined$seqnames == "chrX")

df.joined$seqnames <- as.character(df.joined$seqnames)

# loop over mirs on chrX and change seqnames to
for (i in ind.mirs.chrX) {

    #print(df.joined$uniqueName[i])

    # check that no mir is within 1MB of the PAR boundaries
    if (df.joined$end[i] > PAR.LOWER.BP - 1e6 & df.joined$end[i] < PAR.LOWER.BP + 1e6) {
        print("mir END near LOWER Boundary!")
    }
    if (df.joined$start[i] > PAR.LOWER.BP - 1e6 & df.joined$start[i] < PAR.LOWER.BP + 1e6) {
        print("mir START near LOWER Boundary!")
    }

    if (df.joined$end[i] > PAR.UPPER.BP - 1e6 & df.joined$end[i] < PAR.UPPER.BP + 1e6) {
        print("mir END near UPPER Boundary!")
    }
    if (df.joined$start[i] > PAR.UPPER.BP - 1e6 & df.joined$start[i] < PAR.UPPER.BP + 1e6) {
        print("mir START near UPPER Boundary!")
    }

    # modify seqnames
    if (df.joined$start[i] < PAR.LOWER.BP | df.joined$start[i] > PAR.UPPER.BP) {
        df.joined$seqnames[i] <- "chrXPAR"
    }
}

df.joined$seqnames <- factor(df.joined$seqnames)

# Import Covariate Matrices ############################################################################################
# autosome covariates file
cov.mat.auto <- read_delim(cov.mat.auto.file, delim = " ")
cov.mat.auto %<>%
    dplyr::rename(DonorID = X1, DNAID = X2)

# covariates
batchVars.auto <- c(paste0("PC", 1:10, ".genotype"),
                    paste0("PC", 1:10, ".expression"),
                    "PoolPool2", "PoolPool3", "PoolPool4",
                    "PoolPool5", "PoolPool6", "PoolPool7",
                    "PoolPool8","PurificationMethodmiRNeasy",
                    "PurificationMethodmiRNeasy_mini", "SexM", "RIN", "GestationWeek", "IndexGenos")

# formula strings
formula.string.auto <- paste("expression ~", paste(batchVars.auto, collapse = " + "))

# x chrom covariates file
cov.mat.x <- read_delim(cov.mat.x.file, delim = " ")
cov.mat.x %<>%
    dplyr::rename(DonorID = X1, DNAID = X2)

# covariates
batchVars.x <- c(paste0("PC", 1:10, ".genotype"),
                 paste0("PC", 1:10, ".expression"),
                 "PoolPool2", "PoolPool3", "PoolPool4",
                 "PoolPool5", "PoolPool6", "PoolPool7",
                 "PoolPool8","PurificationMethodmiRNeasy",
                 "PurificationMethodmiRNeasy_mini", "RIN", "GestationWeek")

# formula strings
formula.string.x <- paste("expression ~", paste(batchVars.x, collapse = " + "))

# Export Gene Lists ####################################################################################################

# # export gene list and coords hg19
# write_tsv(dplyr::select(df.lo, uniqename = mcols.uniqename, seqnames, start, end), genes.hg19.output.txt, col_names = FALSE)

# export gene list and coords hg38
df.joined <- df.joined[df.joined$UniqueName == "hsa-mir-4707_hsa-miR-4707-3p" | df.joined$UniqueName == "hsa-mir-342_hsa-miR-342-5p",]

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
        filename <- file.path(path, paste0(chr, ".mirQTLor.", mir, ".pheno"))

        # expression matrix
        df <- data.frame(expression = mat.exprs[mir,], stringsAsFactors = FALSE)
        df$RNAID <- rownames(df)

        # join with samples to get DNAID
        df %<>%
            left_join(samples, by = "RNAID")

        # join with cov matrix
        if (chr == "chrX") {
            df %<>%
                left_join(dplyr::select(cov.mat.x, -DonorID), by = "DNAID")
            formula.string <- formula.string.x
        } else {
            df %<>%
                left_join(dplyr::select(cov.mat.auto, -DonorID), by = "DNAID")
            formula.string <- formula.string.auto
        }

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
rownames(df.assay) <- df.joined$uniqueName

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

