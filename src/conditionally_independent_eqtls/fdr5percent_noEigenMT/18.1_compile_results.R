# Compile results after EMMAX analysis to save as a GRanges Rdata object and flat text file for compression and
# archiving.

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(GenomicRanges)
# library(rtracklayer)
# library(liftOver)
library(mikelaffr)

#date.prefix <- format(Sys.time(), "%Y%m%d")
#date.prefix <- "20200120"

# association results directory name
results.name <- "20200120_mirQTLor_fdr5percent_noEigenMT"

# OUTPUT FILES ########################################################################################################
# directory for compiled results
output.dir <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/quarternary/association_results/compiled/")
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)

# output GRangesList RDS file
output.grl.rds <- paste0(output.dir, results.name, "_quarternary_variants_GRangesList.rds")

# output data frame RDS file
output.df.rds <- paste0(output.dir, results.name, "_quarternary_variants_dataFrame.rds")

# INPUT FILES #########################################################################################################
# emmax result files directory, column labeled
results.dir <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/quarternary/association_results/ps/")

# mirs in this analysis
genes.tsv <- here("results/conditional_eqtls/20200120_mirQTLor_fdr5percent_noEigenMT/tertiary/20200120_mirQTLor_fdr5percent_noEigenMT_VST_miRNA_expression_residual_after_tertiary_eqtl/20200120_mirQTLor_fdr5percent_noEigenMT_genes_hg38.txt")

# dosage genotypes
tfile.dosage.dir.root <- here("results/emmax/tfiles/prefiltered_mirQTLor_dosages/")

# hardcall genotypes
tfile.hardcall.dir.root <- here("results/emmax/tfiles/prefiltered_mirQTLor_hardcalls/")

# eigenMT-BH nominal p-value
nominal.p.value.txt <- here("results/emmax/association_results/20200120_mirQTLor/compiled/20200120_mirQTLor_nomPvalue.txt")

# GLOBALS #############################################################################################################
FDR.THRESHOLD <- as.numeric(read_lines(nominal.p.value.txt))

# Import Data #########################################################################################################

# load genes used in this anaysis
genes <- read_tsv(genes.tsv, col_names = c("uniqueName", "chr", "start", "end"))

# an empty list for GRanges objects, one for each mir
grl.results <- GRangesList()

# an empty data frame for results
df.results <- data.frame()

# number of results files loaded/compiled
nres <- 0

# loop over each gene
#for (i in 1:5) {
for (i in 1:nrow(genes)) {
    printMessage(paste("Importing data for", genes$uniqueName[i]), fillChar = "+", justify = "m")

    uniName <- NULL
    chr <- NULL
    ps.file <- NULL
    pvar.file <- NULL
    acount.file <- NULL
    ps <- NULL
    pvar <- NULL
    acount <- NULL
    snps <- NULL
    df.res <- NULL
    gr.res <- NULL

    uniName <- genes$uniqueName[i]
    chr <- genes$chr[i]

    ps.file <- paste0(results.dir, chr, ".dosage.prefiltered.mirQTLor.quarternary.", uniName, ".ps")
    pvar.file <- paste0(tfile.dosage.dir.root, chr, "/", chr, ".dosage.prefiltered.mirQTLor.", uniName, ".pvar")
    acount.file <- paste0(tfile.dosage.dir.root, chr, "/", chr, ".dosage.prefiltered.mirQTLor.", uniName, ".acount")
    traw.file <- paste0(tfile.dosage.dir.root, chr, "/", chr, ".dosage.prefiltered.mirQTLor.", uniName, ".traw")

    frqx.file <- paste0(tfile.hardcall.dir.root, chr, "/", chr, ".hardcall.prefiltered.mirQTLor.", uniName, ".frqx")


    # check if results file exists
    if (!file.exists(ps.file)) {
        print(paste("No results for", uniName))
        next
    }

    # import results file
    ps <- read_tsv(ps.file)
    ps$UniName <- uniName
    # import pvar file
    pvar <- read_tsv(pvar.file,
                     col_names = c("CHR", "BP.hg38", "SNP", "REF", "ALT", "FILTER", "INFO"),
                     skip = 17,
                     col_types = "ciccccc")
    # import acount file
    acount <- read_tsv(acount.file,
                       col_names = c("CHR", "SNP", "REF", "ALT", "ALT_CTS", "OBS_CT"),
                       skip = 1,
                       col_types = "ccccdd")
    # import traw file
    traw <- read_tsv(traw.file, col_types = cols_only(SNP = col_character(), COUNTED = col_character()))
    # import frqx file
    frqx <- read_tsv(frqx.file,
                     col_names = c("CHR", "SNP", "A1", "A2", "A1.HOM.count", "HET.count", "A2.HOM.count", "A1.HAP.count", "A2.HAP.count", "MISSING.count"),
                     skip = 1,
                     col_types = "cccciiiiii")

    # combine pvar and acount tables
    stopifnot(all(pvar$SNP == acount$SNP))
    snps <- bind_cols(select(pvar, -FILTER, -INFO), select(acount, -SNP, -REF, -ALT, -CHR))
    # combine snps with traw to get counted allele
    stopifnot(all(snps$SNP == traw$SNP))
    snps <- bind_cols(snps, select(traw, EFFECT.ALLELE = COUNTED))
    # combine with frqx file
    snps <- left_join(snps, select(frqx, -CHR, -A1.HAP.count, -A2.HAP.count, -MISSING.count), by = "SNP")
    # change to UCSC chr naming
    snps$CHR <- paste0("chr", snps$CHR)

    # results for this mir
    df.res <- left_join(ps, snps, by = "SNP")

    # compile into total results table
    df.results <- bind_rows(df.results, df.res)

    # make GRanges
    gr.res <- GRanges(seqnames = df.res$CHR,
                      ranges = IRanges(start = df.res$BP.hg38,
                                       end = df.res$BP.hg38,
                                       names = df.res$SNP),
                      strand = "*",
                      REF = df.res$REF,
                      ALT = df.res$ALT,
                      EFFECT.ALLELE = df.res$EFFECT.ALLELE,
                      ALT_CTS = df.res$ALT_CTS,
                      OBS_CT = df.res$OBS_CT,
                      A1 = df.res$A1,
                      A2 = df.res$A2,
                      A1.HOM.count = df.res$A1.HOM.count,
                      HET.count = df.res$HET.count,
                      A2.HOM.count = df.res$A2.HOM.count,
                      BETA = df.res$BETA,
                      PVAL = df.res$P,
                      SNP = df.res$SNP,
                      UniName = df.res$UniName)

    # add GRanges for this mirna to list
    grl.results[[uniName]] <- gr.res

    nres <- nres + 1
}

# report number of results compiled
print(paste(nres, "results files compiled."))

# set genome of GRangesList
#seqlevelsStyle(grl.results) <- "UCSC"
#seqinfo(grl.results) <- Seqinfo(genome = "hg38")

# write RDS files
print(paste("Saving GRangesList object:", output.grl.rds))
saveRDS(grl.results, output.grl.rds)

print(paste("Saving Data Frame object:", output.df.rds))
saveRDS(df.results, output.df.rds)

# Find emiRs
print("Finding emiRs...")
df.results %>%
	group_by(UniName) %>%
	summarise(min.pval = min(P)) %>%
	filter(min.pval < FDR.THRESHOLD) %>%
	pull(UniName) -> emirs

print(paste0("Found ", length(emirs), " emiRs."))

#df.emirs <- data.frame(emir = emirs, covMat = covMat, stringsAsFactors = FALSE)

#print("Saving emiRs summary data frame...")
#saveRDS(df.emirs, paste0(summary.dir, covMat, ".dataFrame.rds"))

# Find nominal p-value
# print("Finding p-value...")
# df.results %>%
# 	mutate(P.adj = p.adjust(P, method = "fdr")) %>%
# 	filter(P.adj < FDR.THRESHOLD) %>%
# 	top_n(1, P) %>%
# 	pull(P) -> p.value.nominal.threshold
#
# p.value.nominal.threshold <- p.value.nominal.threshold[1]
#
# print(p.value.nominal.threshold)
#
# print("Saving p-value...")
# write_lines(p.value.nominal.threshold, output.nom.p.value.file)

print("Finished!")
