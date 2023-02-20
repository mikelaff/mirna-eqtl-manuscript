# create list of snps which pass the homozygous minor sample/het sample filters.
# this list will be used to filter hardcall and dosage files

library(here)
library(readr)
library(dplyr)
library(magrittr)

# input data directory
dir.hardcall <- here("results/genotypes/mirQTL_imputed_genotypes/plink.hardcall/")

# dosage data directory
dir.dosage <- here("results/genotypes/mirQTL_imputed_genotypes/plink.dosage/")

CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# loop over each chrom
for (chr in CHROMS) {

	print(chr)

	snps <- NULL
	df.frqx <- NULL

	# read frqx file
	df.frqx <- read_tsv(paste0(dir.hardcall, chr, ".topmed.hardcall.r2g03.maf001.mirQTL.frqx"))
	#df.frqx <- read_tsv("~/Downloads/chr22.topmed.hardcall.r2g03.maf001.mirQTL.frqx")

	# filter for HOM.A1 != 1 and HET > 1
	df.frqx %>%
		dplyr::filter(`C(HOM A1)` != 1 & `C(HET)` > 1) %>%
		pull(SNP) -> snps

	write_lines(snps, paste0(dir.hardcall, chr, ".topmed.hardcall.r2g03.maf001.mirQTL.HomHetStrict.snps"))
	write_lines(snps, paste0(dir.dosage, chr, ".topmed.dosage.r2g03.maf001.mirQTL.HomHetStrict.snps"))

}


