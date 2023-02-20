# build ExpressionSet from microarray raw data
# June 2022 HNP differentiation experiment


library(here)
library(readr)
library(dplyr)
library(magrittr)
#library(ggplot2)

library(Biobase)
#library(AffyCompatible)
library(oligo)
library(pd.clariom.s.human.ht)

library(mikelaffr)

# OUTPUT ###############################################################################################################
output.se.rds <- paste0(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_es_HNP_Differentiation_Microarray.rds")

# INPUT ################################################################################################################
# sample metadata
metadata.csv <- here("results/microarray/core_rna_samples.csv")

# directory for raw expression files
raw.data.dir <- here("results/microarray/Stein_Mike_Clariom_S_Human_24_08172022/")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# Import Samples #######################################################################################################
df.samples <- read_csv(metadata.csv)

df.samples %<>%
    mutate(CEL = paste0(Core_Sample_ID, ".CEL"))

df.samples$Timepoint <- factor(df.samples$Timepoint, levels = c("Week 1", "Week 2"), labels = c("Week1", "Week2"), ordered = TRUE)

# create AnnotatedDataFrame for use in ExpressionFeatureSet
adf <- AnnotatedDataFrame(data = as.data.frame(df.samples))
sampleNames(adf) <- adf$CEL


# Import Expression Matrix #############################################################################################
# list the raw .CEL files
#cel.files <- list.files(raw.data.dir, pattern = "\\.CEL$")

# import raw data and create ExpressionFeatureSet
raw_data <- read.celfiles(filenames = file.path(raw.data.dir, adf$CEL), phenoData = adf)

stopifnot(validObject(raw_data))

raw_data
# norm_data <- rma(raw_data)
#
# rownames(norm_data)
#
# #head(pData(raw_data))
#
# library(affxparser)
# readCel(filename = file.path(raw.data.dir, adf$CEL[1]))
#
# celfiles <- ReadAffy(celfile.path = raw.data.dir, filenames = adf$CEL)

# Export ExpressionFeatureSet ##########################################################################################
saveRDS(raw_data, output.se.rds)
