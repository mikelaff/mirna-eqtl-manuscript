# test microarray processing


library(here)
library(readr)
library(dplyr)
library(tidyr)
library(magrittr)
#library(ggplot2)

library(ArrayExpress)
library(Biobase)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
#library(AffyCompatible)
library(oligo)
library(arrayQualityMetrics)

library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)

library(gplots)
library(ggplot2)
library(geneplotter)
library(RColorBrewer)
library(pheatmap)
library(enrichplot)

library(stringr)
library(matrixStats)
library(genefilter)
library(openxlsx)

library(mikelaffr)

# OUTPUT ###############################################################################################################

# INPUT ################################################################################################################

# GLOBALS ##############################################################################################################


# Import Samples #######################################################################################################

raw_dat_dir <- "/Users/mikelaff/Desktop/temp/"

#anno_AE <- getAE("E-MTAB-2967", path = raw_dat_dir, type = "raw")


sdrf <- read.delim(file.path(raw_dat_dir, "E-MTAB-2967.sdrf.txt"))

rownames(sdrf) <- sdrf$Array.Data.File
sdrf <- AnnotatedDataFrame(sdrf)

raw_data <- read.celfiles(filenames = file.path(raw_dat_dir, sdrf$Array.Data.File), phenoData = sdrf)

file.exists(raw_dat_dir, sdrf$Array.Data.File)

list.files(raw_dat_dir)

head(pData(raw_data))

raw_data

rownames(exprs(raw_data))

rownames(raw_data)

featureNames(raw_data)

rma_data <- oligo::rma(raw_data, target = "core")

featureNames(rma_data)

rownames(rma_data)

anno_palmieri <- AnnotationDbi::select(hugene10sttranscriptcluster.db,
                                       keys = (featureNames(rma_data)),
                                       columns = c("SYMBOL", "GENENAME"),
                                       keytype = "PROBEID")


anno_palmieri <- subset(anno_palmieri, !is.na(SYMBOL))

