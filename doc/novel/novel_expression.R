# look for novel miRNAs

library(here)
library(dplyr)
library(readr)
library(magrittr)
#library(Gviz)
#library(GenomicRanges)
#library(GenomicAlignments)
library(DESeq2)
#library(biomaRt)
library(rtracklayer)
library(Biostrings)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)

# OUTPUT ################
output_novel_rds <- paste(here("results/rdata_files/"), format(Sys.time(), "%Y%m%d"), "_novel_counts.rds", sep="")


# INPUT ################
# novel count matrix
novel_counts_tsv <- here("results/known_and_novel_mirna/20190421_novel_counts.tsv")
# summarized experiment with mirge2.0 quantified mirna counts from mirbase v22
mirbase_se_rds <- here("results/rdata_files/20180906_mirbase_v22_mirge2.0_counts.rds")
# granges for novel mirnas
mirdeep_granges_rds <- here("results/rdata_files/20190421_mirdeep2_granges.rds")
mirge_granges_rds <- here("results/rdata_files/20190421_mirge_granges.rds")
mirbase_granges_rds <- here("results/rdata_files/20190421_mirbase_granges.rds")

# Import GRanges Objects ###########
gr.mirdeep <- readRDS(mirdeep_granges_rds)
gr.mirge <- readRDS(mirge_granges_rds)
gr.mirbase <- readRDS(mirbase_granges_rds)

#df.mirdeep <- as.data.frame(gr.mirdeep)
#df.mirge <- as.data.frame(gr.mirge)

#df.novel <- bind_rows(df.mirdeep, df.mirge)

gr.novel <- c(gr.mirdeep, gr.mirge)
names(gr.novel) <- gr.novel$Name

# Import Sample Data ##########
se <- readRDS(mirbase_se_rds)
samples <- as.data.frame(colData(se))
rm(se)

# Import Count Data ############
novelCounts <- read_tsv(novel_counts_tsv)
cts <- as.matrix(novelCounts[,2:241])
rownames(cts) <- novelCounts$miRNA

# Create RSE ######################################################################################
# order rowranges on rownames of cts
rowranges <- gr.novel[rownames(cts)]
stopifnot(all(names(rowranges) == rownames(cts)))

# order count matrix
cts <- cts[, rownames(samples)]
stopifnot(all(rownames(samples) == colnames(cts)))

# build RangedSummarizedExperiment
rse <- SummarizedExperiment(assays = list(counts = cts),
                            rowRanges = rowranges,
                            colData = samples)

saveRDS(rse, output_novel_rds)

# Threshold ##########
# build DESeq data set
dds <- DESeqDataSet(rse, design = ~1)
# keep only rows from mature miRNAs
dds <- dds[rowData(dds)$type == "miRNA_putative_mature" | rowData(dds)$type == "miRNA_putative_star", ]
# keep only rows with at least 10 counts in at least 10 samples
dds <- dds[rowSums(counts(dds) >= 10) >= 10, ]

# non overlapping ranges
hits <- findOverlaps(dds, ignore.strand=TRUE,
                     drop.self=TRUE, drop.redundant=TRUE)
dds <- dds[! 1:158 %in% subjectHits(hits),]

gr.thresh <- rowRanges(dds)

hits <- findOverlaps(gr.thresh, gr.mirbase, ignore.strand=TRUE)

dds <- dds[! 1:123 %in% queryHits(hits),]

# Scratch #######
countData <- read_tsv(novel_counts_tsv)
countMat <- as.matrix(countData[,2:241])
rownames(countMat) <- countData$miRNA

# at least 10 counts in at least 10 samples

countMat.trunc <- countMat[rowSums( countMat >= 10 ) >= 10,]

rse.novel <- SummarizedExperiment()

pdfs <- list.files(here("doc/coverage_plots/pdfs/novel/sort/all/"))



df <- data.frame(files = pdfs, stringsAsFactors = FALSE)
df$trunc <- sapply(strsplit(df$files, "\\."), `[`, 1)
df$one <- sapply(strsplit(df$trunc, "_"), `[`, 1)
df$two <- sapply(strsplit(df$trunc, "_"), `[`, 2)
df$three <- sapply(strsplit(df$trunc, "_"), `[`, 3)

df$Derives_from <- NA
df$Derives_from[1:564] <- paste(df$two[1:564], df$three[1:564], sep="_")
df$Derives_from[565:620] <- paste("miRge", df$two[565:620], sep="_")

tmp <- as.data.frame(rowData(dds))

df <- dplyr::filter(df, Derives_from %in% tmp$Derives_from)

file.copy(from=paste(here("doc/coverage_plots/pdfs/novel/sort/all/"), df$files, sep=""), to=here("doc/coverage_plots/pdfs/novel/sort/threshold2/"))

saveRDS(dds, here("results/rdata_files/tmp_novel.rds"))

dds <- readRDS(here("results/rdata_files/tmp_novel.rds"))



