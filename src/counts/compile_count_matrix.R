# compile mirna count matrices from featureCounts against small-rna-seq bams
# all_known_and_novel
# mirdeep_mirge
# mirdeep_mirge_friedlander_nowakowski

library(here)
library(readr)
library(dplyr)
library(magrittr)

# INPUT FILES #########################################################################################################
# counts directory roots
root.dir.all <- here("results/counts/small_rna_seq/all_known_and_novel/")
root.dir.mirdeep.mirge <- here("results/counts/small_rna_seq/mirdeep_mirge/")
root.dir.mirdeep.mirge.fried.nowa <- here("results/counts/small_rna_seq/mirdeep_mirge_friedlander_nowakowski/")

# OUTPUT FILES ########################################################################################################
output.tsv.all <- paste(here("results/counts/small_rna_seq/"), format(Sys.time(), "%Y%m%d"),
                        "_all_known_and_novel_counts.tsv", sep="")
output.tsv.mirdeep.mirge <- paste(here("results/counts/small_rna_seq/"), format(Sys.time(), "%Y%m%d"),
                        "_mirdeep_mirge_novel_counts.tsv", sep="")
output.tsv.mirdeep.mirge.fried.nowa <- paste(here("results/counts/small_rna_seq/"), format(Sys.time(), "%Y%m%d"),
                                  "_mirdeep_mirge_friedlander_nowakowski_counts.tsv", sep="")

# Compile all_known_and_novel Results #################################################################################
# files in directory
files <- list.files(root.dir.all)

# load first count file
countData <- read_tsv(paste(root.dir.all, files[1], sep="/"),
                      skip = 2,
                      col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                    "derives_from", "seq_type", strsplit(files[1], "[.]")[[1]][1]))
# remove duplicated row (mistake in nowakowski gtf file)
countData <- countData[!duplicated(countData$ID),]

#loop through all files in count folder and add count data to countData object
for (i in 2:length(files)) {
  #load count file
  tempData <- read_tsv(paste(root.dir.all, files[i], sep="/"),
                       skip = 2,
                       col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                     "derives_from", "seq_type", strsplit(files[i], "[.]")[[1]][1]))
  tempData <- tempData[!duplicated(tempData$ID),]
  # select only ID and count column
  tempData %<>% dplyr::select(7,11)
  #join tempData with countData
  countData <- left_join(countData, tempData, by="ID")
}

countData.all <- countData

#remove objects
rm(countData, tempData, i, files)

# Compile mirdeep.mirge Results #######################################################################################
# files in directory
files <- list.files(root.dir.mirdeep.mirge)

# load first count file
countData <- read_tsv(paste(root.dir.mirdeep.mirge, files[1], sep="/"),
                      skip = 2,
                      col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                    "derives_from", "seq_type", strsplit(files[1], "[.]")[[1]][1]))
# remove duplicated row (mistake in nowakowski gtf file)
countData <- countData[!duplicated(countData$ID),]

#loop through all files in count folder and add count data to countData object
for (i in 2:length(files)) {
  #load count file
  tempData <- read_tsv(paste(root.dir.mirdeep.mirge, files[i], sep="/"),
                       skip = 2,
                       col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                     "derives_from", "seq_type", strsplit(files[i], "[.]")[[1]][1]))
  tempData <- tempData[!duplicated(tempData$ID),]
  # select only ID and count column
  tempData %<>% dplyr::select(7,11)
  #join tempData with countData
  countData <- left_join(countData, tempData, by="ID")
}

countData.mirdeep.mirge <- countData

#remove objects
rm(countData, tempData, i, files)

# Compile mirdeep.mirge.fried.nowa Results ############################################################################
# files in directory
files <- list.files(root.dir.mirdeep.mirge.fried.nowa)

# load first count file
countData <- read_tsv(paste(root.dir.mirdeep.mirge.fried.nowa, files[1], sep="/"),
                      skip = 2,
                      col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                    "derives_from", "seq_type", strsplit(files[1], "[.]")[[1]][1]))
# remove duplicated row (mistake in nowakowski gtf file)
countData <- countData[!duplicated(countData$ID),]

#loop through all files in count folder and add count data to countData object
for (i in 2:length(files)) {
  #load count file
  tempData <- read_tsv(paste(root.dir.mirdeep.mirge.fried.nowa, files[i], sep="/"),
                       skip = 2,
                       col_names = c("name", "chr", "start", "end", "strand", "length", "ID", "alias",
                                     "derives_from", "seq_type", strsplit(files[i], "[.]")[[1]][1]))
  tempData <- tempData[!duplicated(tempData$ID),]
  # select only ID and count column
  tempData %<>% dplyr::select(7,11)
  #join tempData with countData
  countData <- left_join(countData, tempData, by="ID")
}

countData.mirdeep.mirge.fried.nowa <- countData

#remove objects
rm(countData, tempData, i, files)

# Export Count Matrices ###############################################################################################

# convert to integer
countData.all %<>% mutate_if(is.double, as.integer)
countData.mirdeep.mirge %<>% mutate_if(is.double, as.integer)
countData.mirdeep.mirge.fried.nowa %<>% mutate_if(is.double, as.integer)

# write counts data table
write_tsv(countData.all, output.tsv.all, col_names = TRUE)
write_tsv(countData.mirdeep.mirge, output.tsv.mirdeep.mirge, col_names = TRUE)
write_tsv(countData.mirdeep.mirge.fried.nowa, output.tsv.mirdeep.mirge.fried.nowa, col_names = TRUE)

# Scratch #############################################################################################################
# are the counts from mirdeep.mirge the same as mirdeep.mirge.fried.nowa?
tmp <- filter(countData.mirdeep.mirge.fried.nowa, ID %in% countData.mirdeep.mirge$ID)
all.equal(tmp, countData.mirdeep.mirge)
rm(tmp)

# are the counts of mirbase mirnas in all similar to the mirge2.0 counts?
cts.mirge <- read_tsv(here("results/mirge2.0/20180905_mirbase_v22_mirge2.0_counts.tsv"))
cts.mirge$miRNA <- sapply(strsplit(cts.mirge$miRNA, "/"), `[`, 1)
sum(cts.mirge$miRNA %in% countData.all$name)

tmp <- cts.mirge[!cts.mirge$miRNA %in% countData.all$name,]
mat <- as.matrix(tmp[,2:241])
rownames(mat) <- tmp$miRNA

mat2 <- mat[rowSums(mat >= 10) >= 10,]

mirge <- as.matrix(cts.mirge[,2:241])
rownames(mirge) <- cts.mirge$miRNA

allcts <- as.matrix(countData.all[,11:250])
rownames(allcts) <- countData.all$name

mirge <- mirge[rownames(mirge) %in% rownames(allcts),]
allcts <- allcts[rownames(allcts) %in% rownames(mirge),]

df <- data.frame(names = rownames(allcts),
                 all_rowsums = rowSums(allcts),
                 all_rowmeans = rowMeans(allcts),
                 all_rowvars = rowVars(allcts),
                 all_rowmin = rowMins(allcts),
                 all_rowmax = rowMaxs(allcts))

df2 <- data.frame(names = rownames(mirge),
                  mirge_rowsums = rowSums(mirge),
                  mirge_rowmeans = rowMeans(mirge),
                  mirge_rowvars = rowVars(mirge),
                  mirge_rowmin = rowMins(mirge),
                  mirge_rowmax = rowMaxs(mirge))

df3 <- left_join(df, df2, by = "names")
df3$dup_names <- duplicated(df3$names) | duplicated(df3$names, fromLast = TRUE)

library(ggplot2)

df3 %>%
  filter(!dup_names) %>%
  ggplot(aes(log10(all_rowmeans+1), log10(mirge_rowmeans+1))) +
  geom_point() +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0)

df3 %>%
  ggplot(aes(log10(all_rowsums+1), log10(mirge_rowsums+1), color=dup_names)) +
  geom_point()

df3 %>%
  ggplot(aes(all_rowmeans, mirge_rowmeans)) +
  geom_point() +
  geom_smooth(method = "lm")
