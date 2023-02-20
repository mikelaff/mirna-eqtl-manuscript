# try out miRNA target prediction

library(miRNAtap)
library(topGO)
library(org.Hs.eg.db)
library(readr)

# tutorial from miRNAtap manual
mir <- "miR-128-3p"
predictions <- getPredictedTargets(mir, species = "hsa", method = "geom", min_src = 2)


rankedGenes <- predictions[,'rank_product']
selection <- function(x) TRUE
allGO2genes <- annFUN.org(whichOnto = "BP",
                          feasibleGenes = NULL,
                          mapping = "org.Hs.eg.db",
                          ID = "entrez")

GOdata <- new("topGOdata",
              ontology = "BP",
              allGenes = rankedGenes,
              annot = annFUN.GO2genes,
              GO2genes = allGO2genes,
              geneSel = selection,
              nodeSize = 10)

results.ks <- runTest(GOdata, algorithm = "classic", statistic = "ks")

allRes <- GenTable(GOdata, KS = results.ks, orderBy = "KS", topNodes = 20)


write_lines(dimnames(predictions)[[1]], "~/Desktop/miR-128-3p_targets.txt")



library(readxl)

targets <- read_xlsx("~/Downloads/hsa_MTI.xlsx")

neuro_up <- filter(targets, `miRNA` %in% miRNAs_neurogenesis_up)
neuro_down <- filter(targets, `miRNA` %in% miRNAs_neurogenesis_down)

matur_up <- filter(targets, `miRNA` %in% miRNAs_maturation_up)
matur_down <- filter(targets, `miRNA` %in% miRNAs_maturation_down)

neuro_up <- neuro_up[!duplicated(neuro_up$`Target Gene`),]
neuro_down <- neuro_down[!duplicated(neuro_down$`Target Gene`),]



# 44,310 genes targeted by one or more mirnas in maturation list
# 55,270 genes targeted by one or more mirnas in neurogenesis list

sum(!duplicated(tmp_neuro$`Target Gene`))

sum(!duplicated(tmp_mat$`Target Gene`))

# CTIP2 (BCL11B) ##################
tmp <- filter(targets, `Target Gene` == "BCL11B")

ctip2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# TBR1 ##################
tmp <- filter(targets, `Target Gene` == "TBR1")

tbr1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# FOXP2 ##################
tmp <- filter(targets, `Target Gene` == "FOXP2")

foxp2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# CUX1 ##################
tmp <- filter(targets, `Target Gene` == "CUX1")

cux1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# CUX2 ##################
tmp <- filter(targets, `Target Gene` == "CUX2")

cux2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# SATB2 ##################
tmp <- filter(targets, `Target Gene` == "SATB2")

satb2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# BRN1 ##################
tmp <- filter(targets, `Target Gene` == "POU3F3")

brn1_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# BRN2 ##################
tmp <- filter(targets, `Target Gene` == "POU3F2")

brn2_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]

# RELN ##################
tmp <- filter(targets, `Target Gene` == "RELN")

reln_mirnas <- tmp$miRNA[!duplicated(tmp$miRNA)]
