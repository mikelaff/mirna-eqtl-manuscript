# gene set enrichment analysis with topGO

library(topGO)
library(ALL)

#source("https://bioconductor.org/biocLite.R")
#biocLite("ALL")

data(ALL)
data(geneList)

affyLib <- paste(annotation(ALL), "db", sep=".")
library(package = affyLib, character.only = TRUE)

sampleGOdata <- new("topGOdata",
                    description = "Simple session",
                    ontology = "BP",
                    allGenes = geneList,
                    geneSel = topDiffGenes,
                    nodeSize = 10,
                    annot = annFUN.db,
                    affyLib = affyLib)


resultFisher <- runTest(sampleGOdata,
                        algorithm = "classic",
                        statistic = "fisher")

resultKS <- runTest(sampleGOdata,
                    algorithm = "classic",
                    statistic = "ks")

# elim more conservative than classic
resultKS.elim <- runTest(sampleGOdata,
                         algorithm = "elim",
                         statistic = "ks")


allRes <- GenTable(sampleGOdata,
                   classicFisher = resultFisher,
                   classicKS = resultKS,
                   elimKS = resultKS.elim,
                   orderBy = "elimKS",
                   ranksOf = "classicFisher",
                   topNodes = 10)

pValue.classic <- score(resultKS)
pValue.elim <- score(resultKS.elim)[names(pValue.classic)]

gstat <- termStat(sampleGOdata, names(pValue.classic))
gSize <- gstat$Annotated / max(gstat$Annotated) * 4
gCol <- colMap(gstat$Significant)

plot(pValue.classic, pValue.elim, xlab = "p-value classic", ylab = "p-value elim", pch = 19, cex = gSize)

showSigOfNodes(sampleGOdata, score(resultKS.elim), firstSigNodes = 5, useInfo = 'all')

