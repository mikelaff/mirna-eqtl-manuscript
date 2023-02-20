# GOexpress tutorial

library(GOexpress)
data(AlvMac)

exprs(AlvMac)[1:5,1:5]

head(pData(AlvMac))

is.factor(AlvMac$Treatment)

AlvMac$Treatment

set.seed(4543)
AlvMac_results <- GO_analyse(eSet = AlvMac,
                             f = "Treatment",
                             GO_genes = AlvMac_GOgenes,
                             all_GO = AlvMac_allGO,
                             all_genes = AlvMac_allgenes)
