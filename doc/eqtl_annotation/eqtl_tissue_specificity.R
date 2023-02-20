# look eqtls and blood eqtls

library(here)
library(dplyr)
library(readr)
library(readxl)
library(magrittr)
library(ggplot2)
library(DESeq2)
library(AnnotationHub)
library(GenomicRanges)
library(rtracklayer)
library(liftOver)
library(mikelaffr)

# OUTPUT ###############################################################################################################
dir.pdfs <- here("doc/eqtl_annotation/pdfs/")
dir.create(dir.pdfs, showWarnings = FALSE, recursive = TRUE)

# summarized results of miRNAs and host expression
#df.output.rds <- here("results/rdata_files/2020")

# INPUT ################################################################################################################
# mirQTL eQTLs
eqtls.dataframe.rds <- here("results/conditional_eqtls/20200120_mirQTLor_compiled/20200120_mirQTLor_conditional_eQTLs_compiled_dataFrame.rds")
# blood miRNA-eQTLs
huan.eqtls.xlsx <- here("data/huan_2015/41467_2015_BFncomms7601_MOESM1319_ESM.xlsx")

# target predictions by miRNAtap databases
predictions.2sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_2sources.csv")
predictions.3sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_3sources.csv")
predictions.4sources.csv <- here("data/target_predictions/20190514_miRNAtap_predictions_4sources.csv")

# miRDB 2019 predictions, known mirnas
predictions.mirdb.txt.gz <- here("data/target_predictions/miRDB_v6.0_prediction_result.txt.gz")

# blood overlaps
blood.colocs.rds <- here("results/co-localization/blood_miRNA-eQTL/blood_miRNA-eQTL_mirQTL_overlaps_r2at0.8.rds")

# blood emir sets
blood.emir.set.rds <- here("results/co-localization/blood_miRNA-eQTL/blood_miRNA-eQTLmirQTL_blood-miRNA-eQTL_emiR_sets.rds")






# hg19 to hg38 chain file
chain.file <- here("data/hg19ToHg38.over.chain")

# GRanges of known and novel miRNAs
gr.mirna.rds <- here("data/gtf_and_granges/20200227_all_known_and_novel_mirna_non_overlapping_granges.rds")

# GLOBALS ##############################################################################################################
CHROMS <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

# chromosome lengths
#CHROM.LENGTHS <- seqlengths(Seqinfo(genome = "hg38"))[1:23]

# Annotation Hub ######################################################################################################
ah <- AnnotationHub()
# homo sapiens org db
orgs <- subset(ah, ah$rdataclass == "OrgDb")
orgdb <- query(orgs, "Homo sapiens")[[1]]
rm(ah, orgs)

# Import eQTL Data #####################################################################################################
df.mirQTL <- readRDS(eqtls.dataframe.rds)

df.mirQTL %>%
    filter(SIGNIFICANCE == "eigenMT_fdr5percent") -> df.mirQTL.high

df.blood <- read_xlsx(huan.eqtls.xlsx, skip = 1)

# chain for hg19 to hg38 conversion
ch <- import.chain(chain.file)

# GRanges of huan eqtls
gr <- makeGRangesFromDataFrame(df.blood,
                               seqnames.field = "chr.SNP",
                               start.field = "SNP.pos",
                               end.field = "SNP.pos",
                               strand.field = "SNP.strand",
                               ignore.strand = TRUE,
                               keep.extra.columns = TRUE)
# GRangesList of GRanges conversion
lo <- liftOver(gr, ch)

gr.huan.hg38 <- unlist(lo)

rm(lo, gr, ch)

df.blood <- as_tibble(as.data.frame(gr.huan.hg38))

df.emirs <- readRDS(blood.emir.set.rds)

# all known mirnas
df.mirna <- as_tibble(read_rds(gr.mirna.rds))

df.emirs %<>%
    left_join(df.mirna, by = c("emir" = "Name"))

df.emirs %>%
    filter(set == "blood_only") %>%
    pull(emir) %>%
    unique() %>%
    write_lines("~/Desktop/blood_mirnas.txt")

# Import Target Predictions ############################################################################################
# miRDB predictions
df.mirdb.predictions <- read_tsv(predictions.mirdb.txt.gz, col_names = c("Name", "REFSEQ", "TargetScore", "ENTREZID"))
df.mirdb.predictions %<>%
    filter(Name %in% df.emirs$emir)

# get ensembl and entrezids
rowdata <- AnnotationDbi::select(orgdb,
                                 keys = unique(df.mirdb.predictions$REFSEQ),
                                 keytype = "REFSEQ",
                                 columns = c("ENSEMBL", "SYMBOL", "GENENAME"))

# combine with predictions df
df.mirdb.predictions %<>%
    dplyr::left_join(rowdata, by = "REFSEQ")

df.mirdb.predictions %<>%
    dplyr::filter(!is.na(ENSEMBL))

df.mirdb.predictions %>%
    filter(Name %in% df.emirs$emir[df.emirs$set == "blood_only"]) %>%
    filter(TargetScore >= 90) -> df.targets.blood.only

df.mirdb.predictions %>%
    filter(Name %in% df.emirs$emir[df.emirs$set == "brain_only"]) %>%
    filter(TargetScore >= 90) -> df.targets.brain.only

df.mirdb.predictions %>%
    filter(Name %in% df.emirs$emir[df.emirs$set == "brain_and_blood"]) %>%
    filter(TargetScore >= 90) -> df.targets.both.only



write_lines(unique(df.targets.blood.only$SYMBOL), "~/Desktop/blood_only.txt")

write_lines(unique(df.targets.brain.only$SYMBOL), "~/Desktop/brain_only.txt")

write_lines(unique(df.targets.both.only$SYMBOL), "~/Desktop/both_only.txt")

# GO Results ###############

df.go.blood <- read_tsv(here("results/gene_ontology/blood_brain_targets/blood_only_analysis.txt"), skip = 6)

df.go.brain <- read_tsv(here("results/gene_ontology/blood_brain_targets/brain_only_analysis.txt"), skip = 6)

df.go.both <- read_tsv(here("results/gene_ontology/blood_brain_targets/both_only_analysis.txt"), skip = 6)



df.go.blood %>%
    #top_n(10, -log10(`upload_1 (FDR)`)) %>%
    slice_head(n = 10) %>%
    ggplot(aes(y = -log10(`upload_1 (FDR)`), x = reorder(`GO biological process complete`, log10(`upload_1 (FDR)`)))) +
    geom_col() +
    theme(axis.text.x = element_text(angle = 90))




