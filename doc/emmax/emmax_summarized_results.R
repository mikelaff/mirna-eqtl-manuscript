# Summarize EMMAX results

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)
library(reshape2)
library(lattice)
library(DESeq2)
library(mikelaffr)

date.prefix <- format(Sys.time(), "%Y%m%d")

# OUTPUT FILES ########################################################################################################
dir.pdfs <- here("doc/emmax/pdfs/")
dir.pngs <- here("doc/emmax/pngs/")

# list of eGenes
egenes.tsv <- here("results/emmax/association_results/20190808_mirQTL/eGenes.tsv")
# list of eGene .ps files
egenes.ps.files <- here("results/emmax/association_results/20190808_mirQTL/eGene.ps.files")
# nominal p-value threshold
nom.pval.thresh.txt <- here("results/emmax/association_results/20190808_mirQTL/nominal.pval.threshold")
# list of eSNPs
esnps.txt <- here("results/emmax/association_results/20190808_mirQTL/eSNPs.txt")

# INPUT FILES #########################################################################################################
# GRangesList containing SNP-miR associations within 1MB of each miR (only loaded if needed)
grangeslist.rds <- here("results/rdata_files/20190808_mirQTL.GRangesList.rds")
# Dataframe containing the GRangesList data (computation is faster on dataframe than the GRangesList)
dataframe.rds <- here("doc/emmax/rdata/20190808_mirQTL.PCA10.cov.imputedVarients.dataframe.rds")

# samples file
#samples.txt <- here("src/emmax/20190725_fetalTissue_CW_smallRNAseq_RNAID_DonorID_DNAID.txt")

# mirna ranged summarized experiment
rse.rds <- here("results/rdata_files/20190505_rse_all_known_and_novel_counts.rds")

# index snps
index.snps.txt <- here("results/emmax/association_results/20190808_mirQTL/ld/chrAll.mirQTL.indexSNPs")

# GLOBALS #############################################################################################################
# minimum number of homozygous minor samples to consider associations
MIN.HOM.MINOR.SAMPLES <- 2

# threshold for FDR adjusted p-values for calling significant snp-mir associations
FDR.THRESHOLD <- 0.05

# Import SNP Associations #############################################################################################
# if the convertered dataframe is available, load it
# else import grangeslist and convert to dataframe
if (file.exists(dataframe.rds)) {
    df.snp.mir.assoc <- readRDS(dataframe.rds)
} else {
    grl <- readRDS(grangeslist.rds)

    # unlist GRangesList into GRanges
    gr.snps <- unlist(grl)
    # split name of each range into snpid and mir
    gr.snps$snpid <- sapply(strsplit(names(gr.snps), "\\."), `[`, 2)
    gr.snps$mir <- sapply(strsplit(names(gr.snps), "\\."), `[`, 1)
    # remove names from GRanges
    names(gr.snps) <- NULL

    # convert GRanges into dataframe
    df.snp.mir.assoc <- as.data.frame(gr.snps)

    # select and rename columns
    df.snp.mir.assoc %<>%
        select(chr = seqnames,
               coord.hg19 = start,
               snpid,
               mir, A1, A2.effect, A1.HOM.count, HET.count, A2.HOM.count, beta, pval)

    # save dataframe
    saveRDS(df.snp.mir.assoc, dataframe.rds)

    # clean up
    rm(grl, gr.snps)
}

# Import miRNA GRanges ################################################################################################
# rowranges of rse used for GRanges of miRNAs
rse <- readRDS(rse.rds)
gr.mirs <- rowRanges(rse)

# subset for miRNAs with SNP associations
gr.mirs <- gr.mirs[unique(df.snp.mir.assoc$mir)]

# 864 expressed miRNAs, 7 had no SNPs within 1 Mb = 857 miRs remaining

rm(rse)

# Find eGenes #########################################################################################################

# remove snp-mir associations with less than set number of homozygous minor samples
df.snp.mir.assoc %<>%
    filter(A1.HOM.count >= MIN.HOM.MINOR.SAMPLES)

df.snp.mir.assoc %>%
    filter(A1.HOM.count == 1) -> tmp2
unique(tmp2$mir)

df.snp.mir.assoc %<>% filter(A1.HOM.count != 1, HET.count >= 2)

# compute fdr corrected p-value
df.snp.mir.assoc$pval.adj <- p.adjust(df.snp.mir.assoc$pval, method = "fdr")

# summarize mirs by minimum adjusted p-value, threshold to find eGenes
df.snp.mir.assoc %>%
    group_by(mir) %>%
    summarise(min.pval.adj = min(pval.adj), min.pval = min(pval)) %>%
    filter(min.pval.adj < FDR.THRESHOLD) %>%
    pull(mir) -> eGenes

# find nominal p-value threshold
df.snp.mir.assoc %>%
    filter(pval.adj < FDR.THRESHOLD) %>%
    top_n(1, pval) %>%
    pull(pval) -> p.value.nominal.threshold

# GRanges of eGenes
gr.egenes <- gr.mirs[eGenes]
# find overlaps within the GRanges, drop self overlaps, drop redundant overlaps, subset GRanges
# for those not within the findOverlaps query hits
gr.egenes.nonoverlapping <- gr.egenes[! (1:length(gr.egenes)) %in% queryHits(findOverlaps(gr.egenes, drop.self=TRUE, drop.redundant=TRUE))]

eGenes.nonoverlapping <- names(gr.egenes.nonoverlapping)

# Export eGenes List ##################################################################################################
# filter GRanges and export dataframe as .tsv
if (file.exists(egenes.tsv)) {
    print(paste(egenes.tsv, "already exists. Not exporting data frame."))
} else {
    df <- as.data.frame(gr.egenes.nonoverlapping)
    write_tsv(df, egenes.tsv)
}

# write egene ps file names
if (file.exists(egenes.ps.files)) {
    print(paste(egenes.ps.files, "already exists. Not exporting .ps file list."))
} else {
    write_lines(paste(df$seqnames, ".mirQTL.", rownames(df), ".ps", sep = ""), egenes.ps.files)
}

# write nominal p-value threshold to file
if (file.exists(nom.pval.thresh.txt)) {
    print(paste(nom.pval.thresh.txt, "already exists. Not exporting nominal p-value."))
} else {
    write_lines(p.value.nominal.threshold, nom.pval.thresh.txt)
}

rm(df)

# Import Index SNPs / Get eSNPs #######################################################################################
chrAll.index.snps <- read_lines(index.snps.txt)

df.snp.mir.assoc %>%
    filter(mir %in% eGenes.nonoverlapping) %>%
    filter(snpid %in% chrAll.index.snps) -> eSNPs

sum(!duplicated(eSNPs$mir))
sum(eGenes.nonoverlapping %in% eSNPs$mir)

eGenes.nonoverlapping[!eGenes.nonoverlapping %in% eSNPs$mir]



# Manhattan Plot ######################################################################################################

# remove unused chr levels
df.snp.mir.assoc$chr <- droplevels(df.snp.mir.assoc$chr)
# chr number for plotting
df.snp.mir.assoc$chrLabel <- df.snp.mir.assoc$chr
levels(df.snp.mir.assoc$chrLabel) <- c(1:22,"X")

df.snp.mir.assoc %>%
    filter(mir %in% eGenes.nonoverlapping) %>%
    group_by(mir) %>%
    filter(pval == min(pval)) %>%
    filter(1:n() == 1) -> df.egenes

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")

chrlengths <- data.frame(chr = factor(names(seqlengths(gr.egenes)[1:23]), levels = chroms), coord.hg19 = seqlengths(gr.egenes)[1:23])
chrlengths$chrLabel <- chrlengths$chr
levels(chrlengths$chrLabel) <- c(1:22,"X")

chrlengths$pval <- 0.09

tmp <- chrlengths
tmp$coord.hg19 <- 1

#chrlengths$coord.hg19[23] <- 500e6

pdf("~/Desktop/eqtl_man.pdf", height = 8, width = 16, useDingbats = FALSE)
df.snp.mir.assoc %>%
    bind_rows(chrlengths, tmp) %>%
    filter(pval < 0.1) %>%
    filter(chr != "chrX") %>%
    ggplot() +
    geom_point(aes(x = coord.hg19, y = -log10(pval), color = chrLabel), alpha = 0.7, size = 0.2) +
    facet_grid(~chrLabel, scales = "free_x", space = "free_x", switch = "x") +
    geom_hline(yintercept = -log10(p.value.nominal.threshold), color = "black", linetype = "solid") +
    labs(x = "Chromosome",
         y = expression(paste("-lo",g[10],"(",italic("p-value"),")"))) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks = element_blank(),
          panel.spacing.x = unit(0, "mm"),
          panel.border = element_blank(),
          panel.grid = element_blank(),
          panel.background = element_blank(),
          strip.background = element_blank(),
          axis.title.x = element_blank(),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 14),
          strip.text = element_text(size = 12)) +
    scale_color_manual(values = c(rep(c("navy","cornflowerblue"),12))) +
    scale_y_continuous(expand=expand_scale(mult = c(0,0.05))) +
    scale_x_continuous(expand=c(0,0)) +
    geom_point(mapping = aes(x = coord.hg19, y = -log10(pval)), data = filter(df.egenes, chr != "chrX"), size=1, color="red")
    #geom_text_repel(mapping = aes(x = coord.hg19, y = -log10(pval), label = mir), data = df.egenes, size=2, angle=90)
dev.off()

ggsave("~/Desktop/eqtl_man_rast.png", height = 4, width = 10, units = "in", dpi = "retina")

ggsave("~/Desktop/eqtl_man.pdf", height = 8, width = 16, useDingbats = FALSE)

# Stuff #######
stop()
overlappping <- findOverlaps(gr, drop.self=TRUE, drop.redundant=TRUE)

gr[subjectHits(overlappping)]
grQ <- gr[queryHits(overlappping)]

gr2 <- gr[! (1:length(gr)) %in% subjectHits(findOverlaps(gr, drop.self=TRUE, drop.redundant=TRUE))]



findOverlaps(gr2, drop.self=TRUE, drop.redundant=TRUE)

%>%
    filter(mir %in% eGenes) -> tmp

summary(tmp$pval.adj)









# plot p-value dist?
#1. compute LambdaGC for each miR -> plot dist
#2. compute lambdaGC for all p-values -> plot qq
#3. get

df.snp.mir.assoc %>%
    ggplot(aes(x=pval)) +
    geom_histogram(bins = 20)

df.snp.mir.assoc %>%
    sample_n(1e4) %>%
    arrange(pval.adj) %>%
    filter(pval.adj < 0.05) -> tmp



    ggplot(aes(x=pval.adj, y=pval)) +
    geom_point()








# a function to compute LambdaGC (Genomic Control)
genomicControl <- function(pvalues) {
    chisq <- qchisq(pvalues, 1, lower.tail = FALSE)
    median(chisq) / qchisq(0.5, 1)
}


df.snps %>%
    filter(mir %in% tmp$mir) -> df.snps.emirs

df.snps.emirs %>%
    group_by(mir) %>%
    filter(pval.adj == min(pval.adj)) -> tmp2

df.snp.mir.assoc %>%
    filter(mir == "hsa-miR-9-3p") %>%
    ggplot(aes(x = coord.hg19, y = -log10(pval))) +
    geom_point()




df.snps %>%
    ggplot(aes(x = pval.adj)) +
    geom_density()



gr.mirs$lambdaGC <- NA
gr.mirs$min.adj.pval <- NA

# total number of associations
numAssociations <- sum(sapply(grl, length))

#pdf("~/Desktop/qqplots.pdf")
for (i in 1:length(grl)) {
#for (i in 1:3) {

    mir <- names(grl)[i]

    pvalues <- grl[[i]]$pval
    expected <- (rank(pvalues, ties.method="first")+.5)/(length(pvalues)+1)

    lambdaGC <- genomicControl(pvalues)

    # evalues <- ppoints(pvalues)

    # fdr adjusted pvalues
    grl[[i]]$pval.fdr <- p.adjust(grl[[i]]$pval, method = "fdr", n = numAssociations)

    gr.mirs[mir]$lambdaGC <- lambdaGC
    gr.mirs[mir]$min.adj.pval <- min(grl[[i]]$pval.fdr)

    print(paste(mir, ": lambdaGC =", lambdaGC, "min adj pval:", gr.mirs[mir]$min.adj.pval))

    # plot(-log10(expected), -log10(pvalues))

    df <- data.frame(observed = pvalues, expected)

    df %>%
        ggplot(aes(x=-log10(expected), y=-log10(observed))) +
        geom_point(size=.5, color="blue") +
        geom_abline(intercept = 0, slope = 1) +
        labs(title = mir) +
        annotate("text", x = 0.5, y = 2, label = as.expression(bquote(lambda[GC] == .(round(lambdaGC, 2)))), size=5) +
        plotTheme() -> p

    plot(p)

}

#dev.off()

stop()

df.grl <- as.data.frame(grl)

gr <- unlist(grl)

pvals.all <- df.snps$pval

for (i in 1:length(grl)) {
    print(i)
    pvals.all <- c(pvals.all, grl[[i]]$pval)
}

pvals.all <- df.snp.mir.assoc$pval
expected <- (rank(pvals.all, ties.method="first")+.5)/(length(pvals.all)+1)
df <- data.frame(observed = pvals.all, expected)

df.samp <- sample_n(df, 1e5)

lambdaGC <- genomicControl(pvals.all)

df.samp %>%
    ggplot(aes(x=-log10(expected), y=-log10(observed))) +
    geom_point(size=.5, color="blue") +
    geom_abline(intercept = 0, slope = 1) +
    labs(title = "All p-values") +
    plotTheme(theme = 'presentation') +
    annotate("text", x = 1, y = 8, label = as.expression(bquote(lambda[GC] == .(round(lambdaGC, 2)))), size=10)

ggsave(paste(dir.pngs, "qq_all_pvals_pres.png", sep="/"), height=6, width=6)

df.snp.mir.assoc$chrLabel <- as.factor(sapply(strsplit(as.character(df.snp.mir.assoc$chr), "chr"), `[`, 2))

df.snp.mir.assoc %>%
    filter(pval < 0.1) %>%
    mutate(ord = letters[chr]) %>%
    ggplot() +
    geom_point(aes(x = coord.hg19, y = -log10(pval), color = ord), alpha = 0.7, size = 0.2) +
    facet_grid(~chrLabel, scales = "free_x", space = "free_x", switch = "x") +
    geom_hline(yintercept = -log10(p.value.nominal.threshold), color = "grey", linetype = "dotted") +
    labs(x = "Chromosome",
         y = expression(paste("-lo",g[10],"(",italic("P"),")"))) +
    theme_bw() +
    theme(axis.ticks.y = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          panel.spacing.x = unit(0.0,"cm"),
          strip.background = element_blank(),
          panel.border = element_blank(),
          axis.text.x  = element_blank(),
          legend.position = "na",
          axis.title = element_text(size=25),
          axis.text.y = element_text(size = 14),
          strip.text.x = element_text(size=14)) +
    scale_color_manual(values = c(rep(c("navy","cornflowerblue"),12))) +
    scale_y_continuous(expand=c(0,0.5))+scale_x_continuous(expand=c(0,0))

as.data.frame(mcols(gr.mirs)) %>%
    ggplot(aes(lambdaGC)) +
    geom_histogram(color="white", bins = 25) +
    labs(title = "Lambda GC for Associations within 1Mb of each miRNA") +
    scale_y_continuous(expand = expand_scale(mult = c(0,.1))) +
    scale_x_continuous(limits = c(0,3)) +
    plotTheme() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

ggsave(paste(dir.pdfs, "lambdaGC_dist.pdf", sep="/"), height=8, width=8)


gr.mir.sub <- gr.mirs[mcols(gr.mirs)$min.adj.pval < 0.05]

as.data.frame(mcols(gr.mir.sub)) %>%
    ggplot(aes(lambdaGC)) +
    geom_histogram(color="white", bins = 25) +
    labs(title = "Lambda GC for Associations within 1Mb of each miRNA") +
    scale_y_continuous(expand = expand_scale(mult = c(0,.1))) +
    scale_x_continuous(limits = c(0,3)) +
    plotTheme() +
    theme(panel.border = element_rect(color = "black", size = 1, fill = NA))

# Import Data #########################################################################################################

chroms <- c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
            "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX")




# Scratch ######################

stop()


genomicControl <- function(pvalues) {
    chisq <- qchisq(1 - pvalues, 1)
    median(chisq) / qchisq(0.5, 1)
}


mean(sapply(grl, length))

pvalues <- grl[[1]]$pval

hist(pvalues, breaks = 50)

chisq <- qchisq(1-pvalues,1)

genomicControl(pvalues)

lambda <- median(chisq)/qchisq(0.5,1)

qqmath(~-log10(pvalues),
       distribution=function(x){-log10(qunif(1-x))}, )


tmp <- dplyr::select(df, id, `hsa-miR-1306-3p.pval`, `hsa-miR-130b-3p.pval`)

melted <- melt(tmp, id.vars = "id")

observed <- -log10(melted$value)
expected <- -log10(ppoints(length(melted$value)))

plot(order(expected), order(observed))

hist((rank(pvalues, ties.method="first")+.5)/(length(pvalues)+1))
hist(ppoints(pvalues))



df.mirna <-

df <- read_table("~/Desktop/chr22.FetalTissueCWSmallRNASeqSamples.chr22_41896_mat.clumped")

df <- read.table("~/Desktop/chr22.FetalTissueCWSmallRNASeqSamples.chr22_41896_mat.clumped", header = TRUE, sep = "")

# get geno pheno for example mirnas
mir <- "hsa-miR-1304-3p"

df.snp.mir.assoc %>%
    dplyr::mutate(log10p = -log10(pval)) %>%
    top_n(1, log10p)


dplyr::filter(mir == mir) %>%
    dplyr::mutate(log10p = log10(pval)) %>%
    top_n(1, log10p) %>%
    pull(log10p)

df.snp.mir.assoc %>%
    group_by(mir) %>%
    summarise(min.pval.adj = min(pval.adj), min.pval = min(pval), snpid = ) %>%
    filter(min.pval.adj < FDR.THRESHOLD) -> tmp
