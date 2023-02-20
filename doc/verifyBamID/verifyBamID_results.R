# plot verify bam id results

library(here)
library(readr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(mikelaffr)


# OUTPUT FILES ########################################################################################################
dir.pdf <- here("doc/verifyBamID/pdfs/")

# INPUT FILES #########################################################################################################
smallRNAseq.verifybamid.tsv <- here("results/verifyBamID/20190724_smallRNAseq_verifyBamID_bestSM.tsv")
totalRNAseq.verifybamid.tsv <- here("results/verifyBamID/20190724_totalRNAseq_verifyBamID_bestSM.tsv")

# total metadata
total.meta.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_totalRNAseq_metadata.tsv")
# small metadata
small.meta.tsv <- here("data/metadata/compiled_for_SteinLabHNPs/20190725_smallRNAseq_metadata.tsv")


# Import Data #####################################3
df.small.verifyBamID <- read_tsv(smallRNAseq.verifybamid.tsv)
df.total.verifyBamID <- read_tsv(totalRNAseq.verifybamid.tsv)


# metadata
total.meta <- read_tsv(total.meta.tsv)
small.meta <- read_tsv(small.meta.tsv)


df.total.verifyBamID %<>%
    dplyr::filter(SEQ_ID %in% total.meta$RNAID) %>%
    dplyr::transmute(rnaid = SEQ_ID,
                     donor_id.VBID = sapply(strsplit(CHIP_ID, "_"), `[`, 1),
                     dnaid.VBID = sapply(strsplit(CHIP_ID, "_"), `[`, 2),
                     VBID.FREEMIX = FREEMIX,
                     VBID.CHIPMIX = CHIPMIX,
                     VBID.FREEMIX_CAUTION = FREEMIX > 0.02,
                     VBID.CHIPMIX_CAUTION = CHIPMIX > 0.02)

df.total.verifyBamID$VerifiedDNAID <- ifelse(test = df.total.verifyBamID$VBID.FREEMIX_CAUTION,
                                             yes = "MIXTURE",
                                             no = ifelse(test = df.total.verifyBamID$VBID.CHIPMIX_CAUTION,
                                                         yes = "MISSING",
                                                         no = df.total.verifyBamID$dnaid.VBID))

df.total.verifyBamID$DonorID.via.VerifiedDNAID <- ifelse(test = df.total.verifyBamID$VerifiedDNAID == df.total.verifyBamID$dnaid.VBID,
                                                         yes = df.total.verifyBamID$donor_id.VBID,
                                                         no = NA)
df.small.verifyBamID %<>%
    dplyr::transmute(rnaid = SEQ_ID,
                     donor_id.VBID = sapply(strsplit(CHIP_ID, "_"), `[`, 1),
                     dnaid.VBID = sapply(strsplit(CHIP_ID, "_"), `[`, 2),
                     VBID.FREEMIX = FREEMIX,
                     VBID.CHIPMIX = CHIPMIX,)


df.total.verifyBamID %>%
    ggplot(aes(x=reorder(rnaid, VBID.FREEMIX), y=100*VBID.FREEMIX, fill=VBID.FREEMIX_CAUTION)) +
    geom_col() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    plotTheme("presentation") +
    labs(title = "VerifyBamID Mixture Estimation (sequence only)",
         y = "FREEMIX (Estimated % Mixture)",
         x = "260 Fetal Cortical Wall Samples",
         fill = "CAUTION") +
    scale_y_continuous(expand=expand_scale(mult = c(0.0,0.08))) +
    scale_fill_manual(values = c("grey50", "red")) +
    geom_hline(yintercept = 2)

ggsave(paste(dir.pdf, "vbid_freemix.pdf", sep="/"), height = 6, width = 14)

df.total.verifyBamID %>%
    ggplot(aes(x=reorder(rnaid, VBID.CHIPMIX), y=VBID.CHIPMIX, fill=VBID.CHIPMIX_CAUTION)) +
    geom_col() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "bottom") +
    plotTheme("presentation") +
    labs(title = "VerifyBamID Mixture/Swap Estimation (genotypes + sequence)",
         y = "CHIPMIX (Prob. of Mixture/Swap)",
         x = "260 Fetal Cortical Wall Samples",
         fill = "CAUTION") +
    scale_y_continuous(expand=expand_scale(mult = c(0.0,0.08))) +
    scale_fill_manual(values = c("grey50", "red")) +
    geom_hline(yintercept = 0.02)

ggsave(paste(dir.pdf, "vbid_chipmix.pdf", sep="/"), height = 6, width = 14)




df.small.verifyBamID %>%
    ggplot(aes(x=reorder(rnaid, VBID.FREEMIX), y=100*VBID.FREEMIX)) +
    geom_col() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    plotTheme("presentation") +
    labs(title = "VerifyBamID Mixture Estimation (sequence only)",
         y = "FREEMIX (Estimated % Mixture)",
         x = "240 Fetal Cortical Wall Samples") +
    scale_y_continuous(expand=expand_scale(mult = c(0.0,0.08))) +
    geom_hline(yintercept = 2)

ggsave(paste(dir.pdf, "vbid_freemix_small.pdf", sep="/"), height = 6, width = 14)

df.small.verifyBamID %>%
    ggplot(aes(x=reorder(rnaid, VBID.CHIPMIX), y=VBID.CHIPMIX)) +
    geom_col() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank()) +
    plotTheme("presentation") +
    labs(title = "VerifyBamID Mixture/Swap Estimation (genotypes + sequence)",
         y = "CHIPMIX (Prob. of Mixture/Swap)",
         x = "240 Fetal Cortical Wall Samples") +
    scale_y_continuous(expand=expand_scale(mult = c(0.0,0.08))) +
    geom_hline(yintercept = 0.02)

ggsave(paste(dir.pdf, "vbid_chipmix_small.pdf", sep="/"), height = 6, width = 14)







