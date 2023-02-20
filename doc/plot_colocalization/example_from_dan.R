## This code is for ploting the overlap of caQTLs and GWAS hits
library(GenomicRanges);
library(GenomicFeatures);
library(Gviz);
library(rtracklayer);
library(GenomicAlignments);
library(BSgenome.Hsapiens.UCSC.hg38);
library(biomaRt);
options(stringsAsFactors=FALSE);
##Contains the locations of the BAM files
fsheet = '/proj/steinlab/projects/R00/ATACpreprocess/DifferentChromAccess/ATACsampleinfo_noblacklist_CC4toCC20_remvdup.csv';
samples = read.csv(fsheet);
nsamples = dim(samples)[1];
##Load cqn normalizations
##Contains the normalization factors of the ATAC-seq data
fcqn = '/proj/steinlab/projects/R00/ATACpreprocess/DifferentChromAccess/CSAW/Invitro/DESeq2/CSAW_1L_cqnNormFactors.Rdata';
load(fcqn);
normfactors = colMeans(cqnNormFactors);
## verified Bam IDs
verifiedBamIDs=read.csv("/proj/steinlab/projects/R00/ATACpreprocess/VerifyBamID/VerifyBamIDList.csv");
verifiedBamIDs=verifiedBamIDs[which(verifiedBamIDs$FREEMIX<=0.02),];
NverifiedBamIDs=verifiedBamIDs[which(verifiedBamIDs$CellType=="Neuron"),];
PverifiedBamIDs=verifiedBamIDs[which(verifiedBamIDs$CellType=="Progenitor"),];
## Alzheimer's disease
load("OverlappedSNPs_Neuron.Rdata");
resultGR=GRanges(output$chrnames,IRanges(output$peakstart,output$peakend));
mcols(resultGR)=data.frame(CHR_A=output$CHR_A,BP_A=output$BP_A,SNP_A=output$SNP_A);
resultGR=unique(resultGR);
##find rsid for caSNPs
NcaQTLs=read.csv("/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Neuron/CSAW/PCAs_7/NoClumpedGraph/Neuron_all_caQTLs_PC1_7_wA1A2.csv");
resultGR$rsid=NcaQTLs$rsid[match(resultGR$SNP_A,NcaQTLs$SNP)];
## dir genotype data
Ngenodir="/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/TruncatedFiles/Neuron/";
Pgenodir="/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/TruncatedFiles/Progenitor/";
## dir EMMAX caQTL results
Nemmax="/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Neuron/CSAW/PCAs_7/Neuron/RawData/";
Pemmax="/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Progenitor/CSAW/PCAs_4/Progenitor/RawData/";
##GWAS
GWASdata_1=read.table("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/AD/Kunkle_etal_Stage1_results.txt?file=1",header=T);
GWASdata_2=read.table("/proj/steinlab/projects/R00/ATACPartitionedHeritability/diseasetraitmunged/origdownloads/AD/Kunkle_etal_Stage2_results.txt?file=1",header=T);
GWASdata=rbind(GWASdata_1,GWASdata_2);
## make data a GRanges
GWASdataGR = GRanges(paste0("chr",GWASdata$Chromosome),IRanges(GWASdata$Position,GWASdata$Position));
mcols(GWASdataGR)=GWASdata;
## convert hg19 to hg38
path="/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/DiseaseColocalization/CoordinatesColocalization/ENIGMA3Traits/hg19ToHg38.over.chain";
ch = import.chain(path);
GWASsnpGR=unlist(liftOver(GWASdataGR,ch));
seqlevels(GWASsnpGR) = sub("chr", "", seqlevels(GWASsnpGR));
GWASsnpGR$P=-log10(as.numeric(GWASsnpGR$Pvalue));
totalGWASsnpGR=GWASsnpGR;
genemart=useEnsembl(biomart = "ensembl",dataset = "hsapiens_gene_ensembl");
pdf('AD_NcaQTLs.pdf');
for (n in 1:length(resultGR)){
    ##genotype IDs
    NgenoIDs=read.table(paste0(Ngenodir,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],".tfam"),header=F);
    PgenoIDs=read.table(paste0(Pgenodir,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],".tfam"),header=F);
    NIDs=paste0(NgenoIDs$V1,"_",NgenoIDs$V2);
    PIDs=paste0(PgenoIDs$V1,"_",PgenoIDs$V2);
    NverifiedBamIDs=NverifiedBamIDs[match(NIDs,NverifiedBamIDs$MatchedDNAID),];
    PverifiedBamIDs=PverifiedBamIDs[match(PIDs,PverifiedBamIDs$MatchedDNAID),];
    #uniqbeta=unique(thissnpGR$Beta);
    ##Set parameters to read in only the chromosome of interest for ATAC
    ##Set a region extend to upstream and downstream extendsize bp
    upextendsize = 100000;
    downextendsize = 100000;
    param = ScanBamParam(what=c("pos","qwidth"),which=GRanges(paste0("chr",seqnames(resultGR)[n]), IRanges(start(resultGR)[n]-upextendsize,end(resultGR)[n]+downextendsize)),flag=scanBamFlag(isUnmappedQuery=FALSE));
    areawidth = width(GRanges(paste0("chr",seqnames(resultGR)[n]), IRanges(start(resultGR)[n]-upextendsize,end(resultGR)[n]+downextendsize)));

    ##Set an empty values file
    values = matrix(NA,nrow=nsamples,ncol=(end(resultGR)[n]-start(resultGR)[n]+1+upextendsize+downextendsize));
    for (i in 1:nsamples){
        cat('Working on loading bamfile...',i,'\n');
        ##Read in the bam file
        bamFile = samples$bamReads[i];
        bam = readGAlignments(bamFile,param=param);
        ##change bam file to GRanges to shift the reads
        bamgrange=granges(bam);
        #all reads aligning to + strand offset by +4bp
        bamgrange[strand(bamgrange)=="+"]=shift(bamgrange[strand(bamgrange)=="+"],4);
        #all reads aligning to - strand offset by -5bp
        bamgrange[strand(bamgrange)=="-"]=shift(bamgrange[strand(bamgrange)=="-"],-5);
        ##Get the coverage values
        bamcoverage = coverage(bamgrange);
        chrind = which(names(bamcoverage) == paste0("chr",seqnames(resultGR)[n]));
        ##Save coverage values
        values[i,] = as.numeric(bamcoverage[chrind][[1]])[(start(resultGR)[n]-upextendsize):(end(resultGR)[n]+downextendsize)];

    }
    ##Take an average of every 40 bp to reduce amount of data
    chunksize = 50;
    areawidthreduced = floor(areawidth/chunksize);
    valuesreduced = matrix(NA,nrow=nsamples,ncol=areawidthreduced);
    for (j in 1:areawidthreduced) {
        indstoavg = (j*chunksize-chunksize+1):(j*chunksize);
        valuesreduced[,j] = rowMeans(values[,indstoavg]);
    }
    ##Divide the values by the normalization factor
    valuesreduced = apply(valuesreduced,2,function (x) x / normfactors);
    ## record every track in a track list
    ## Load chromosome track
    cat('Working on adding itrack...\n');
#    itrack = IdeogramTrack(genome = "hg38", chromosome = paste0("chr",seqnames(resultGR)[n]));
#    trackplot = c(itrack);
    trackplot = c();
    ## Load axis track
    cat('Working on adding axistrack...\n');
    gtrack = GenomeAxisTrack();
    trackplot = append(trackplot, gtrack);
    genehg38 = BiomartGeneRegionTrack(genome="hg38", biomart=genemart, chromosome=paste0('chr',seqnames(resultGR[n])), start= start(resultGR)[n]-upextendsize, end= end(resultGR)[n]+downextendsize, sshowId=TRUE, geneSymbols=TRUE,transcriptAnnotation="symbol",name="ENSEMBL_hg38");
    trackplot = append(trackplot,genehg38);

    ## SNP pvalue track
    ## Emmax caQTL -log10(p-value)
    Nsnplist=read.table(paste0(Nemmax,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],"_Neuron.ps"),header=T);
    Ndonor_n=read.table(paste0("/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/",seqnames(resultGR)[n],".dose.R2g03.QC.Neuron.frqx"),skip=1);
    ##keep SNPs with at least 2 donors in minor allele hom or het
    NkeepSNPs=Ndonor_n$V2[which(Ndonor_n$V5 >=2 | Ndonor_n$V6 >=2)];
    Nsnplist=Nsnplist[which(Nsnplist$SNP %in% NkeepSNPs),];
    SNPseqs=sapply(Nsnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[1]);
    SNPcoods=sapply(Nsnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[2]);
    NSNPGR = GRanges(seqnames=SNPseqs,IRanges(as.numeric(SNPcoods),as.numeric(SNPcoods)));
    mcols(NSNPGR) = data.frame(SNP=Nsnplist$SNP,P=-log10(Nsnplist$P));

    Psnplist=read.table(paste0(Pemmax,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],"_Progenitor.ps"),header=T);
    Pdonor_n=read.table(paste0("/proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/",seqnames(resultGR)[n],".dose.R2g03.QC.Progenitor.frqx"),skip=1);
    ##keep SNPs with at least 2 donors in minor allele hom or het
    PkeepSNPs=Pdonor_n$V2[which(Pdonor_n$V5 >=2 | Pdonor_n$V6 >=2)];
    Psnplist=Psnplist[which(Psnplist$SNP %in% PkeepSNPs),];
    SNPseqs=sapply(Psnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[1]);
    SNPcoods=sapply(Psnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[2]);
    PSNPGR = GRanges(seqnames=SNPseqs,IRanges(as.numeric(SNPcoods),as.numeric(SNPcoods)));
    mcols(PSNPGR) = data.frame(SNP=Psnplist$SNP,P=-log10(Psnplist$P));
    ##get the max -log10(P)
    SNPymax=max(mcols(PSNPGR)$P,mcols(NSNPGR)$P);

    ##Neuron LD
    system(paste0("module add plink/1.90b3;plink --bfile /proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/",paste0(seqnames(resultGR)[n]),".dose.R2g03.QC.IDfixed --keep /proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/TruncatedFiles/KeepIds_Neuron.txt --r2 --ld-window-r2 0.2 --ld-window-kb 10000 --ld-window 2000 --ld-snp ",resultGR$SNP_A[n]," --out Nplink"));
    LD=read.table("Nplink.ld",header=TRUE);
    NLD=LD;
    ##In order to color the different SNP dots differently, need to make overlay tracks for each possibility
    ldmatchind = match(LD$SNP_B,NSNPGR$SNP);
    ldsnpsnocolor = setdiff(1:length(NSNPGR),ldmatchind[which(!is.na(ldmatchind))]);
    NSNPGRnocolor = NSNPGR[ldsnpsnocolor];
    ##Add -log10P
    mcols(NSNPGRnocolor) = data.frame(log10P=mcols(NSNPGRnocolor)$P);    
    ##Make a data track for each SNP in the region as a p-value
    ntracknocolor <- DataTrack(NSNPGRnocolor,chromosome=paste0(seqnames(resultGR)[n]), name = "Neuron caQTLs \n -log10(P)",type="p",legend=FALSE,col="#192752",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.8)],NSNPGR$SNP);
    NSNPGRred = NSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add -log10P
    mcols(NSNPGRred) = data.frame(log10P=mcols(NSNPGRred)$P);
    ntrackred <- DataTrack(NSNPGRred, chromosome=paste0(seqnames(resultGR)[n]),name = "Neuron caQTLs \n -log10(P)",type="p",legend=FALSE,col="red",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.6 & LD$R2 <0.8)],NSNPGR$SNP);
    NSNPGRorange = NSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add -log10P
    mcols(NSNPGRorange) = data.frame(log10P=mcols(NSNPGRorange)$P);
    ntrackorange <- DataTrack(NSNPGRorange, chromosome=paste0(seqnames(resultGR)[n]),name = "Neuron caQTLs \n -log10(P)",type="p",legend=FALSE,col="orange",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.4 & LD$R2 <0.6)],NSNPGR$SNP);
    NSNPGRgreen = NSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(NSNPGRgreen) = data.frame(log10P=mcols(NSNPGRgreen)$P);
    ntrackgreen <- DataTrack(NSNPGRgreen, chromosome=paste0(seqnames(resultGR)[n]),name = "Neuron caQTLs \n -log10(P)",type="p",legend=FALSE,col="green",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.2 & LD$R2 <0.4)],NSNPGR$SNP);
    NSNPGRlightblue = NSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(NSNPGRlightblue) = data.frame(log10P=mcols(NSNPGRlightblue)$P);
    ntracklightblue <- DataTrack(NSNPGRlightblue, chromosome=paste0(seqnames(resultGR)[n]),name = "Neuron caQTLs \n -log10(P)",type="p",legend=FALSE,col="lightblue",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");
    ##start making a plot in Gviz
    NSNP_track = OverlayTrack(trackList = list(ntracknocolor,ntracklightblue,ntrackgreen,ntrackorange,ntrackred));

    ##ProgenitorLD
    system(paste0("module add plink/1.90b3;plink --bfile /proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/",paste0(seqnames(resultGR)[n]),".dose.R2g03.QC.IDfixed --keep /proj/steinlab/projects/R00/atac-qtl/EMMAXResult/Genotype_Hg38/TruncatedFiles/KeepIds_Progenitor.txt --r2 --ld-window-r2 0.2 --ld-window-kb 10000 --ld-window 2000 --ld-snp ",resultGR$SNP_A[n]," --out Nplink.ld"));
    LD=read.table("Nplink.ld",header=TRUE);
    ##In order to color the different SNP dots differently, need to make overlay tracks for each possibility
    ldmatchind = match(LD$SNP_B,PSNPGR$SNP);
    ldsnpsnocolor = setdiff(1:length(PSNPGR),ldmatchind[which(!is.na(ldmatchind))]);
    PSNPGRnocolor = PSNPGR[ldsnpsnocolor];
    ##Add chr to seqnames
    mcols(PSNPGRnocolor) = data.frame(log10P=mcols(PSNPGRnocolor)$P);
    ##Make a data track for each SNP in the region as a p-value
    ptracknocolor <- DataTrack(PSNPGRnocolor, chromosome=paste0(seqnames(resultGR)[n]),name = "Progenitor caQTLs \n -log10(P)",type="p",legend=FALSE,col="#192752",ylim=c(0,SNPymax*1.1),baseline=-log10(2.9e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.8)],PSNPGR$SNP);
    PSNPGRred = PSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(PSNPGRred) = data.frame(log10P=mcols(PSNPGRred)$P);
    ptrackred <- DataTrack(PSNPGRred,chromosome=paste0(seqnames(resultGR)[n]) ,name = "Progenitor caQTLs \n -log10(P)",type="p",legend=FALSE,col="red",ylim=c(0,SNPymax*1.1),baseline=-log10(2.9e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.6 & LD$R2 <0.8)],PSNPGR$SNP);
    PSNPGRorange = PSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(PSNPGRorange) = data.frame(log10P=mcols(PSNPGRorange)$P);
    ptrackorange <- DataTrack(PSNPGRorange,chromosome=paste0(seqnames(resultGR)[n]) ,name = "Progenitor caQTLs \n -log10(P)",type="p",legend=FALSE,col="orange",ylim=c(0,SNPymax*1.1),baseline=-log10(2.9e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.4 & LD$R2 <0.6)],PSNPGR$SNP);
    PSNPGRgreen = PSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(PSNPGRgreen) = data.frame(log10P=mcols(PSNPGRgreen)$P);
    ptrackgreen <- DataTrack(PSNPGRgreen, chromosome=paste0(seqnames(resultGR)[n]),name = "Progenitor caQTLs \n -log10(P)",type="p",legend=FALSE,col="green",ylim=c(0,SNPymax*1.1),baseline=-log10(2.9e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = match(LD$SNP_B[which(LD$R2>=0.2 & LD$R2 <0.4)],PSNPGR$SNP);
    PSNPGRlightblue = PSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    mcols(PSNPGRlightblue) = data.frame(log10P=mcols(PSNPGRlightblue)$P);
    ptracklightblue <- DataTrack(PSNPGRlightblue,chromosome=paste0(seqnames(resultGR)[n]) ,name = "Progenitor caQTLs \n -log10(P)",type="p",legend=FALSE,col="lightblue",ylim=c(0,SNPymax*1.1),baseline=-log10(2.9e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ##start making a plot in Gviz
    PSNP_track = OverlayTrack(trackList = list(ptracknocolor,ptracklightblue,ptrackgreen,ptrackorange,ptrackred));
    ##GWAS track
    thisresultGR=GRanges(seqnames(resultGR)[n],IRanges(start(resultGR)[n]-upextendsize,end(resultGR)[n]+downextendsize));
    GWASsnpGR=totalGWASsnpGR[queryHits(findOverlaps(totalGWASsnpGR,thisresultGR))];
    pymax=1.2*max(GWASsnpGR$P);
    ##1000 Genome LD
    system(paste0("module add plink/1.90b3;plink --bfile /proj/steinlab/projects/1000genomes/phase3EURhg38/ALL.chr",paste0(seqnames(resultGR)[n]),"_GRCh38.genotypes.20170504.EUR --r2 --ld-window-r2 0 --ld-window-kb 10000 --ld-window 2000 --ld-snp ",resultGR$rsid[n]," --out Nplink"));
    G1000LD=read.table("Nplink.ld",header=TRUE);
    SNP_BGR=GRanges(G1000LD$CHR_B,IRanges(G1000LD$BP_B,G1000LD$BP_B));
    mcols(SNP_BGR)=data.frame(rsid=G1000LD$SNP_B,R2=G1000LD$R2);
    ldsnpsna=setdiff(1:length(GWASsnpGR),queryHits(findOverlaps(GWASsnpGR,SNP_BGR)));
    GWASsnpGRna=GWASsnpGR[ldsnpsna];
    ##Add chr to seqnames
    GWASsnpGRna = renameSeqlevels(GWASsnpGRna, paste0("chr",seqlevels(GWASsnpGRna)));
    ##Make a data track for each SNP in the region as a p-value
    mcols(GWASsnpGRna) = data.frame(log10P=mcols(GWASsnpGRna)$P);
    dtrackna <- DataTrack(GWASsnpGRna,chromosome=G1000LD$CHR_B[1],name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="#192752",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");
    ldsnpsnocolor = unique(queryHits(findOverlaps(GWASsnpGR,SNP_BGR[which(SNP_BGR$R2<0.2)])));
    GWASsnpGRnocolor = GWASsnpGR[ldsnpsnocolor];
    ##Add chr to seqnames
    GWASsnpGRnocolor = renameSeqlevels(GWASsnpGRnocolor, paste0("chr",seqlevels(GWASsnpGRnocolor)));
    ##Make a data track for each SNP in the region as a p-value
    mcols(GWASsnpGRnocolor) = data.frame(log10P=mcols(GWASsnpGRnocolor)$P);

    dtracknocolor <- DataTrack(GWASsnpGRnocolor,chromosome=G1000LD$CHR_B[1],name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="#192752",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = unique(queryHits(findOverlaps(GWASsnpGR,SNP_BGR[which(SNP_BGR$R2 >= 0.8)])));
    GWASsnpGRred = GWASsnpGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    GWASsnpGRred = renameSeqlevels(GWASsnpGRred, paste0("chr",seqlevels(GWASsnpGRred)));
    mcols(GWASsnpGRred) = data.frame(log10P=mcols(GWASsnpGRred)$P);
    dtrackred <- DataTrack(GWASsnpGRred,chromosome=G1000LD$CHR_B[1], name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="red",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = unique(queryHits(findOverlaps(GWASsnpGR,SNP_BGR[which(SNP_BGR$R2 < 0.8 & SNP_BGR$R2 >=0.6)])));
    GWASsnpGRorange = GWASsnpGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    GWASsnpGRorange = renameSeqlevels(GWASsnpGRorange, paste0("chr",seqlevels(GWASsnpGRorange)));
    mcols(GWASsnpGRorange) = data.frame(log10P=mcols(GWASsnpGRorange)$P);
    dtrackorange <- DataTrack(GWASsnpGRorange, chromosome=G1000LD$CHR_B[1],name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="orange",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = unique(queryHits(findOverlaps(GWASsnpGR,SNP_BGR[which(SNP_BGR$R2 < 0.6 & SNP_BGR$R2 >=0.4)])));
    GWASsnpGRgreen = GWASsnpGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    GWASsnpGRgreen = renameSeqlevels(GWASsnpGRgreen, paste0("chr",seqlevels(GWASsnpGRgreen)));
    mcols(GWASsnpGRgreen) = data.frame(log10P=mcols(GWASsnpGRgreen)$P);
    dtrackgreen <- DataTrack(GWASsnpGRgreen,chromosome=G1000LD$CHR_B[1], name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="green",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ldmatchind = unique(queryHits(findOverlaps(GWASsnpGR,SNP_BGR[which(SNP_BGR$R2 < 0.4 & SNP_BGR$R2 >=0.2)])));
    GWASsnpGRlightblue = GWASsnpGR[ldmatchind[which(!is.na(ldmatchind))]];
    ##Add chr to seqnames
    GWASsnpGRlightblue = renameSeqlevels(GWASsnpGRlightblue, paste0("chr",seqlevels(GWASsnpGRlightblue)));
    mcols(GWASsnpGRlightblue) = data.frame(log10P=mcols(GWASsnpGRlightblue)$P);
    dtracklightblue <- DataTrack(GWASsnpGRlightblue,chromosome=G1000LD$CHR_B[1], name = "Alzheimer's Disease \n -log10(P)",type="p",legend=FALSE,col="lightblue",ylim=c(0,pymax),baseline=-log10(5e-8),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

    ##start making a plot in Gviz
    GWAS_track = OverlayTrack(trackList = list(dtrackna,dtracknocolor,dtracklightblue,dtrackgreen,dtrackorange,dtrackred));
    cat('Working on adding differential peaks annotation track...\n');
    op_annottrack = AnnotationTrack(start=start(resultGR)[n],width=width(resultGR)[n],name="Peak",genome = "hg38", chromosome=paste0("chr",seqnames(resultGR)[n]),fill="deepskyblue1");
    ##Create a series of scaled bamtracks for plotting averaging across all donors of the same condition
    snpGR=GRanges(resultGR$CHR_A[n],IRanges(resultGR$BP_A[n],resultGR$BP_A[n]));
    mcols(snpGR)=mcols(resultGR)[n,];
    snpop_annottrack = AnnotationTrack(start=start(snpGR),width=width(snpGR),name="SNP",genome = "hg38", chromosome=paste0(seqnames(snpGR)),fill="darkolivegreen3");
    ## Genotypes
    Ngeno=read.table(paste0(Ngenodir,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],".tped"),header=F);
    Ngeno = Ngeno[which(Ngeno$V2==snpGR$SNP_A[1]),];
    NA1A2=paste0(Ngeno[1,seq(5,dim(Ngeno)[2],by=2)],"/",Ngeno[1,seq(6,dim(Ngeno)[2],by=2)]);
    Pgeno=read.table(paste0(Pgenodir,"chr",seqnames(resultGR)[n],"/chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],".tped"),header=F);
    Pgeno = Pgeno[which(Pgeno$V2==snpGR$SNP_A[1]),];
    PA1A2=paste0(Pgeno[1,seq(5,dim(Pgeno)[2],by=2)],"/",Pgeno[1,seq(6,dim(Pgeno)[2],by=2)]);
    ## find the max of y-axis
    startseq = seq(start(resultGR)[n]-upextendsize,(areawidthreduced-1)*chunksize+start(resultGR)[n]-upextendsize,by=chunksize);
    endseq = seq(start(resultGR)[n]+chunksize-upextendsize,areawidthreduced*chunksize+start(resultGR)[n]-upextendsize,by=chunksize);
    progvals = valuesreduced[match(PverifiedBamIDs$Names,samples$SampleID),];
    progvals1 = colMeans(progvals[which(PA1A2=="1/1"), ,drop = FALSE]);
    progvals2 = colMeans(progvals[which(PA1A2=="1/2"|PA1A2=="2/1"), ,drop = FALSE]);
    progvals3 = colMeans(progvals[which(PA1A2=="2/2"), ,drop = FALSE]);
    neurvals = valuesreduced[match(NverifiedBamIDs$Names,samples$SampleID),];
    neurvals1 = colMeans(neurvals[which(NA1A2=="1/1"), ,drop = FALSE]);
    neurvals2 = colMeans(neurvals[which(NA1A2=="1/2"|NA1A2=="2/1"), ,drop = FALSE]);
    neurvals3 = colMeans(neurvals[which(NA1A2=="2/2"), ,drop = FALSE]);
    ymax = max(c(progvals1,progvals2,progvals3,neurvals1,neurvals2,neurvals3)[which(!is.na(c(progvals1,progvals2,progvals3,neurvals1,neurvals2,neurvals3)))]);
    ##progenitor tracks
    progtrackatac1 = DataTrack(start=startseq,end=endseq,data=progvals1,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)), lwd.baseline=2);
    progtrackatac2 = DataTrack(start=startseq,end=endseq,data=progvals2,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)),lwd.baseline=2);
    progtrackatac3 = DataTrack(start=startseq,end=endseq,data=progvals3,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)),lwd.baseline=2);
    displayPars(progtrackatac1) <- list(alpha.title = 1, alpha = 0.75,col="coral4");
    displayPars(progtrackatac2) <- list(alpha.title = 1, alpha = 0.75,col="green3");
    displayPars(progtrackatac3) <- list(alpha.title = 1, alpha = 0.75,col="darkblue");
    progot <- OverlayTrack(trackList = list(progtrackatac1,progtrackatac2,progtrackatac3));
    ##neuron tracks
    neurtrackatac1 = DataTrack(start=startseq,end=endseq,data=neurvals1,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
    neurtrackatac2 = DataTrack(start=startseq,end=endseq,data=neurvals2,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
    neurtrackatac3 = DataTrack(start=startseq,end=endseq,data=neurvals3,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
    displayPars(neurtrackatac1) <- list(alpha.title = 1, alpha = 0.75,col="coral4");
    displayPars(neurtrackatac2) <- list(alpha.title = 1, alpha = 0.75,col="green3");
    displayPars(neurtrackatac3) <- list(alpha.title = 1, alpha = 0.75,col="darkblue");
    neurot <- OverlayTrack(trackList = list(neurtrackatac1,neurtrackatac2,neurtrackatac3));
    ##In order to color the different SNP dots differently, need to make overlay tracks for each possibility
    condNsnpfiles=list.files(".",pattern=glob2rx(paste0("chr",seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],"*_Neuron.ps")));
    for (m in 1:length(condNsnpfiles)){
        condrsid=unlist(strsplit(condNsnpfiles[m],"_",fixed=TRUE))[4];
        condSNPGR=totalGWASsnpGR[which(totalGWASsnpGR$MarkerName==condrsid)];
        ##conditional SNP track
        condsnp_annottrack = AnnotationTrack(start=start(condSNPGR),width=width(condSNPGR),name="conditional \n SNP",genome = "hg38", chromosome=paste0(seqnames(condSNPGR)));    
        ##conditinal caQTL results
        condNsnplist=read.table(condNsnpfiles[m],header=F);
        colnames(condNsnplist)=c("SNP","Beta","P");
        condNsnplist=condNsnplist[which(condNsnplist$SNP %in% NkeepSNPs),];
        SNPseqs=sapply(condNsnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[1]);
        SNPcoods=sapply(condNsnplist$SNP, function(x) unlist(strsplit(x,":",fixed="TRUE"))[2]);
        condNSNPGR = GRanges(seqnames=SNPseqs,IRanges(as.numeric(SNPcoods),as.numeric(SNPcoods)));
        mcols(condNSNPGR) = data.frame(SNP=condNsnplist$SNP,P=-log10(condNsnplist$P));

        ldmatchind = match(NLD$SNP_B,condNSNPGR$SNP);
        ldsnpsnocolor = setdiff(1:length(condNSNPGR),ldmatchind[which(!is.na(ldmatchind))]);
        condNSNPGRnocolor = condNSNPGR[ldsnpsnocolor];
        ##Add -log10P
        mcols(condNSNPGRnocolor) = data.frame(log10P=mcols(condNSNPGRnocolor)$P);
        ##Make a data track for each SNP in the region as a p-value
        ntracknocolor <- DataTrack(condNSNPGRnocolor,chromosome=paste0(seqnames(resultGR)[n]), name = "Cond NcaQTLs \n -log10(P)",type="p",legend=FALSE,col="#192752",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

        ldmatchind = match(NLD$SNP_B[which(NLD$R2>=0.8)],condNSNPGR$SNP);
        condNSNPGRred = condNSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
        ##Add -log10P
        mcols(condNSNPGRred) = data.frame(log10P=mcols(condNSNPGRred)$P);
        ntrackred <- DataTrack(condNSNPGRred, chromosome=paste0(seqnames(resultGR)[n]),name = "Cond NcaQTLs \n -log10(P)",type="p",legend=FALSE,col="red",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

        ldmatchind = match(NLD$SNP_B[which(NLD$R2>=0.6 & NLD$R2 <0.8)],condNSNPGR$SNP);
        condNSNPGRorange = condNSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
        ##Add -log10P
        mcols(condNSNPGRorange) = data.frame(log10P=mcols(condNSNPGRorange)$P);
        ntrackorange <- DataTrack(condNSNPGRorange, chromosome=paste0(seqnames(resultGR)[n]),name = "Cond NcaQTLs \n -log10(P)",type="p",legend=FALSE,col="orange",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

        ldmatchind = match(NLD$SNP_B[which(NLD$R2>=0.4 & NLD$R2 <0.6)],condNSNPGR$SNP);
        condNSNPGRgreen = condNSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
        ##Add chr to seqnames
        mcols(condNSNPGRgreen) = data.frame(log10P=mcols(condNSNPGRgreen)$P);
        ntrackgreen <- DataTrack(condNSNPGRgreen, chromosome=paste0(seqnames(resultGR)[n]),name = "Cond NcaQTLs \n -log10(P)",type="p",legend=FALSE,col="green",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

        ldmatchind = match(NLD$SNP_B[which(NLD$R2>=0.2 & NLD$R2 <0.4)],condNSNPGR$SNP);
        condNSNPGRlightblue = condNSNPGR[ldmatchind[which(!is.na(ldmatchind))]];
        ##Add chr to seqnames
        mcols(condNSNPGRlightblue) = data.frame(log10P=mcols(condNSNPGRlightblue)$P);
        ntracklightblue <- DataTrack(condNSNPGRlightblue, chromosome=paste0(seqnames(resultGR)[n]),name = "Cond NcaQTLs \n -log10(P)",type="p",legend=FALSE,col="lightblue",ylim=c(0,SNPymax*1.1),baseline=-log10(1.5e-5),col.baseline="grey",lty.baseline=2,lwd.baseline=1,genome="hg38");

        ##start making a plot in Gviz
        condNSNP_track = OverlayTrack(trackList = list(ntracknocolor,ntracklightblue,ntrackgreen,ntrackorange,ntrackred));
        ## plot chromosome track, axis track, giene track and data track in one pdf figure

        cat('Working on ploting...\n');
        upextendsize = 100000;
        downextendsize = 100000;
        plotname = paste0(seqnames(resultGR)[n],"_",start(resultGR)[n],"_",end(resultGR)[n],"\n rsid=",snpGR$SNP_A,"\n Condition on ",condrsid); 
        plotTracks(c(trackplot,condsnp_annottrack,GWAS_track,PSNP_track,progot,NSNP_track,neurot,snpop_annottrack,op_annottrack,condNSNP_track),main=plotname,cex.main=1,from = start(resultGR)[n]-upextendsize,to = end(resultGR)[n]+downextendsize,transcriptAnnotation="symbol",add53=TRUE,showBandID=TRUE,cex.bands=0.7,stackHeight=0.8,background.title = "white",col.axis="black",col.title="black",cex.title=0.4,cex.axis=0.3,just.group="below",collapseTranscripts="meta");
#collapseTranscripts="meta"

        plotname = paste0(seqnames(resultGR)[n],"_",start(resultGR)[n]-upextendsize,"_",end(resultGR)[n]+downextendsize,"\n rsid=",snpGR$SNP_A,"\n Zoomed In "," Condition on ",condrsid);
        upextendsize = 10000;
        downextendsize = 10000;
        valuecen=length(progvals1)/2;
        ind=(valuecen-upextendsize/chunksize):(valuecen+downextendsize/chunksize);
        ymax = max(c(progvals1[ind],progvals2[ind],progvals3[ind],neurvals1[ind],neurvals2[ind],neurvals3[ind])[which(!is.na(c(progvals1[ind],progvals2[ind],progvals3[ind],neurvals1[ind],neurvals2[ind],neurvals3[ind])))]);
        ##progenitor track
        progtrackatac1 = DataTrack(start=startseq,end=endseq,data=progvals1,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)), lwd.baseline=2);
        progtrackatac2 = DataTrack(start=startseq,end=endseq,data=progvals2,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)),lwd.baseline=2);
        progtrackatac3 = DataTrack(start=startseq,end=endseq,data=progvals3,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Progenitor ATAC Avg \n','N= ', length(PA1A2)),lwd.baseline=2);
        displayPars(progtrackatac1) <- list(alpha.title = 1, alpha = 0.75,col="coral4");
        displayPars(progtrackatac2) <- list(alpha.title = 1, alpha = 0.75,col="green3");
        displayPars(progtrackatac3) <- list(alpha.title = 1, alpha = 0.75,col="darkblue");
        progot <- OverlayTrack(trackList = list(progtrackatac1,progtrackatac2,progtrackatac3));
        ##neuron track
        neurtrackatac1 = DataTrack(start=startseq,end=endseq,data=neurvals1,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
        neurtrackatac2 = DataTrack(start=startseq,end=endseq,data=neurvals2,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
        neurtrackatac3 = DataTrack(start=startseq,end=endseq,data=neurvals3,chromosome=paste0("chr",seqnames(resultGR)[n]),genome="hg38",type="a",ylim=c(0,ymax),name=paste0('Neuron ATAC Avg \n','N= ', length(NA1A2)),lwd.baseline=2);
        displayPars(neurtrackatac1) <- list(alpha.title = 1, alpha = 0.75,col="coral4");
        displayPars(neurtrackatac2) <- list(alpha.title = 1, alpha = 0.75,col="green3");
        displayPars(neurtrackatac3) <- list(alpha.title = 1, alpha = 0.75,col="darkblue");
        neurot <- OverlayTrack(trackList = list(neurtrackatac1,neurtrackatac2,neurtrackatac3));
        ##plot
        plotTracks(c(trackplot,condsnp_annottrack,GWAS_track,PSNP_track,progot,NSNP_track,neurot,snpop_annottrack,op_annottrack,condNSNP_track),main=plotname,cex.main=1,from = start(resultGR)[n]-upextendsize,to = end(resultGR)[n]+downextendsize,transcriptAnnotation="symbol",add53=TRUE,showBandID=TRUE,cex.bands=0.7,stackHeight=0.8,background.title = "white",col.axis="black",col.title="black",cex.title=0.4,cex.axis=0.3,just.group="below");
}
}
#}
dev.off();

