# library(BiocManager)
# BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments"))
# BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
### # BiocManager::install(c("phastCons100way.UCSC.mm10"))????

library(ATACseqQC)
library(stringr)
## input the bamFile from the ATACseqQC package 
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/6_ATACseqQC/")
source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))
# Estimate the library complexity
bamfile <- list.files()
bamfile <- str_subset(bamfile, "bai", negate = TRUE)
bamfile <- str_subset(bamfile, "preMarkdup")

bamfileLabels <- gsub(".bam","",bamfile)

possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                  "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                  "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                  "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                  "U2"))

for(i in 1:length(bamfile)){
  bamTop100 <- scanBam(BamFile(bamfile[2], yieldSize = 100),
                       param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
  tags <- names(bamTop100)[lengths(bamTop100)>0]
  print(tags)
  ## files will be output into outPath
  outPath <- "shifted_chr1/"
  # dir.create(outPath)
  
  ## shift the coordinates of 5'ends of alignments in the bam file
  seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)
  seqlev <- "chr1"
  which <- as(seqinformation[seqlev], "GRanges")
  gal <- readBamFile(bamfile[2], tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
  shiftedBamfile <- file.path(outPath, paste0(bamfileLabels[2],".shifted.bam"))
  gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
  txs <- transcripts(TxDb.Mmusculus.UCSC.mm10.knownGene)
  
  pdf(paste0(bamfileLabels[2],"_plots.pdf"))
  pt <- PTscore(gal1, txs)
  plot(pt$log2meanCoverage, pt$PT_score, 
       xlab="log2 mean coverage",
       ylab="Promoter vs Transcript")
  
  nfr <- NFRscore(gal1, txs)
  plot(nfr$log2meanCoverage, nfr$NFR_score, 
       xlab="log2 mean coverage",
       ylab="Nucleosome Free Regions score",
       main="NFRscore for 200bp flanking TSSs",
       xlim=c(-10, 0), ylim=c(-5, 5))
  
  tsse <- TSSEscore(gal1, txs)
  tsse$TSSEscore
  dev.off
}







for (i in 1:length(bamfile)){ 
  print(i)
  bamQC <- bamQC(paste0("../5_removeDuplicates/",bamfile[i]), 
        outPath=NULL) 
  saveRDS(bamQC,paste0(bamfileLabels[i],".rds"))
  }
# took about 3 hours to run
for (i in 1:length(bamfile)){ 
  jpeg(file=paste0(bamfileLabels[i],"_fragSize.jpeg"))
  fragSize <- fragSizeDist(bamfile[i], bamfile.labels[i])
  dev.off() 
}

for (i in 1:length(bamfile)){ 
  jpeg(file=paste0(bamfileLabels[i],"_mapQ.jpeg"))
  tmp <- readRDS(paste0(bamfileLabels[i],".rds")) 
  mapQ <- tmp$MAPQ
  mapQ$Var1<- as.numeric(mapQ$Var1)
  mapQ<- mapQ[order(mapQ$Var1),]
  barplot(mapQ$Freq, names = mapQ$Var1, log = "y",
          xlab = "MAPQ Score", ylab = "Number of Reads",
          col="#69b3a2", main = "Read Frequency per MAPQ Score",
          sub = bamfileLabels[i])
  dev.off() 
}
samplesPercentMito <- c()
for (i in 1:length(bamfile)){ 
  tmp <- readRDS(paste0("ATACseqQC/",bamfileLabels[i],".rds")) 
  idxstats<- `tmp`$idxstats
  NuclearReadsMapped<- sum(idxstats$mapped[1:20])
  mitoReadsMapped <- idxstats$mapped[22]
  percentMito <- (mitoReadsMapped/(NuclearReadsMapped+mitoReadsMapped))*100
  samplesPercentMito <- append(percentMito,samplesPercentMito)
}

samplesPercentMito <- readRDS("ATACseqQC/samplesPercentMT.rds")
samplesPercentMito <- as.data.frame(samplesPercentMito)
treatments <- c("5F", "5F", "BG", "BG", 
                "Fen", "Fen", "Fen", 
                "Fluni", "Fluni", "Fluni", 
                "HC",  "HC",  "HC", 
                "HQ", "HQ",
                "Myr", "Myr", "Myr", 
                "Nerol", "Nerol", "Nerol",
                "PBS", "PBS", "PBS", "PBS", "PBS", "PBS")
samplesPercentMito <- cbind(treatments, samplesPercentMito)
colnames(samplesPercentMito) <- c("treatments", "percent")
jpeg(file="percentMT.jpeg")

# Libraries
library(tidyverse)
library(hrbrthemes)
library(viridis)

# create a dataset
samplesPercentMito %>%
  ggplot( aes(x=treatments, y=percent, fill=treatments)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    theme_ipsum() +
    theme(
      legend.position="none",
     plot.title = element_text(size=11)
    ) +
    ggtitle("Percent Mitochondrial DNA per Sample") +
    xlab("") +
    ylim(0,10)


dev.off() 