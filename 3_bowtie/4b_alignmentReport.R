######
# NGmerge Summary Report
# Ellen Ketter
# March 2023
#######
library(tidyr)
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/4_bowtie/slurm")
slurmFiles <- list.files()
NGmergeReport <- data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:length(slurmFiles)){
  content <- readLines(slurmFiles[i])
  info <- c(slurmFiles[i],content[15],NA)
  NGmergeReport<- rbind(info,NGmergeReport)
}
colnames(NGmergeReport)<- c("slurm","alignmentRate","sample")
NGmergeReport$alignmentRate <- as.numeric(gsub("% overall alignment rate", "",NGmergeReport$alignmentRate))
NGmergeReport$sample[12:17] <- c("Fen2", "HC3","5F2", "PBSa3", "PBSa2",
                                 "PBSa1","PBS2","PBS3","PBS1","Nerol2",
                                 "Nerol3","Fluni3","Nerol1","Myr3","Myr2",
                                 "HQ2","HC1","Myr1","Fluni1","BG3","Fen3",
                                 "Fluni3","HC2","BG1","HQ1","Fen1","5F1")
orderedSlurmFiles <- c("874","739","255","389","388","387",
                       "385","386","384","382","304","383",
                       "381","380","377","379","305","302",
                       "378","268","271","374","303","267",
                       "376","269","855","762","470")  

jpeg(file="../alignmentRate.jpeg")
par(mar=c(5,11,1,4))
ratioRemovedPlot <- barplot(NGmergeReport$alignmentRate, horiz = T, xlab = "Overall Alignment Rate Per Sample",
        names = NGmergeReport$slurm, las=1, col="dark blue")
dev.off()