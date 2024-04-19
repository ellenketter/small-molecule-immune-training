######
# NGmerge Summary Report
# Ellen Ketter
# March 2023
#######

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/3_NGmerge/slurm/")
slurmFiles <- list.files()
NGmergeReport <- data.frame(matrix(ncol = 3, nrow = 0))
for (i in 1:length(slurmFiles)){
  content <- readLines(slurmFiles[i])
  info <- c(substring(content[1],19,37), substring(content[2],40,50), substring(content[3],21,30),NA)
  NGmergeReport<- rbind(info,NGmergeReport)
}
colnames(NGmergeReport)<- c("sample","readPairs","adaptersRemoved","ratioRemoved")
NGmergeReport$readPairs <- as.numeric(NGmergeReport$readPairs)
NGmergeReport$adaptersRemoved <- as.numeric(NGmergeReport$adaptersRemoved)
NGmergeReport$ratioRemoved = NGmergeReport$adaptersRemoved / NGmergeReport$readPairs
NGmergeReport <- NGmergeReport[-c(1,2),]
median(NGmergeReport$readPairs)
jpeg(file="../readPairs.jpeg")
par(mar=c(5,11,1,4))
ratioRemovedPlot<- barplot(NGmergeReport$readPairs, horiz = T, xlab = "Total Fragment Pairs",
        names = NGmergeReport$sample,las=1, col="dark green", log = "x")
dev.off()
