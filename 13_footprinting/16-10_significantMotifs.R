# BiocManager::install("seqLogo")
###################################################
### code chunk number 1: ImportLibrary
###################################################
library(DiffLogo)


###################################################
### code chunk number 2: ImportMotifsFromMotifDb
###################################################
library(MotifDb)

## import motifs
fpStats <- readRDS("analysis/mergedStatsPadj.rds")
sigMotifs <- unique(fpStats$Motif[fpStats$pAdj <= 0.1])

originalHitIndex <- c()
allJaspar <- list()
for (i in 1:length(sigMotifs)){
  hitIndeces = grep (sigMotifs[i], values (MotifDb)$geneSymbol, ignore.case=TRUE)
  list= as.list(MotifDb[hitIndeces])
  chooseJaspar <- grep("Hsapiens-jaspar2022", names(list))
  originalHitIndex <- append(originalHitIndex, hitIndeces[chooseJaspar])
  jaspar <- as.list(list[chooseJaspar])
  allJaspar <- append(allJaspar, jaspar)
  print(sigMotifs[i])
  print(names(jaspar))}


# MA0838.1.CEBPG was the significant one, so remove MA1636.1
# sigMotifs Atf1 and NRF1 not in MotifDb

sigMotifs <- sigMotifs[-c(12,18)]
originalHitIndex <- originalHitIndex[-30]
allJaspar <- allJaspar[-30]

# sequence counts 
sequenceCounts = as.numeric(values (MotifDb)$sequenceCount[originalHitIndex])
names(sequenceCounts) = names(allJaspar)

# Manually add info from jaspar for Nrf1 (2020 https://jaspar2020.genereg.net/matrix/MA0506.1/)
# no human jaspar motif for Atf1
allJaspar$NRF1 <- getPwmFromPfmOrJasparFile("analysis/MA0506.1.jaspar")
allJaspar$NRF1 <- allJaspar$NRF1[,-1]
colnames(allJaspar$NRF1) <- seq(1:11) 
rownames(allJaspar$NRF1) <- rownames(allJaspar$KLF6)
sigMotifs <- append(sigMotifs, "NRF1")
names(allJaspar) <- sigMotifs

###################################################
### code chunk number 3: PlotSequenceLogo
###################################################
## plot classic sequence logos
# install.packages("ggseqlogo")
library(ggplot2)
library(ggseqlogo)
for(i in 1:length(sigMotifs)){
  p <- ggseqlogo(allJaspar[[i]]) + ylab('Information Content') +
    xlab('Nucleotide Position') + 
    ggtitle(paste0(sigMotifs[i]," - JASPAR Motif Sequence (2022)")) + 
    theme_classic()
  ggsave(paste0(sigMotifs[[i]],".jpg"), p, width = 4, height = 3)}


###################################################
### code chunk number 4: PlotDiffLogoTable
###################################################
## plot table of difference logos for CTFC motifs (DNA)

pdf('DiffLogoTable.pdf', height = 75, width = 55) 
par(mar= c(1,1,1,1))
diffLogoTable(PWMs = allJaspar)
dev.off()


# export table of significant motifs

library(reshape2)
library(dplyr)
library(openxlsx)
sigMotifs <- unique(fpStats$Motif[fpStats$pAdj <= 0.1])
dfSigMotifs <- mergedStatsAdjP %>% filter(Motif %in% sigMotifs)
dfSigMotifs_pAdj <- dcast(dfSigMotifs, Motif ~ molecule, value.var = "pAdj")
dfSigMotifs_foldChange <- dcast(dfSigMotifs, Motif ~ molecule, value.var = "foldChange")
write.xlsx(dfSigMotifs_pAdj, "adjustedPvalues_significantMotifs.xlsx")
write.xlsx(dfSigMotifs_foldChange, "foldChangeVsPBS_significantMotifs.xlsx")
