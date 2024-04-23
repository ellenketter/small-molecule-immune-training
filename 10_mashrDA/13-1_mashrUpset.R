library(UpSetR)
library(dplyr)
library(ggplot2)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "upsetPlots/"


# Call results from text files ---------------------------------------------

mashSigResults <- readRDS("mashSigResults_dd.rds") 

listInput <- list(BetaGlucan = names(mashSigResults$bg),
                  FiveF = names(mashSigResults$'5f'),
                  Fenoterol = names(mashSigResults$fen),
                  Flunisolide = names(mashSigResults$fluni),
                  Hydrocortisone = names(mashSigResults$hc),
                  Hydroquinone = names(mashSigResults$hq),
                  Myricetin = names(mashSigResults$myr),
                  Nerol = names(mashSigResults$nerol))

tiff("upsetPlot_sigMash_top100.tiff", units="in", width=20, height=10, res=500)
upset(fromList(listInput),order.by = "freq", nsets = 8, nintersects = 100, 
      mainbar.y.label = "Shared Mash Peaks", sets.x.label = "Mash Peak Count")
dev.off()

# call only upregulated results from mash
mashMeans <- readRDS("mashMeans_dd.rds") 

treatments <- unique(mashMeans$treatment)
getPeaks <- function(x){mashMeans[mashMeans$treatment == x,]}
mashPeaks <- lapply(treatments, getPeaks)
names(mashPeaks) <- treatments

upregulatedMash <- mashMeans[mashMeans$significant == "UP",] # mean is positive, lfsr <0.05

treatments <- c(unique(upregulatedMash$treatment))
getUpPeaks <- function(x){upregulatedMash[upregulatedMash$treatment == x,]}
mashUpPeaks <- lapply(treatments, getUpPeaks)
new <- lapply(mashUpPeaks,"[",1)
names(new) <- treatments



listInput2 <- list(BetaGlucan = names(unlist(new$bg)),
                  FiveF = names(unlist(new$'5f')),
                  Fenoterol = names(unlist(new$fen)),
                  Flunisolide = names(unlist(new$fluni)),
                  Hydrocortisone = names(unlist(new$hc)),
                  Hydroquinone = names(unlist(new$hc)),
                  Myricetin = names(unlist(new$myr)),
                  Nerol = names(unlist(new$nerol)))

tiff("upsetPlot_sigMash_upregulated.tiff", units="in", width=10, height=10, res=500)
upset(fromList(listInput2),order.by = "freq", nsets=8, nintersects = 100,
      mainbar.y.label = "Shared Mash Peaks", sets.x.label = "Mash Peak Count")
dev.off()

# Plot Mash Peaks
ggplot(mashPeaks$myr, aes(x=posteriorMean, y=-log10(lfsr))) + 
  geom_point()

# List peaks which are DA vs PBS across all conditions
allIntersect <- which(rowSums(fromList(listInput)) == 8)
allPeakList<- c(mashSigResults$`5f`, mashSigResults$bg, mashSigResults$fen, mashSigResults$fluni,
  mashSigResults$hc, mashSigResults$hq, mashSigResults$myr, mashSigResults$nerol)
allPeakList <- allPeakList[!duplicated(allPeakList)]
allIntersect <- names(allPeakList)[allIntersect]

saveRDS(allIntersect, "allTreatmentIntersectPeaks.rds")


# Make a big table of sharing

# rownames = topGenePeak, colnames = stims, countif sig = UP

mash_res_list <- readRDS("mash_res_list.rds")

bg <- data.frame(mash_res_list$bg$Gene.Name, mash_res_list$bg$significant)
names(bg) <- c("gene", "bg")
`5f` <- data.frame(mash_res_list$`5f`$Gene.Name, mash_res_list$`5f`$significant)
names(`5f`) <- c("gene", "5f")
fen <- data.frame(mash_res_list$fen$Gene.Name, mash_res_list$fen$significant)
names(fen) <- c("gene", "fen")
fluni <- data.frame(mash_res_list$fluni$Gene.Name, mash_res_list$fluni$significant)
names(fluni) <- c("gene", "fluni")
hc <- data.frame(mash_res_list$hc$Gene.Name, mash_res_list$hc$significant)
names(hc) <- c("gene", "hc")
hq <- data.frame(mash_res_list$hq$Gene.Name, mash_res_list$hq$significant)
names(hq) <- c("gene", "hq")
myr <- data.frame(mash_res_list$myr$Gene.Name, mash_res_list$myr$significant)
names(myr) <- c("gene", "myr")
nerol <- data.frame(mash_res_list$nerol$Gene.Name, mash_res_list$nerol$significant)
names(nerol) <- c("gene", "nerol")

library(tidyverse)
list_df = list(bg, `5f`, fen, fluni, hc, hq, myr, nerol)
df2 <- list_df %>% reduce(inner_join)

df2[df2 == "NO"] <- 0
df2[df2 == "DOWN"] <- 0
df2[df2 == "UP"] <- 1

rownames(df2) <- df2$gene
df2 <- df2[-1]

# vectors of interest
# everything shared
all <-  c(1,1,1,1, 1,1,1,1)
geneNames_all <- names(which(colSums(t(df2) == all) == ncol(df2)))
geneNames_all <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_all,]
write.csv(geneNames_all, "./upsetPlots/geneNames_all.csv")
# everything but beta glucan
noBG <-  c(0,1,1,1, 1,1,1,1)
geneNames_noBg <- names(which(colSums(t(df2) == noBG) == ncol(df2)))
geneNames_noBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noBg,]
write.csv(geneNames_noBg, "./upsetPlots/geneNames_noBg.csv")
# hq-fluni-bg
hcFluniBG <-  c(1,0,0,1, 1,0,0,0)
geneNames_hcFluniBg <- names(which(colSums(t(df2) == hcFluniBG) == ncol(df2)))
geneNames_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluniBg,]
write.csv(geneNames_hcFluniBg, "./upsetPlots/geneNames_hcFluniBg.csv")
# hc-fluni only 
hcFluni <-  c(0,0,0,1, 1,0,0,0)
geneNames_hcFluni <- names(which(colSums(t(df2) == hcFluni) == ncol(df2)))
geneNames_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluni,]
write.csv(geneNames_hcFluni, "./upsetPlots/geneNames_hcFluni.csv")
# not hc-fluni
no_hcFluni <-  c(1,1,1,0, 0,1,1,1)
geneNames_no_hcFluni <- names(which(colSums(t(df2) == no_hcFluni) == ncol(df2)))
geneNames_no_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluni,]
write.csv(geneNames_no_hcFluni, "./upsetPlots/geneNames_no_hcFluni.csv")
# myr-nerol-hq-5f-fen only
no_hcFluniBg <-  c(0,1,1,0, 0,1,1,1)
geneNames_no_hcFluniBg <- names(which(colSums(t(df2) == no_hcFluniBg) == ncol(df2)))
geneNames_no_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluniBg,]
write.csv(geneNames_no_hcFluniBg, "./upsetPlots/geneNames_no_hcFluniBg.csv")





which(rownames(df2) == "Slc2a9")
df2[10991,]
