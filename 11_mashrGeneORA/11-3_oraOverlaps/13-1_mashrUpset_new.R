library(UpSetR)
library(dplyr)
library(ggplot2)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "upsetPlots/"
# dir.create(OUT_DIR)

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
# names(list_df) <- c("bg", "5f", "fen", "fluni", "hc", "hq", "myr", "nerol")
df2 <- list_df %>% reduce(inner_join)

listInput_up <- list(BetaGlucan = df2$gene[df2$bg == "UP"],
                  FiveF = df2$gene[df2$`5f` == "UP"],
                  Fenoterol = df2$gene[df2$fen == "UP"],
                  Flunisolide = df2$gene[df2$fluni == "UP"],
                  Hydrocortisone = df2$gene[df2$hc == "UP"],
                  Hydroquinone = df2$gene[df2$hq == "UP"],
                  Myricetin = df2$gene[df2$myr == "UP"],
                  Nerol = df2$gene[df2$nerol == "UP"])

png("upsetPlot_sigMash_top12_up.png", units="in", width=6, height=5, res=500)
upset(fromList(listInput_up),order.by = "freq", nsets = 8, nintersects = 12, 
      mainbar.y.label = "Shared Mash Peaks: Upregulated", sets.x.label = "Mash Peak Count")
dev.off()

allUpPeaks <- unique(unname(unlist(listInput_up)))


df2[df2 == "NO"] <- 0
df2[df2 == "DOWN"] <- 0
df2[df2 == "UP"] <- 1

rownames(df2) <- df2$gene
df2 <- df2[-1]

# vectors of interest

# all upregulated peaks from all treatment conditions
geneNames_allUpPeaks <- allUpPeaks
geneNames_allUpPeaks <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_allUpPeaks,]
write.csv(geneNames_allUpPeaks, "./upsetPlots/geneNames_allUpPeaks.csv")

# Beta Glucan Alone
bgAlone_1 <-  c(1,0,0,0, 0,0,0,0)
geneNames_bgAlone_1 <- names(which(colSums(t(df2) == bgAlone_1) == ncol(df2)))
geneNames_bgAlone_1 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_bgAlone_1,]
write.csv(geneNames_bgAlone_1, "./upsetPlots/geneNames_bgAlone_1.csv")

# not hc-fluni
no_hcFluni <-  c(1,1,1,0, 0,1,1,1)
geneNames_no_hcFluni <- names(which(colSums(t(df2) == no_hcFluni) == ncol(df2)))
geneNames_no_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluni,]
write.csv(geneNames_no_hcFluni, "./upsetPlots/geneNames_no_hcFluni_2.csv")

# myr + bg
myrBg_3 <- c(1,0,0,0, 0,0,1,0)
geneNames3 <- names(which(colSums(t(df2) == myrBg_3) == ncol(df2)))
geneNames3 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames3,]
write.csv(geneNames3, "./upsetPlots/geneNames_myrBg_3.csv")

# everything shared
all <-  c(1,1,1,1, 1,1,1,1)
geneNames_all <- names(which(colSums(t(df2) == all) == ncol(df2)))
geneNames_all <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_all,]
write.csv(geneNames_all, "./upsetPlots/geneNames_all_4.csv")

# hc-fluni only 
hcFluni <-  c(0,0,0,1, 1,0,0,0)
geneNames_hcFluni <- names(which(colSums(t(df2) == hcFluni) == ncol(df2)))
geneNames_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluni,]
write.csv(geneNames_hcFluni, "./upsetPlots/geneNames_hcFluni_5.csv")

# hc-fluni-bg
hcFluniBG <-  c(1,0,0,1, 1,0,0,0)
geneNames_hcFluniBg <- names(which(colSums(t(df2) == hcFluniBG) == ncol(df2)))
geneNames_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluniBg,]
write.csv(geneNames_hcFluniBg, "./upsetPlots/geneNames_hcFluniBg_6.csv")

# fluni-bg
fluniBG <-  c(1,0,0,1, 0,0,0,0)
geneNames_fluniBg <- names(which(colSums(t(df2) == fluniBG) == ncol(df2)))
geneNames_fluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fluniBg,]
write.csv(geneNames_fluniBg, "./upsetPlots/geneNames_fluniBg_7.csv")


# myr-nerol-hq-5f-fen only
no_hcFluniBg <-  c(0,1,1,0, 0,1,1,1)
geneNames_no_hcFluniBg <- names(which(colSums(t(df2) == no_hcFluniBg) == ncol(df2)))
geneNames_no_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluniBg,]
write.csv(geneNames_no_hcFluniBg, "./upsetPlots/geneNames_no_hcFluniBg_8.csv")

# everything but hc
noHC <-  c(1,1,1,1, 0,1,1,1)
geneNames_noHC <- names(which(colSums(t(df2) == noHC) == ncol(df2)))
geneNames_noHC <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noHC,]
write.csv(geneNames_noHC, "./upsetPlots/geneNames_noHC_9.csv")

# only fluni
fluni <-  c(0,0,0,1, 0,0,0,0)
geneNames_fluni <- names(which(colSums(t(df2) == fluni) == ncol(df2)))
geneNames_fluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fluni,]
write.csv(geneNames_fluni, "./upsetPlots/geneNames_fluniOnly_10.csv")

# only myr
myr <-  c(0,0,0,0, 0,0,1,0)
geneNames_myr <- names(which(colSums(t(df2) == myr) == ncol(df2)))
geneNames_myr <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_myr,]
write.csv(geneNames_myr, "./upsetPlots/geneNames_myrOnly_11.csv")

# only HC
hc <-  c(0,0,0,0, 1,0,0,0)
geneNames_hc <- names(which(colSums(t(df2) == hc) == ncol(df2)))
geneNames_hc <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hc,]
write.csv(geneNames_hc, "./upsetPlots/geneNames_hcOnly_12.csv")

# everything but beta glucan
noBG <-  c(0,1,1,1, 1,1,1,1)
geneNames_noBg <- names(which(colSums(t(df2) == noBG) == ncol(df2)))
geneNames_noBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noBg,]
write.csv(geneNames_noBg, "./upsetPlots/geneNames_noBg_13.csv")

## Additional from dd + commonbaseline------------
# NoHcFluni5f

noHcFluni5f <- c(1,0,1,0, 0,1,1,1)
geneNames_noHcFluni5f <- names(which(colSums(t(df2) == noHcFluni5f) == ncol(df2)))
geneNames_noHcFluni5f <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noHcFluni5f,]
write.csv(geneNames_noHcFluni5f, "./upsetPlots/geneNames_noHcFluni5f_14.csv")

fenBg <- c(1,0,1,0, 0,0,0,0)
geneNames_fenBg <- names(which(colSums(t(df2) == fenBg) == ncol(df2)))
geneNames_fenBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fenBg,]
write.csv(geneNames_fenBg, "./upsetPlots/geneNames_fenBg_15.csv")

hcBg <- c(1,0,0,0, 1,0,0,0)
geneNames_hcBg <- names(which(colSums(t(df2) == hcBg) == ncol(df2)))
geneNames_hcBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcBg,]
write.csv(geneNames_hcBg, "./upsetPlots/geneNames_hcBg_16.csv")


which(rownames(df2) == "Slc2a9")
df2[10991,]
