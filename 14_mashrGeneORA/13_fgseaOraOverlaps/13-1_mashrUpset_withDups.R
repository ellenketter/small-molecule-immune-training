library(UpSetR)
library(dplyr)
library(ggplot2)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "upsetPlots/withDups/"
# dir.create(OUT_DIR)

# Make a big table of sharing

# rownames = topGenePeak, colnames = stims, countif sig = UP

mash_res_list <- readRDS("mash_res_list.rds")

bg <- data.frame(mash_res_list$bg$Gene.Name, mash_res_list$bg$significant,
                 mash_res_list$bg$peak)
names(bg) <- c("gene", "bg","peak")
`5f` <- data.frame(mash_res_list$`5f`$Gene.Name, mash_res_list$`5f`$significant,
                   mash_res_list$`5f`$peak)
names(`5f`) <- c("gene", "5f", "peak")
fen <- data.frame(mash_res_list$fen$Gene.Name, mash_res_list$fen$significant,
                  mash_res_list$fen$peak)
names(fen) <- c("gene", "fen", "peak")
fluni <- data.frame(mash_res_list$fluni$Gene.Name, mash_res_list$fluni$significant,
                    mash_res_list$fluni$peak)
names(fluni) <- c("gene", "fluni", "peak")
hc <- data.frame(mash_res_list$hc$Gene.Name, mash_res_list$hc$significant,
                 mash_res_list$hc$peak)
names(hc) <- c("gene", "hc","peak")
hq <- data.frame(mash_res_list$hq$Gene.Name, mash_res_list$hq$significant,
                 mash_res_list$hq$peak)
names(hq) <- c("gene", "hq","peak")
myr <- data.frame(mash_res_list$myr$Gene.Name, mash_res_list$myr$significant,
                  mash_res_list$myr$peak)
names(myr) <- c("gene", "myr","peak")
nerol <- data.frame(mash_res_list$nerol$Gene.Name, mash_res_list$nerol$significant,
                    mash_res_list$nerol$peak)
names(nerol) <- c("gene", "nerol","peak")

library(tidyverse)
list_df = list(bg, `5f`, fen, fluni, hc, hq, myr, nerol)
# names(list_df) <- c("bg", "5f", "fen", "fluni", "hc", "hq", "myr", "nerol")
df2 <- list_df %>% reduce(inner_join)


listInput_up <- list(BetaGlucan = df2$peak[df2$bg == "UP"],
                  FiveF = df2$peak[df2$`5f` == "UP"],
                  Fenoterol = df2$peak[df2$fen == "UP"],
                  Flunisolide = df2$peak[df2$fluni == "UP"],
                  Hydrocortisone = df2$peak[df2$hc == "UP"],
                  Hydroquinone = df2$peak[df2$hq == "UP"],
                  Myricetin = df2$peak[df2$myr == "UP"],
                  Nerol = df2$peak[df2$nerol == "UP"])

png("upsetPlot_sigMash_top12_upPeaks_withDups.png", units="in", width=6, height=5, res=500)
upset(fromList(listInput_up),order.by = "freq", nsets = 8, nintersects = 12, 
      mainbar.y.label = "Shared Mash Peaks: Upregulated", sets.x.label = "Mash Peak Count")
dev.off()

dfs_up <- list(BetaGlucan = data.frame(gene = df2$gene[df2$bg == "UP"], bg = 1),
                     FiveF = data.frame(gene =df2$gene[df2$`5f` == "UP"], `5f`= 1),
                     Fenoterol = data.frame(gene =df2$gene[df2$fen == "UP"],fen= 1),
                     Flunisolide = data.frame(gene =df2$gene[df2$fluni == "UP"],fluni=1),
                     Hydrocortisone = data.frame(gene =df2$gene[df2$hc == "UP"],hc=1),
                     Hydroquinone = data.frame(gene =df2$gene[df2$hq == "UP"],hq=1),
                     Myricetin = data.frame(gene =df2$gene[df2$myr == "UP"],myr=1),
                     Nerol = data.frame(gene =df2$gene[df2$nerol == "UP"],nerol=1))

combined <- dfs_up %>% reduce(full_join)

# remove gene with empty name
combined <- combined[!combined$gene == "",]

# replace NAs with 0s
combined[is.na(combined)] <- 0



rownames(combined) <- combined$gene
# remove gene column
combined <- combined[-1]

# load mash_res_list without duplicated genes
mash_res_list <- readRDS("mash_res_list.rds")

# vectors of interest
allBgPeaksUp <- mash_res_list$bg[mash_res_list$bg$significant == "UP",]
write.csv(allBgPeaksUp, paste0(OUT_DIR, "allBgPeaksUp.csv"))

# Beta Glucan Alone
bgAlone_1 <-  c(1,0,0,0, 0,0,0,0)
geneNames_bgAlone_1 <- names(which(colSums(t(combined) == bgAlone_1) == ncol(combined)))
geneNames_bgAlone_1 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_bgAlone_1,]
write.csv(geneNames_bgAlone_1, paste0(OUT_DIR, "geneNames_bgAlone_1.csv"))

# not hc-fluni
no_hcFluni <-  c(1,1,1,0, 0,1,1,1)
geneNames_no_hcFluni <- names(which(colSums(t(combined) == no_hcFluni) == ncol(combined)))
geneNames_no_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluni,]
write.csv(geneNames_no_hcFluni, paste0(OUT_DIR, "geneNames_no_hcFluni_2.csv"))

# myr + bg
myrBg_3 <- c(1,0,0,0, 0,0,1,0)
geneNames3 <- names(which(colSums(t(combined) == myrBg_3) == ncol(combined)))
geneNames3 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames3,]
write.csv(geneNames3, paste0(OUT_DIR, "geneNames_myrBg_3.csv"))

# everything shared
all <-  c(1,1,1,1, 1,1,1,1)
geneNames_all <- names(which(colSums(t(combined) == all) == ncol(combined)))
geneNames_all <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_all,]
write.csv(geneNames_all, paste0(OUT_DIR, "geneNames_all_4.csv"))

# hc-fluni only 
hcFluni <-  c(0,0,0,1, 1,0,0,0)
geneNames_hcFluni <- names(which(colSums(t(combined) == hcFluni) == ncol(combined)))
geneNames_hcFluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluni,]
write.csv(geneNames_hcFluni, paste0(OUT_DIR, "geneNames_hcFluni_5.csv"))

# hc-fluni-bg
hcFluniBG <-  c(1,0,0,1, 1,0,0,0)
geneNames_hcFluniBg <- names(which(colSums(t(combined) == hcFluniBG) == ncol(combined)))
geneNames_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcFluniBg,]
write.csv(geneNames_hcFluniBg, paste0(OUT_DIR, "geneNames_hcFluniBg_6.csv"))

# fluni-bg
fluniBG <-  c(1,0,0,1, 0,0,0,0)
geneNames_fluniBg <- names(which(colSums(t(combined) == fluniBG) == ncol(combined)))
geneNames_fluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fluniBg,]
write.csv(geneNames_fluniBg, paste0(OUT_DIR, "geneNames_fluniBg_7.csv"))


# myr-nerol-hq-5f-fen only
no_hcFluniBg <-  c(0,1,1,0, 0,1,1,1)
geneNames_no_hcFluniBg <- names(which(colSums(t(combined) == no_hcFluniBg) == ncol(combined)))
geneNames_no_hcFluniBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_no_hcFluniBg,]
write.csv(geneNames_no_hcFluniBg, paste0(OUT_DIR, "geneNames_no_hcFluniBg_8.csv"))

# everything but hc
noHC <-  c(1,1,1,1, 0,1,1,1)
geneNames_noHC <- names(which(colSums(t(combined) == noHC) == ncol(combined)))
geneNames_noHC <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noHC,]
write.csv(geneNames_noHC, paste0(OUT_DIR, "geneNames_noHC_9.csv"))

# only fluni
fluni <-  c(0,0,0,1, 0,0,0,0)
geneNames_fluni <- names(which(colSums(t(combined) == fluni) == ncol(combined)))
geneNames_fluni <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fluni,]
write.csv(geneNames_fluni, paste0(OUT_DIR, "geneNames_fluniOnly_10.csv"))

# only myr
myr <-  c(0,0,0,0, 0,0,1,0)
geneNames_myr <- names(which(colSums(t(combined) == myr) == ncol(combined)))
geneNames_myr <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_myr,]
write.csv(geneNames_myr, paste0(OUT_DIR, "geneNames_myrOnly_11.csv"))

# only HC
hc <-  c(0,0,0,0, 1,0,0,0)
geneNames_hc <- names(which(colSums(t(combined) == hc) == ncol(combined)))
geneNames_hc <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hc,]
write.csv(geneNames_hc, paste0(OUT_DIR, "geneNames_hcOnly_12.csv"))

# everything but beta glucan
noBG <-  c(0,1,1,1, 1,1,1,1)
geneNames_noBg <- names(which(colSums(t(combined) == noBG) == ncol(combined)))
geneNames_noBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noBg,]
write.csv(geneNames_noBg, paste0(OUT_DIR, "geneNames_noBg_13.csv"))

## Additional from dd + commonbaseline------------
# NoHcFluni5f

noHcFluni5f <- c(1,0,1,0, 0,1,1,1)
geneNames_noHcFluni5f <- names(which(colSums(t(combined) == noHcFluni5f) == ncol(combined)))
geneNames_noHcFluni5f <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_noHcFluni5f,]
write.csv(geneNames_noHcFluni5f, paste0(OUT_DIR, "geneNames_noHcFluni5f_14.csv"))

fenBg <- c(1,0,1,0, 0,0,0,0)
geneNames_fenBg <- names(which(colSums(t(combined) == fenBg) == ncol(combined)))
geneNames_fenBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_fenBg,]
write.csv(geneNames_fenBg, paste0(OUT_DIR, "geneNames_fenBg_15.csv"))

hcBg <- c(1,0,0,0, 1,0,0,0)
geneNames_hcBg <- names(which(colSums(t(combined) == hcBg) == ncol(combined)))
geneNames_hcBg <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_hcBg,]
write.csv(geneNames_hcBg, paste0(OUT_DIR, "geneNames_hcBg_16.csv"))


which(rownames(combined) == "Slc2a9")
combined[9441,]

# library(partitions)
# myComps <- t(as.matrix(compositions(6, 8)))
# head(myComps)
# myComps2<- myComps[!rowSums(myComps>1),]

