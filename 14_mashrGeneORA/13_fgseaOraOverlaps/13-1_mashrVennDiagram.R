library(UpSetR)
library(dplyr)
library(ggplot2)
library(VennDiagram)
library(RColorBrewer)
library(tidyverse)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "upsetPlots/vennDiagram/"
# dir.create(OUT_DIR)

# Make a big table of sharing

# rownames = topGenePeak, colnames = stims, countif sig = UP

mash_res_list <- readRDS("mash_res_list.rds")


# beta glucan
bg <- data.frame(mash_res_list$bg$Gene.Name, mash_res_list$bg$significant)
names(bg) <- c("gene", "bg")

# steroids
fluni <- data.frame(mash_res_list$fluni$Gene.Name, mash_res_list$fluni$significant)
names(fluni) <- c("gene", "fluni")
hc <- data.frame(mash_res_list$hc$Gene.Name, mash_res_list$hc$significant)
names(hc) <- c("gene", "hc")

list_steroids = list(fluni, hc)
steroids <- list_steroids %>% reduce(inner_join)

# other drugs
`5f` <- data.frame(mash_res_list$`5f`$Gene.Name, mash_res_list$`5f`$significant)
names(`5f`) <- c("gene", "5f")
fen <- data.frame(mash_res_list$fen$Gene.Name, mash_res_list$fen$significant)
names(fen) <- c("gene", "fen")
hq <- data.frame(mash_res_list$hq$Gene.Name, mash_res_list$hq$significant)
names(hq) <- c("gene", "hq")
myr <- data.frame(mash_res_list$myr$Gene.Name, mash_res_list$myr$significant)
names(myr) <- c("gene", "myr")
nerol <- data.frame(mash_res_list$nerol$Gene.Name, mash_res_list$nerol$significant)
names(nerol) <- c("gene", "nerol")

list_otherDrugs = list(`5f`, fen, hq, myr, nerol)
otherDrugs <- list_otherDrugs %>% reduce(inner_join)

# find upregulated, shared genes
bgUp <- bg$gene[bg$bg == "UP"]
steroidsUp <- steroids$gene[steroids$fluni == "UP" & steroids$hc == "UP"]



otherDrugsUp <- otherDrugs$gene[otherDrugs$`5f` == "UP" & otherDrugs$fen == "UP" & 
                                  otherDrugs$hq == "UP" & otherDrugs$myr == "UP" & 
                                  otherDrugs$nerol == "UP"]

listAll <- list(bgUp = bgUp,
                steroidsUp = steroidsUp,
                otherDrugsUp = otherDrugsUp)

png("upsetPlot_sigMash_venn_up.png", units="in", width=6, height=5, res=500)
upset(fromList(listAll),order.by = "freq", nsets = 8, nintersects = 12, 
      mainbar.y.label = "Shared Mash Peaks: Upregulated", sets.x.label = "Mash Peak Count")
dev.off()

# Venn Diagram visualization
myCol <- brewer.pal(3, "Pastel2")

# Chart
venn.diagram(
  x = list(bgUp, steroidsUp, otherDrugsUp),
  category.names = c("BG" , "Steroids " , "OtherDrugs"),
  filename = 'venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .6,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.6,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)





dfAll <- list(bgUp = data.frame(gene = bgUp, bg = rep(1)),
                steroidsUp = data.frame(gene = steroidsUp,steroids = rep(1)),
                otherDrugsUp = data.frame(gene = otherDrugsUp,otherDrugs = rep(1)))



all <- dfAll %>% reduce(full_join)
rownames(all) <- all$gene
all <- all[-1]

all[is.na(all)] <- 0




# vectors of interest

# Beta Glucan Alone

bgOnly_1 <-  c(1,0,0)
geneNames_bgOnly_1 <- names(which(colSums(t(all) == bgOnly_1) == ncol(all)))
geneNames_bgOnly_1 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_bgOnly_1,]
write.csv(geneNames_bgOnly_1, "./upsetPlots/vennDiagram/geneNames_bgOnly_1.csv")

bgOtherDrugs_2 <- c(1,0,1)
geneNames_bgOtherDrugs_2 <- names(which(colSums(t(all) == bgOtherDrugs_2) == ncol(all)))
geneNames_bgOtherDrugs_2 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_bgOtherDrugs_2,]
write.csv(geneNames_bgOtherDrugs_2, "./upsetPlots/vennDiagram/geneNames_bgOtherDrugs_2.csv")


bgSteroids_3 = c(1,1,0)
geneNames_bgSteroids_3 <- names(which(colSums(t(all) == bgSteroids_3) == ncol(all)))
geneNames_bgSteroids_3 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_bgSteroids_3,]
write.csv(geneNames_bgSteroids_3, "./upsetPlots/vennDiagram/geneNames_bgSteroids_3.csv")

steroidsOnly_4 = c(0,1,0)
geneNames_steroidsOnly_4 <- names(which(colSums(t(all) == steroidsOnly_4) == ncol(all)))
geneNames_steroidsOnly_4 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_steroidsOnly_4,]
write.csv(geneNames_steroidsOnly_4, "./upsetPlots/vennDiagram/geneNames_steroidsOnly_4.csv")

all_5 = c(1,1,1)
geneNames_all_5 <- names(which(colSums(t(all) == all_5) == ncol(all)))
geneNames_all_5 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_all_5,]
write.csv(geneNames_all_5, "./upsetPlots/vennDiagram/geneNames_all_5.csv")


otherDrugsOnly_6 = c(0,0,1)
geneNames_otherDrugsOnly_6 <- names(which(colSums(t(all) == otherDrugsOnly_6) == ncol(all)))
geneNames_otherDrugsOnly_6 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_otherDrugsOnly_6,]
write.csv(geneNames_otherDrugsOnly_6, "./upsetPlots/vennDiagram/geneNames_otherDrugsOnly_6.csv")

steroidsOtherDrugs_7 = c(0,1,1)
geneNames_steroidsOtherDrugs_7 <- names(which(colSums(t(all) == steroidsOtherDrugs_7) == ncol(all)))
geneNames_steroidsOtherDrugs_7 <- mash_res_list$bg[mash_res_list$bg$Gene.Name %in% geneNames_steroidsOtherDrugs_7,]
write.csv(geneNames_steroidsOtherDrugs_7, "./upsetPlots/vennDiagram/geneNames_steroidsOtherDrugs_7.csv")



