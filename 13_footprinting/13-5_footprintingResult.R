library(ggplot2)
library(ggrepel)
library(dplyr)
library(tidyr)
library(ComplexHeatmap)
library(UpSetR)
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(cowplot))

setwd("~/Downloads/")
list <- readRDS("mergedStats.rds") # TF Activity is NOT fold change, it's difference between PBS - test
columns<- colnames(list$BG)
columns[c(3,5)]<-c("Protection_Score","TC") 
for (i in 1:length(list)){
  colnames(list[[i]]) <- columns
}


# make dataframe
new <- bind_rows(list, .id = "molecule")
motif <- str_split_fixed(new$Motif, '\\.', 3)
new$Motif <- motif[,3]

# calculate TF activity scores
new$testScore <- new$Protection_Score + new$TC
new$pbsScore <- new$Protection_Score_PBS + new$TC_PBS
new$foldChange<- new$testScore/new$pbsScore  # = foldchange

# Make a dataframe of molecule, fold change, p value, and motif
molecule <- unique(new$molecule)
dfLong <- subset(new, select = c("molecule","Motif", "P_values", "foldChange"))
dfList <- split(dfLong, dfLong$molecule)


dfList$`5F`$pAdj <- p.adjust(list$`5F`$P_values, method = "BH")
dfList$BG$pAdj <- p.adjust(list$BG$P_values, method = "BH")
dfList$Fen$pAdj <- p.adjust(list$Fen$P_values, method = "BH")
dfList$HC$pAdj <- p.adjust(list$HC$P_values, method = "BH")
dfList$HQ$pAdj <- p.adjust(list$HQ$P_values, method = "BH")
dfList$Fluni$pAdj <- p.adjust(list$Fluni$P_values, method = "BH")
dfList$Myr$pAdj <- p.adjust(list$Myr$P_values, method = "BH")
dfList$Nerol$pAdj <- p.adjust(list$Nerol$P_values, method = "BH")

figureDF <- bind_rows(dfList, .id = "molecule")
figureDF$Num <- new$Num
figureDF <- subset(figureDF, Num > 1000) # filter out motifs with lower sample rate

# determine significance based on adjusted p values
figureDF$significance <- "NotSignificant"
figureDF$significance[log2(figureDF$foldChange) > 0 & figureDF$pAdj < 0.05] <- "IncreasedVsPBS"
figureDF$significance[log2(figureDF$foldChange) < 0 & figureDF$pAdj < 0.05] <- "DecreasedVsPBS"

saveRDS(figureDF, "mergedStatsPadj.rds")

# heatmap by fold change ------------

df <- subset(figureDF, select = c(Motif,molecule,foldChange))
df <- df %>%  spread(key = molecule, value = foldChange)
rownames(df) <- df$Motif
df <- df[-1]


df$Var <- apply(df, 1, sd)

df <- df %>%
  top_n(20, wt = Var) %>%
  select(-Var)

df <- t(scale(t(df)))

options(repr.plot.width = 3, repr.plot.height = 12)
pdf("heatmapTop20.pdf", width = 6, height = 8)
p <- Heatmap(as.matrix(df),
             name = "scaledFoldChangeVsPBS",
             cluster_columns = TRUE,
             cluster_rows = TRUE,
             show_row_names = TRUE,
             rect_gp = gpar(col = "black", lwd = 0.5)
)
p
dev.off()
# heatmap by adjusted p value ------------
# heatmap by fold change ------------

df <- subset(figureDF, select = c(Motif,molecule,pAdj))
df <- df %>%  spread(key = molecule, value = pAdj)
rownames(df) <- df$Motif
df <- df[-1]


df$Var <- apply(df, 1, sd)

df <- df %>%
  top_n(60, wt = Var) %>%
  select(-Var)

# df <- t(scale(t(df)))

options(repr.plot.width = 3, repr.plot.height = 12)
pdf("heatmapTop60_pAdj.pdf", width = 5, height = 12)
p <- Heatmap(as.matrix(df),
             name = "adjustedPvalue",
             cluster_columns = TRUE,
             cluster_rows = TRUE,
             show_row_names = TRUE,
             rect_gp = gpar(col = "black", lwd = 0.5)
)
p
dev.off()



## Summary of significant results ---------
sigList<- split(figureDF, figureDF$molecule)

BetaGlucan <- c(sum(sigList$BG$significance == "IncreasedVsPBS"), -sum(sigList$BG$significance == "DecreasedVsPBS"))
FluoroindoleCarboxylicAcid <- c(sum(sigList$'5F'$significance == "IncreasedVsPBS"), -sum(sigList$'5F'$significance == "DecreasedVsPBS"))
Fenoterol <- c(sum(sigList$Fen$significance == "IncreasedVsPBS"), -sum(sigList$Fen$significance == "DecreasedVsPBS"))
Flunisolide <- c(sum(sigList$Fluni$significance == "IncreasedVsPBS"), -sum(sigList$Fluni$significance == "DecreasedVsPBS"))
Hydrocortisone <- c(sum(sigList$HC$significance == "IncreasedVsPBS"), -sum(sigList$HC$significance == "DecreasedVsPBS"))
Hydroquinone <- c(sum(sigList$HQ$significance == "IncreasedVsPBS"), -sum(sigList$HQ$significance == "DecreasedVsPBS"))
Myricetin <- c(sum(sigList$Myr$significance == "IncreasedVsPBS"), -sum(sigList$Myr$significance == "DecreasedVsPBS"))
Nerol <- c(sum(sigList$Nerol$significance == "IncreasedVsPBS"), -sum(sigList$Nerol$significance == "DecreasedVsPBS"))

countsMash <- data.frame(Flunisolide,  Hydrocortisone, Myricetin, Fenoterol, 
                         Nerol, FluoroindoleCarboxylicAcid, Hydroquinone,  BetaGlucan)
rownames(countsMash) <- c("Upregulated", "Downregulated")
countsMash <- reshape2::melt(as.matrix(countsMash)) 
colnames(countsMash) <- c("Direction","Treatment","SignificantFootprints")

svg(file= "summaryFootprintsDR_0.05.svg", 
    width = 5, height = 4)
ggplot(countsMash, aes(fill=Treatment, y=SignificantFootprints, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  #geom_text(aes(label = SignificantFootprints), position = position_stack(), size = 3, colour = "black") + 
  scale_fill_manual(values = c("#2BC4E4", "#F9987A", "#485695","#7C543E", 
                               "thistle4", "#11A290","#7BD1C1", "#B3B7B9") ) +
  geom_hline(aes(yintercept = 0)) + 
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust =1, vjust=1)) +
  ylab("Number of Transcription Factors 
Differentially Regulated vs PBS

Downregulated      Upregulated") +
  labs(subtitle = '    Steroids                      Unrelated              
                                   Small Molecules       ') +
  xlab("") +
  ylim(-10,30) 

dev.off()


## TF Families DR in stimulated vs PBS ------
makeSigLists<- function(x) {sigList%>% list(toupper(x$Motif[x$significance == "IncreasedVsPBS"]), 
                   toupper(x$Motif[x$significance == "DecreasedVsPBS"]))}
sigListPerMolecule<- lapply(molecule, makeSigLists)

names(BetaGlucan) <- c("Upregulated", "Downregulated")
# BG_up <- unique(gsub("\\d", "",BetaGlucan$Upregulated))
# for (i in 1:length(BetaGlucan$Upregulated)) {
#   BG_up[i] <- sum(str_count(BetaGlucan$Upregulated, (gsub("\\d", "",BetaGlucan$Upregulated))[i]))
# } # try to group by familes

FluoroindoleCarboxylicAcid <- data.frame(sigList$'5F'$significance == "IncreasedVsPBS", sigList$'5F'$significance == "DecreasedVsPBS")
Fenoterol <- data.frame(sigList$Fen$significance == "IncreasedVsPBS", sigList$Fen$significance == "DecreasedVsPBS")
Flunisolide <- data.frame(sigList$Fluni$significance == "IncreasedVsPBS", sigList$Fluni$significance == "DecreasedVsPBS")
Hydrocortisone <- data.frame(sigList$HC$significance == "IncreasedVsPBS", sigList$HC$significance == "DecreasedVsPBS")
Hydroquinone <- data.frame(sigList$HQ$significance == "IncreasedVsPBS", sigList$HQ$significance == "DecreasedVsPBS")
Myricetin <- data.frame(sigList$Myr$significance == "IncreasedVsPBS", sigList$Myr$significance == "DecreasedVsPBS")
Nerol <- data.frame(sigList$Nerol$significance == "IncreasedVsPBS", sigList$Nerol$significance == "DecreasedVsPBS")




## Venn Diagram ----------------------
# Make lists of significant changed footprints
# 
# For every molecule, list motif if significance is not "NotSignificant
molecule <- unique(new$molecule)

subsetUpregulated <- function(x) {
  new$Motif[new$molecule == x & new$significance == "IncreasedVsPBS"]
}
upFootprintFC<- lapply(molecule, subsetUpregulated)
names(significantFootprintFC) <- molecule
subsetDownregulated <- function(x) {
  new$Motif[new$molecule == x & new$significance == "DecreasedVsPBS"]
}
downFootprintFC<- lapply(molecule, subsetDownregulated)
names(downFootprintFC) <- molecule

library(UpSetR)
pdf("upsetPlot.pdf")
upset(fromList(significantFootprintFC),order.by = "freq", nsets = 8, nintersects = 11, 
      mainbar.y.label = "Shared Upregulated TFs Vs PBS (Top 11 Intersects)", 
      sets.x.label = "Frequency of Sharing")
dev.off()

unrelatedMolecules<-list(significantFootprintFC$`5F`, significantFootprintFC$Fen, 
                                significantFootprintFC$HQ, significantFootprintFC$Myr, 
                                significantFootprintFC$Nerol, significantFootprintFC$BG)
names(unrelatedMolecules) <- c("5F", "Fen", "HQ", "Myr", "Nerol", "BG")

intersectUnrelatedMolecules <- Reduce(intersect, unrelatedMolecules)

# Output signifcance of several key molecules
# make a heatmap of the following families
# SP, KLF, EGR, MAX, FOS, JUN, CEBP
# and all steroid markers 
# Up: "PHOX2B"  "ARGFX"   "NKX6-2"  "NKX6-1"  "Lhx3"    "MNX1"    "NKX2-8"  "ALX3"    "ONECUT1" "Arid3a"  "ZNF384" 
# Down: 

# volcano plots
df <- split(new, new$molecule)

pdf("volcanoBG.pdf")
BG <- ggplot(data=df$BG, aes(x=log2(foldChange), y=-log10(P_values), 
                             col=significance, label = Motif)) +
  ggtitle("Beta Glucan: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel(max.overlaps=30) +
  scale_color_manual(values=c("blue", "red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
BG
dev.off()

pdf("volcano5F.pdf")
`5F` <- ggplot(data=df$`5F`, aes(x=log2(foldChange), y=-log10(P_values),
                                 col=significance, label = Motif)) +
  ggtitle("5F: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
`5F`
dev.off()

pdf("volcanoFen.pdf")
Fen <- ggplot(data=df$Fen, aes(x=log2(foldChange), y=-log10(P_values), 
              col=significance, label = Motif)) +
  ggtitle("Fenoterol: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 
Fen
dev.off()

pdf("volcanoFluni.pdf")
Fluni <- ggplot(data=df$Fluni, aes(x=log2(foldChange), y=-log10(P_values), 
                                   col=significance, label = Motif)) +
  ggtitle("Flunisolide: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("blue", "red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
Fluni
dev.off()

pdf("volcanoHC.pdf")
HC <- ggplot(data=df$HC, aes(x=log2(foldChange), y=-log10(P_values),
                             col=significance, label = Motif)) +
  ggtitle("Hydrocortisone: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel(max.overlaps = 30) +
  scale_color_manual(values=c("blue", "red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
HC
dev.off()

pdf("volcanoHQ.pdf")
HQ <- ggplot(data=df$HQ, aes(x=log2(foldChange), y=-log10(P_values),
                             col=significance, label = Motif)) +
  ggtitle("Hydroquinone: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel(max.overlaps = 25) +
  scale_color_manual(values=c("blue", "red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red") 
HQ
dev.off()

pdf("volcanoMyr.pdf")
Myr <- ggplot(data=df$Myr, aes(x=log2(foldChange), y=-log10(P_values),
                               col=significance, label = Motif)) +
  ggtitle("Myricetin: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel(max.overlaps = 25) +
  scale_color_manual(values=c("blue", "red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
Myr
dev.off()

pdf("volcanoNerol.pdf")
Nerol <- ggplot(data=df$Nerol, aes(x=log2(foldChange), y=-log10(P_values),
                                   col=significance, label = Motif)) +
  ggtitle("Nerol: TF Fold Change Vs PBS by Footprinting") + 
  geom_point(size=0.5) +
  xlim(-0.5,0.5) +
  ylim(0,6) +
  theme_minimal() +
  geom_text_repel(max.overlaps = 20) +
  scale_color_manual(values=c("red",  "black")) +
  geom_vline(xintercept=c(-0.1, 0.1), col="red") +
  geom_hline(yintercept=-log10(0.05), col="red")
Nerol
dev.off()

