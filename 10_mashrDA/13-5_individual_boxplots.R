#R version 4.1.0 loaded

library(ggplot2)
library(DESeq2)
library(statmod)
library(RColorBrewer)
library(edgeR)
library(devtools)
library(cluster)
library(ggrepel)
library(dplyr)


setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/11_annotatePeaks/")

OUT_DIR <- "11_individual_boxplots/"
# dir.create(OUT_DIR)
expression <- read.csv("../10_differentialAccessibility/expressionCPMs_allPBS_processOrder.txt")
annotate <- read.delim("old/background_annotate_out.txt")
annotate_names <- colnames(annotate)
annotate_names[1] <- "PeakID"
colnames(annotate) <- annotate_names

annotate[annotate$Gene.Name == "Zeb2",]



# load all limma results files to associate gene names with peak IDs


sampleCpm <- list.files("../10_differentialAccessibility/background_vs_DA_peaks/" ,pattern = "\\_")

cpms <- c()  
readFiles <- function(x){
  cpms <- read.delim(paste0("../10_differentialAccessibility/background_vs_DA_peaks/", x))
}

gcUnique <- lapply (sampleCpm, readFiles)

finalNames <- gsub(".txt", "", sampleCpm)
finalNames <- gsub("results_","",finalNames)
names(gcUnique) <- finalNames

peaksScarb1 <- annotate[annotate$Gene.Name =="Scarb1",]
saveRDS(peaksScarb1, "../13_mashr/peaksScarb1.rds")
peaksSh3glb1 <- annotate[annotate$Gene.Name =="Sh3glb1",]
saveRDS(peaksSh3glb1, "../13_mashr/peaksSh3glb1.rds")
expressionScarb1 <- data.frame()
for (i in 1:length(peaksScarb1$PeakID)){
  expressionScarb1 <- rbind(expressionScarb1, expression[rownames(expression)== peaksScarb1$PeakID[i],])
}

# make into a list of averages
treatments <- unique(gsub(pattern = "\\.\\d", replacement = "", colnames(expressionScarb1)))

scarb1Cpms <- c()  
splitDF <- function(x){
  scarb1Cpms <- expressionScarb1[grep(x, colnames(expressionScarb1))]
}

scarb1List <- lapply(treatments, splitDF)
names(scarb1List) <- treatments


treatmentMeans <- lapply(scarb1List, function(x){rowMeans(x)})
dfMeans <- bind_rows(treatmentMeans, .id = "treatments")
dfMeans <- as.data.frame(dfMeans)
rownames(dfMeans) <- dfMeans$treatments
dfMeans <- subset(dfMeans, select = -c(1))
deltas <- data.frame()
for (i in 1:8){
  deltas <- rbind(deltas, dfMeans[i,]-dfMeans[9,])
}
# 52283 is the winner
dfNew <- expressionScarb1[3,]
newRow <- gsub(pattern = "\\.\\d", replacement = "", colnames(expressionScarb1))
dfNew <- as.data.frame(t(rbind(dfNew, newRow)))
colnames(dfNew) <- c("cpm", "treatment")


svg(paste0(OUT_DIR, "Scarb1.svg"), width=8, height=4)
ggplot(dfNew, aes(y=cpm, x=treatment)) + 
  geom_boxplot(aes(fill=treatment), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  # facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Limma Normalized Peak Counts (CPM)", title=paste0("Scarb1, Chr5 125339059-125341372")) + 
  # scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()


# Sh3glb1 ----------------------------------------------
peaksSh3glb1 <- annotate[annotate$Gene.Name =="Sh3glb1",]

expressionSh3glb1 <- data.frame()
for (i in 1:length(peaksSh3glb1$PeakID)){
  expressionSh3glb1 <- rbind(expressionSh3glb1, expression[rownames(expression)== peaksSh3glb1$PeakID[i],])
}

# make into a list of averages
treatments <- unique(gsub(pattern = "\\.\\d", replacement = "", colnames(expressionSh3glb1)))

Sh3glb1Cpms <- c()  
splitDF <- function(x){
  Sh3glb1Cpms <- expressionSh3glb1[grep(x, colnames(expressionSh3glb1))]
}

Sh3glb1List <- lapply(treatments, splitDF)
names(Sh3glb1List) <- treatments


treatmentMeans <- lapply(Sh3glb1List, function(x){rowMeans(x)})
dfMeans <- bind_rows(treatmentMeans, .id = "treatments")
dfMeans <- as.data.frame(dfMeans)
rownames(dfMeans) <- dfMeans$treatments
dfMeans <- subset(dfMeans, select = -c(1))
deltas <- data.frame()
for (i in 1:8){
  deltas <- rbind(deltas, dfMeans[i,]-dfMeans[9,])
}
# 61455 is the winner
dfNew <- expressionSh3glb1[3,]
newRow <- gsub(pattern = "\\.\\d", replacement = "", colnames(expressionSh3glb1))
dfNew <- as.data.frame(t(rbind(dfNew, newRow)))
colnames(dfNew) <- c("cpm", "treatment")
dfNew$cpm <- as.numeric(dfNew$cpm)

svg(paste0(OUT_DIR, "Sh3glb1.svg"), width=8, height=4)
ggplot(dfNew, aes(y=cpm, x=treatment)) + 
  geom_boxplot(aes(fill=treatment), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  # facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Limma Normalized Peak Counts (CPM)", title=paste0("Sh3glb1, Chr3 144718853-144719722")) + 
  # scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()
