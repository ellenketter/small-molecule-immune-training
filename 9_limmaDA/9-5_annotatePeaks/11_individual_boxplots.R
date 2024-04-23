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
annotate <- read.delim("background_annotate_out.txt")
annotate_names <- colnames(annotate)
annotate_names[1] <- "PeakID"
colnames(annotate) <- annotate_names

# Beta Glucan -------------------------------------------------------------
#load all limma results files

limma_bg_up <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_bg_up.txt")
limma_bg_down <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_bg_down.txt")
top10_up_bg <- head(limma_bg_up[order(limma_bg_up$adj.P.Val),],10)
top10_down_bg <- head(limma_bg_down[order(limma_bg_down$adj.P.Val),],10)

# Beta Glucan - UP --------------------------
# sampleN, peakID, expression (CPM), gene.name, gene.description
bg_peaks_up <- rownames(top10_up_bg)
bg_expression <- data.frame()
bg_peak2gene <- data.frame()
for (i in 1:length(bg_peaks_up)){
  bg_expression <- rbind(bg_expression, expression[rownames(expression)== bg_peaks_up[i],])
  bg_peak2gene <- rbind(bg_peak2gene, annotate[annotate$PeakID == bg_peaks_up[i],])
}
bg_expression <- bg_expression[c("BG.1","BG.3","PBS.1","PBS.2","PBS.3")]
# bg_expression <- t(bg_expression)
bg_expression <- cbind(bg_peaks_up, bg_expression)
colnames(bg_expression)[1] <- "PeakID"
bg_peak2gene <- bg_peak2gene[c(1,16)]
betaGlucan <- merge(bg_expression, bg_peak2gene, 
                    by = "PeakID")

betaGlucan2 <- reshape2::melt(betaGlucan, id="Gene.Name")
betaGlucan3<- mutate(betaGlucan2, variable = gsub(pattern = ".[0-9]", "",variable)) 
betaGlucan3 <- betaGlucan3[-c(1:10),]
betaGlucan3$value <- as.numeric(betaGlucan3$value)

tiff(paste0(OUT_DIR, "bgUpBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(betaGlucan3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Beta Glucan: Upregulated")) + 
  scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Beta Glucan - Down --------------------------
bg_peaks_down <- rownames(top10_down_bg)
bg_expression <- data.frame()
bg_peak2gene <- data.frame()
for (i in 1:length(bg_peaks_down)){
  bg_expression <- rbind(bg_expression, expression[rownames(expression)== bg_peaks_down[i],])
  bg_peak2gene <- rbind(bg_peak2gene, annotate[annotate$PeakID == bg_peaks_down[i],])
}
bg_expression <- bg_expression[c("BG.1","BG.3","PBS.1","PBS.2","PBS.3")]
# bg_expression <- t(bg_expression)
bg_expression <- cbind(bg_peaks_down, bg_expression)
colnames(bg_expression)[1] <- "PeakID"
bg_peak2gene <- bg_peak2gene[c(1,16)]
betaGlucan <- merge(bg_expression, bg_peak2gene, 
                    by = "PeakID")

betaGlucan2 <- reshape2::melt(betaGlucan, id="Gene.Name")
betaGlucan3<- mutate(betaGlucan2, variable = gsub(pattern = ".[0-9]", "",variable)) 
betaGlucan3 <- betaGlucan3[-c(1:10),]
betaGlucan3$value <- as.numeric(betaGlucan3$value)

tiff(paste0(OUT_DIR, "bgDownBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(betaGlucan3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Beta Glucan: Downregulated")) + 
  scale_fill_manual(values=c("#9983BD", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()



# Fenoterol -------------------------------------------------------------
#load all limma results files

limma_fen_up <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_fen_up.txt")
limma_fen_down <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_fen_down.txt")
top10_up_fen <- head(limma_fen_up[order(limma_fen_up$adj.P.Val),],10)
top10_down_fen <- head(limma_fen_down[order(limma_fen_down$adj.P.Val),],10)

# Fenoterol - UP --------------------------
# sampleN, peakID, expression (CPM), gene.name, gene.description
fen_peaks_up <- rownames(top10_up_fen)
fen_expression <- data.frame()
fen_peak2gene <- data.frame()
for (i in 1:length(fen_peaks_up)){
  fen_expression <- rbind(fen_expression, expression[rownames(expression)== fen_peaks_up[i],])
  fen_peak2gene <- rbind(fen_peak2gene, annotate[annotate$PeakID == fen_peaks_up[i],])
}
fen_expression <- fen_expression[c("Fen.1","Fen.2","Fen.3","PBS.1","PBS.2","PBS.3")]

fen_expression <- cbind(fen_peaks_up, fen_expression)
colnames(fen_expression)[1] <- "PeakID"
fen_peak2gene <- fen_peak2gene[c(1,16)]
fenoterol <- merge(fen_expression, fen_peak2gene, 
                    by = "PeakID")

fenoterol2 <- reshape2::melt(fenoterol, id="Gene.Name")
fenoterol3<- mutate(fenoterol2, variable = gsub(pattern = ".[0-9]", "",variable)) 
fenoterol3 <- fenoterol3[-c(1:10),]
fenoterol3$value <- as.numeric(fenoterol3$value)

tiff(paste0(OUT_DIR, "fenUpBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(fenoterol3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Fenoterol: Upregulated")) + 
  scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Fenoterol - Down --------------------------
fen_peaks_down <- rownames(top10_down_fen)
fen_expression <- data.frame()
fen_peak2gene <- data.frame()
for (i in 1:length(fen_peaks_down)){
  fen_expression <- rbind(fen_expression, expression[rownames(expression)== fen_peaks_down[i],])
  fen_peak2gene <- rbind(fen_peak2gene, annotate[annotate$PeakID == fen_peaks_down[i],])
}
fen_expression <- fen_expression[c("Fen.1","Fen.2","Fen.3","PBS.1","PBS.2","PBS.3")]
# fen_expression <- t(fen_expression)
fen_expression <- cbind(fen_peaks_down, fen_expression)
colnames(fen_expression)[1] <- "PeakID"
fen_peak2gene <- fen_peak2gene[c(1,16)]
fenoterol <- merge(fen_expression, fen_peak2gene, 
                    by = "PeakID")

fenoterol2 <- reshape2::melt(fenoterol, id="Gene.Name")
fenoterol3<- mutate(fenoterol2, variable = gsub(pattern = ".[0-9]", "",variable)) 
fenoterol3 <- fenoterol3[-c(1:10),]
fenoterol3$value <- as.numeric(fenoterol3$value)

tiff(paste0(OUT_DIR, "fenDownBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(fenoterol3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Fenoterol: Downregulated")) + 
  scale_fill_manual(values=c("#9983BD", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Flunisolide -------------------------------------------------------------
#load all limma results files

limma_fluni_up <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_fluni_up.txt")
limma_fluni_down <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_fluni_down.txt")
top10_up_fluni <- head(limma_fluni_up[order(limma_fluni_up$adj.P.Val),],10)
top10_down_fluni <- head(limma_fluni_down[order(limma_fluni_down$adj.P.Val),],10)

# Flunisolide - UP --------------------------
# sampleN, peakID, expression (CPM), gene.name, gene.description
fluni_peaks_up <- rownames(top10_up_fluni)
fluni_expression <- data.frame()
fluni_peak2gene <- data.frame()
for (i in 1:length(fluni_peaks_up)){
  fluni_expression <- rbind(fluni_expression, expression[rownames(expression)== fluni_peaks_up[i],])
  fluni_peak2gene <- rbind(fluni_peak2gene, annotate[annotate$PeakID == fluni_peaks_up[i],])
}
fluni_expression <- fluni_expression[c("Fluni.1","Fluni.2","PBS.1","PBS.2","PBS.3")]
# fluni_expression <- t(fluni_expression)
fluni_expression <- cbind(fluni_peaks_up, fluni_expression)
colnames(fluni_expression)[1] <- "PeakID"
fluni_peak2gene <- fluni_peak2gene[c(1,16)]
flunisolide <- merge(fluni_expression, fluni_peak2gene, 
                    by = "PeakID")

flunisolide2 <- reshape2::melt(flunisolide, id="Gene.Name")
flunisolide3<- mutate(flunisolide2, variable = gsub(pattern = ".[0-9]", "",variable)) 
flunisolide3 <- flunisolide3[-c(1:10),]
flunisolide3$value <- as.numeric(flunisolide3$value)

tiff(paste0(OUT_DIR, "fluniUpBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(flunisolide3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Flunisolide: Upregulated")) + 
  scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Flunisolide - Down --------------------------
fluni_peaks_down <- rownames(top10_down_fluni)
fluni_expression <- data.frame()
fluni_peak2gene <- data.frame()
for (i in 1:length(fluni_peaks_down)){
  fluni_expression <- rbind(fluni_expression, expression[rownames(expression)== fluni_peaks_down[i],])
  fluni_peak2gene <- rbind(fluni_peak2gene, annotate[annotate$PeakID == fluni_peaks_down[i],])
}
fluni_expression <- fluni_expression[c("Fluni.1","Fluni.2","PBS.1","PBS.2","PBS.3")]

fluni_expression <- cbind(fluni_peaks_down, fluni_expression)
colnames(fluni_expression)[1] <- "PeakID"
fluni_peak2gene <- fluni_peak2gene[c(1,16)]
flunisolide <- merge(fluni_expression, fluni_peak2gene, 
                    by = "PeakID")

flunisolide2 <- reshape2::melt(flunisolide, id="Gene.Name")
flunisolide3<- mutate(flunisolide2, variable = gsub(pattern = ".[0-9]", "",variable)) 
flunisolide3 <- flunisolide3[-c(1:10),]
flunisolide3$value <- as.numeric(flunisolide3$value)

tiff(paste0(OUT_DIR, "fluniDownBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(flunisolide3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Flunisolide: Downregulated")) + 
  scale_fill_manual(values=c("#9983BD", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Hydrocortisone -------------------------------------------------------------
#load all limma results files

limma_hc_up <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_hc_up.txt")
limma_hc_down <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_hc_down.txt")

top10_up_hc <- head(limma_hc_up[order(limma_hc_up$adj.P.Val),],10)
top10_down_hc <- head(limma_hc_down[order(limma_hc_down$adj.P.Val),],10)

# Hydrocortisone - UP --------------------------
# sampleN, peakID, expression (CPM), gene.name, gene.description
hc_peaks_up <- rownames(top10_up_hc)
hc_expression <- data.frame()
hc_peak2gene <- data.frame()
for (i in 1:length(hc_peaks_up)){
  hc_expression <- rbind(hc_expression, expression[rownames(expression)== hc_peaks_up[i],])
  hc_peak2gene <- rbind(hc_peak2gene, annotate[annotate$PeakID == hc_peaks_up[i],])
}
hc_expression <- hc_expression[c("HC.1","HC.2","PBS.1","PBS.2","PBS.3")]
# hc_expression <- t(hc_expression)
hc_expression <- cbind(hc_peaks_up, hc_expression)
colnames(hc_expression)[1] <- "PeakID"
hc_peak2gene <- hc_peak2gene[c(1,16)]
hydrocortisone <- merge(hc_expression, hc_peak2gene, 
                    by = "PeakID")

hydrocortisone2 <- reshape2::melt(hydrocortisone, id="Gene.Name")
hydrocortisone3<- mutate(hydrocortisone2, variable = gsub(pattern = ".[0-9]", "",variable)) 
hydrocortisone3 <- hydrocortisone3[-c(1:10),]
hydrocortisone3$value <- as.numeric(hydrocortisone3$value)

tiff(paste0(OUT_DIR, "hcUpBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(hydrocortisone3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Hydrocortisone: Upregulated")) + 
  scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Hydrocortisone - Down --------------------------
hc_peaks_down <- rownames(top10_down_hc)
hc_expression <- data.frame()
hc_peak2gene <- data.frame()
for (i in 1:length(hc_peaks_down)){
  hc_expression <- rbind(hc_expression, expression[rownames(expression)== hc_peaks_down[i],])
  hc_peak2gene <- rbind(hc_peak2gene, annotate[annotate$PeakID == hc_peaks_down[i],])
}
hc_expression <- hc_expression[c("HC.1","HC.2","PBS.1","PBS.2","PBS.3")]
# hc_expression <- t(hc_expression)
hc_expression <- cbind(hc_peaks_down, hc_expression)
colnames(hc_expression)[1] <- "PeakID"
hc_peak2gene <- hc_peak2gene[c(1,16)]
hydrocortisone <- merge(hc_expression, hc_peak2gene, 
                    by = "PeakID")

hydrocortisone2 <- reshape2::melt(hydrocortisone, id="Gene.Name")
hydrocortisone3<- mutate(hydrocortisone2, variable = gsub(pattern = ".[0-9]", "",variable)) 
hydrocortisone3 <- hydrocortisone3[-c(1:10),]
hydrocortisone3$value <- as.numeric(hydrocortisone3$value)

tiff(paste0(OUT_DIR, "hcDownBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(hydrocortisone3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Hydrocortisone: Downregulated")) + 
  scale_fill_manual(values=c("#9983BD", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Myricetin -------------------------------------------------------------
#load all limma results files

limma_myr_up <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_myr_up.txt")
limma_myr_down <- read.delim("../10_differentialAccessibility/background_vs_DA_peaks/results_myr_down.txt")

top10_up_myr <- head(limma_myr_up[order(limma_myr_up$adj.P.Val),],10)
top10_down_myr <- head(limma_myr_down[order(limma_myr_down$adj.P.Val),],10)

# Myricetin - UP --------------------------
# sampleN, peakID, expression (CPM), gene.name, gene.description
myr_peaks_up <- rownames(top10_up_myr)
myr_expression <- data.frame()
myr_peak2gene <- data.frame()
for (i in 1:length(myr_peaks_up)){
  myr_expression <- rbind(myr_expression, expression[rownames(expression)== myr_peaks_up[i],])
  myr_peak2gene <- rbind(myr_peak2gene, annotate[annotate$PeakID == myr_peaks_up[i],])
}
myr_expression <- myr_expression[c("Myr.1","Myr.2","Myr.3","PBS.1","PBS.2","PBS.3")]
# myr_expression <- t(myr_expression)
myr_expression <- cbind(myr_peaks_up, myr_expression)
colnames(myr_expression)[1] <- "PeakID"
myr_peak2gene <- myr_peak2gene[c(1,16)]
myricetin <- merge(myr_expression, myr_peak2gene, 
                    by = "PeakID")

myricetin2 <- reshape2::melt(myricetin, id="Gene.Name")
myricetin3<- mutate(myricetin2, variable = gsub(pattern = ".[0-9]", "",variable)) 
myricetin3 <- myricetin3[-c(1:10),]
myricetin3$value <- as.numeric(myricetin3$value)

tiff(paste0(OUT_DIR, "myrUpBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(myricetin3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Myricetin: Upregulated")) + 
  scale_fill_manual(values=c("red", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

# Myricetin - Down --------------------------
myr_peaks_down <- rownames(top10_down_myr)
myr_expression <- data.frame()
myr_peak2gene <- data.frame()
for (i in 1:length(myr_peaks_down)){
  myr_expression <- rbind(myr_expression, expression[rownames(expression)== myr_peaks_down[i],])
  myr_peak2gene <- rbind(myr_peak2gene, annotate[annotate$PeakID == myr_peaks_down[i],])
}
myr_expression <- myr_expression[c("Myr.1","Myr.2","Myr.3","PBS.1","PBS.2","PBS.3")]
# myr_expression <- t(myr_expression)
myr_expression <- cbind(myr_peaks_down, myr_expression)
colnames(myr_expression)[1] <- "PeakID"
myr_peak2gene <- myr_peak2gene[c(1,16)]
myricetin <- merge(myr_expression, myr_peak2gene, 
                    by = "PeakID")

myricetin2 <- reshape2::melt(myricetin, id="Gene.Name")
myricetin3<- mutate(myricetin2, variable = gsub(pattern = ".[0-9]", "",variable)) 
myricetin3 <- myricetin3[-c(1:10),]
myricetin3$value <- as.numeric(myricetin3$value)

tiff(paste0(OUT_DIR, "myrDownBoxplots.tiff"), units="in", width=8, height=4, res=250)
ggplot(myricetin3, aes(y=value, x=variable)) + 
  geom_boxplot(aes(fill=variable), alpha=0.5) +
  geom_jitter(shape=21, size=1, color="black",alpha=1, position=position_jitter(0.1)) + 
  facet_wrap(~Gene.Name, scales = 'free', nrow = 2) + 
  labs(x="", y="Normalized Peak Counts", title=paste0("Top 10 DA Peaks in Myricetin: Downregulated")) + 
  scale_fill_manual(values=c("#9983BD", "grey"))+
  theme(legend.position = "none") +
  theme_classic()
dev.off()

