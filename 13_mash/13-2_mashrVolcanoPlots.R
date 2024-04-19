library(UpSetR)
library(dplyr)
library(ggplot2)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "volcanoPlots/"
# dir.create(OUT_DIR)

# Call results from text files ---------------------------------------------

mashMeans <- readRDS("mashMeans_dd.rds") 

treatments <- unique(mashMeans$treatment)
getPeaks <- function(x){mashMeans[mashMeans$treatment == x,]}
mashPeaks <- lapply(treatments, getPeaks)
names(mashPeaks) <- treatments

png(file= "volcanoPlots/results_5f_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$`5f`, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: 5F') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_bg_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$bg, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Beta Glucan') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_fen_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$fen, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Fenoterol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_fluni_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$fluni, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Flunisolide') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_hc_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$hc, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydrocortisone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_hq_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$hq, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydroquinone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

svg(file= "volcanoPlots/results_myr_0.05.svg")
ggplot(data = mashPeaks$myr, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Myricetin') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "volcanoPlots/results_nerol_0.05.png", width = 600, height = 400)
ggplot(data = mashPeaks$nerol, aes(x = posteriorMean, y = -log10(lfsr), col = significant)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Nerol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

# Barplot of DA Peaks: Up vs Down -----------------------------------------


bg <- c(sum(mashPeaks$bg$significant == "UP"), sum(mashPeaks$bg$significant == "DOWN"))
fiveF <- c(sum(mashPeaks$'5f'$significant == "UP"), sum(mashPeaks$'5f'$significant == "DOWN"))
fen <- c(sum(mashPeaks$fen$significant == "UP"), sum(mashPeaks$fen$significant == "DOWN"))
fluni <- c(sum(mashPeaks$fluni$significant == "UP"), sum(mashPeaks$fluni$significant == "DOWN"))
hc <- c(sum(mashPeaks$hc$significant == "UP"), sum(mashPeaks$hc$significant == "DOWN"))
hq <- c(sum(mashPeaks$hq$significant == "UP"), sum(mashPeaks$hq$significant == "DOWN"))
myr <- c(sum(mashPeaks$myr$significant == "UP"), sum(mashPeaks$myr$significant == "DOWN"))
nerol <- c(sum(mashPeaks$nerol$significant == "UP"), sum(mashPeaks$nerol$significant == "DOWN"))

countsMash <- data.frame(bg, fluni, myr, hc, fen, hq, nerol, fiveF)
rownames(countsMash) <- c("Upregulated", "Downregulated")
countsMash <- reshape2::melt(as.matrix(countsMash)) 
colnames(countsMash) <- c("Direction","Treatment","PeakCount")

png(file= "volcanoPlots/totalPeakCount_lfsr_0.05.png", width = 500, height = 400)
ggplot(countsMash, aes(fill=Direction, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle('Number of Differentially Accessible Peaks by Mash (LFSR < 0.05) ') +
  ylim(0,50000)
dev.off()

# edited for figure 4

BetaGlucan <- c(sum(mashPeaks$bg$significant == "UP"), -sum(mashPeaks$bg$significant == "DOWN"))
FluoroindoleCarboxylicAcid <- c(sum(mashPeaks$'5f'$significant == "UP"), -sum(mashPeaks$'5f'$significant == "DOWN"))
Fenoterol <- c(sum(mashPeaks$fen$significant == "UP"), -sum(mashPeaks$fen$significant == "DOWN"))
Flunisolide <- c(sum(mashPeaks$fluni$significant == "UP"), -sum(mashPeaks$fluni$significant == "DOWN"))
Hydrocortisone <- c(sum(mashPeaks$hc$significant == "UP"), -sum(mashPeaks$hc$significant == "DOWN"))
Hydroquinone <- c(sum(mashPeaks$hq$significant == "UP"), -sum(mashPeaks$hq$significant == "DOWN"))
Myricetin <- c(sum(mashPeaks$myr$significant == "UP"), -sum(mashPeaks$myr$significant == "DOWN"))
Nerol <- c(sum(mashPeaks$nerol$significant == "UP"), -sum(mashPeaks$nerol$significant == "DOWN"))

countsMash <- data.frame(Flunisolide,  Hydrocortisone, Myricetin, Fenoterol, 
                         Nerol, FluoroindoleCarboxylicAcid, Hydroquinone,  BetaGlucan)
rownames(countsMash) <- c("Upregulated", "Downregulated")
countsMash <- reshape2::melt(as.matrix(countsMash)) 
colnames(countsMash) <- c("Direction","Treatment","PeakCount")

svg(file= "volcanoPlots/summaryPeaksDA.svg", 
    width = 500, height = 400)
ggplot(countsMash, aes(fill=Treatment, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  geom_text(aes(label = PeakCount), position = position_stack(vjust = 0.5), size = 3, colour = "white") + 
  scale_fill_manual(values = c("#2BC4E4", "#F9987A", "#485695","#7C543E", 
                               "thistle4", "#11A290","#7BD1C1", "#B3B7B9") ) +
  geom_hline(aes(yintercept = 0)) + 
  theme_classic() +
  theme(legend.position="none",
        axis.text.x = element_text(angle = 45, hjust =1, vjust=1)) +
  ylab("Number of Differentially 
Accessible Peaks by Mashr

Downregulated      Upregulated") +
  labs(subtitle = '    Steroids                      Unrelated              Pathogen-
                                   Small Molecules       Associated') +
  xlab("") +
  ylim(-20000,20000) 
  
dev.off()

