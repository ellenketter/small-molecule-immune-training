library("ggplot2")
library ("DESeq2")
library("pheatmap")
library("wesanderson")
library(statmod)
library(goseq)
library("RColorBrewer")
library(edgeR)
library(ggfortify)
library(devtools)
library("FactoMineR")
library("factoextra")
library(tidyverse)
library(cluster)
library(factoextra)
library(ggrepel)
library(data.table)
library(VennDiagram)


setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks/")
OUT_DIR <- "../10_differentialAccessibility/"

# import raw counts -------------------------------------------------------
expression <- read.table(paste0(OUT_DIR,"expressionCPMs_allPBS_processOrder.txt"))
fit <- readRDS(paste0(OUT_DIR,"fit_allPBS_processOrder.rds"))

fitPeaks <- fit$genes
# write.table(fitPeaks, file="backgroundPeaks_allPBS_processOrder.txt")


# Make treatment-specific results table-----------------------------------
# stims <- colnames(fit)
# stims <- stims[stims != "intercept"]

# res_list <- lapply(stims, function(x){
#   topTable(fit, coef = x, n = Inf, sort = "p")
#  })




results_5f <- topTable(fit, coef = "5f", n = Inf, sort = "p")
results_bg <- topTable(fit, coef = "bg", n = Inf, sort = "p")
results_fen <- topTable(fit, coef = "fen", n = Inf, sort = "p")
results_fluni <- topTable(fit, coef = "fluni", n = Inf, sort = "p")
results_hc <- topTable(fit, coef = "hc", n = Inf, sort = "p")
results_hq <- topTable(fit, coef = "hq", n = Inf, sort = "p")
results_myr <- topTable(fit, coef = "myr", n = Inf, sort = "p")
results_nerol <- topTable(fit, coef = "nerol", n = Inf, sort = "p")

results_bg$DA <- "NO"
results_bg$DA[results_bg$logFC > 0 & results_bg$adj.P.Val < 0.05] <- "UP"
results_bg$DA[results_bg$logFC < 0 & results_bg$adj.P.Val < 0.05] <- "DOWN"

results_5f$DA <- "NO"
results_5f$DA[results_5f$logFC > 0 &results_5f$adj.P.Val < 0.05] <- "UP"
results_5f$DA[results_5f$logFC < 0 &results_5f$adj.P.Val < 0.05] <- "DOWN"

results_fen$DA <- "NO"
results_fen$DA[results_fen$logFC > 0 & results_fen$adj.P.Val < 0.05] <- "UP"
results_fen$DA[results_fen$logFC < 0 & results_fen$adj.P.Val < 0.05] <- "DOWN"

results_fluni$DA <- "NO"
results_fluni$DA[results_fluni$logFC > 0 & results_fluni$adj.P.Val < 0.05] <- "UP"
results_fluni$DA[results_fluni$logFC < 0 & results_fluni$adj.P.Val < 0.05] <- "DOWN"

results_hc$DA <- "NO"
results_hc$DA[results_hc$logFC > 0 & results_hc$adj.P.Val < 0.05] <- "UP"
results_hc$DA[results_hc$logFC < 0 & results_hc$adj.P.Val < 0.05] <- "DOWN"

results_hq$DA <- "NO"
results_hq$DA[results_hq$logFC > 0 & results_hq$adj.P.Val < 0.05] <- "UP"
results_hq$DA[results_hq$logFC < 0 & results_hq$adj.P.Val < 0.05] <- "DOWN"

results_myr$DA <- "NO"
results_myr$DA[results_myr$logFC > 0 & results_myr$adj.P.Val < 0.05] <- "UP"
results_myr$DA[results_myr$logFC < 0 & results_myr$adj.P.Val < 0.05] <- "DOWN"

results_nerol$DA <- "NO"
results_nerol$DA[results_nerol$logFC > 0 & results_nerol$adj.P.Val < 0.05] <- "UP"
results_nerol$DA[results_nerol$logFC < 0 & results_nerol$adj.P.Val < 0.05] <- "DOWN"

head(results_bg[order(results_bg$adj.P.Val) & results_bg$DA == 'NO', ])

write.table(results_5f,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_5f.txt")
write.table(results_bg,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_bg.txt")
write.table(results_fen,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fen.txt")
write.table(results_fluni,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fluni.txt")
write.table(results_hc,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hc.txt")
write.table(results_hq,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hq.txt")
write.table(results_myr,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_myr.txt")
write.table(results_nerol,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_nerol.txt")

results_bg_up <- results_bg[results_bg$DA == "UP",]
results_5f_up <- results_5f[results_5f$DA == "UP",]
results_fen_up <- results_fen[results_fen$DA == "UP",]
results_fluni_up <- results_fluni[results_fluni$DA == "UP",]
results_hc_up <- results_hc[results_hc$DA == "UP",]
results_hq_up <- results_hq[results_hq$DA == "UP",]
results_myr_up <- results_myr[results_myr$DA == "UP",]
results_nerol_up <- results_nerol[results_nerol$DA == "UP",]

write.table(results_5f_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_5f_up.txt")
write.table(results_bg_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_bg_up.txt")
write.table(results_fen_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fen_up.txt")
write.table(results_fluni_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fluni_up.txt")
write.table(results_hc_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hc_up.txt")
write.table(results_hq_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hq_up.txt")
write.table(results_myr_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_myr_up.txt")
write.table(results_nerol_up,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_nerol_up.txt")

results_bg_down <- results_bg[results_bg$DA == "DOWN",]
results_5f_down <- results_5f[results_5f$DA == "DOWN",]
results_fen_down <- results_fen[results_fen$DA == "DOWN",]
results_fluni_down <- results_fluni[results_fluni$DA == "DOWN",]
results_hc_down <- results_hc[results_hc$DA == "DOWN",]
results_hq_down <- results_hq[results_hq$DA == "DOWN",]
results_myr_down <- results_myr[results_myr$DA == "DOWN",]
results_nerol_down <- results_nerol[results_nerol$DA == "DOWN",]

write.table(results_5f_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_5f_down.txt")
write.table(results_bg_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_bg_down.txt")
write.table(results_fen_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fen_down.txt")
write.table(results_fluni_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_fluni_down.txt")
write.table(results_hc_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hc_down.txt")
write.table(results_hq_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_hq_down.txt")
write.table(results_myr_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_myr_down.txt")
write.table(results_nerol_down,sep = "\t",  "../10_differentialAccessibility/background_vs_DA_peaks/results_nerol_down.txt")

# Volcano Plots -------------------------------------
png(file= "../10_differentialAccessibility/results_bg_0.05.png", width = 600, height = 400)
ggplot(data = results_bg, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = 0, col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Beta Glucan') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "../10_differentialAccessibility/results_5f_0.05.png", width = 600, height = 400)
ggplot(data = results_5f, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: 5F') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_fen_0.05.png", width = 600, height = 400)
ggplot(data = results_fen, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Fenoterol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "../10_differentialAccessibility/results_fluni_0.05.png", width = 600, height = 400)
ggplot(data = results_fluni, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Flunisolide') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_hc_0.05.png", width = 600, height = 400)
ggplot(data = results_hc, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydrocortisone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_hq_0.05.png", width = 600, height = 400)
ggplot(data = results_hq, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydroquinone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_myr_0.05.png", width = 600, height = 400)
ggplot(data = results_myr, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Myricetin') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_nerol_0.05.png", width = 600, height = 400)
ggplot(data = results_nerol, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Nerol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

# with help from this plot tutorial https://biostatsquid.com/volcano-plots-r-tutorial/

# Barplot of DA Peaks: Up vs Down -----------------------------------------

bg <- c(sum(results_bg$DA == "UP"),sum(results_bg$DA == "DOWN"))
fiveF <- c(sum(results_5f$DA == "UP"),sum(results_5f$DA == "DOWN"))
fen <- c(sum(results_fen$DA == "UP"),sum(results_fen$DA == "DOWN"))
fluni <- c(sum(results_fluni$DA == "UP"),sum(results_fluni$DA == "DOWN"))
hc <- c(sum(results_hc$DA == "UP"),sum(results_hc$DA == "DOWN"))
hq <- c(sum(results_hq$DA == "UP"),sum(results_hq$DA == "DOWN"))
myr <- c(sum(results_myr$DA == "UP"),sum(results_myr$DA == "DOWN"))
nerol <- c(sum(results_nerol$DA == "UP"),sum(results_nerol$DA == "DOWN"))

countsDA <- data.frame(bg, fluni, myr, hc, fen, hq, nerol, fiveF)
rownames(countsDA) <- c("Upregulated", "Downregulated")
countsDA <- reshape2::melt(as.matrix(countsDA)) 
colnames(countsDA) <- c("Direction","Treatment","PeakCount")

png(file= "../10_differentialAccessibility/totalPeakCount_p_0.05.png", width = 500, height = 400)
ggplot(countsDA, aes(fill=Direction, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  # scale_y_log10() +
  ggtitle('Number of Differentially Accessible Peaks Per Condition (p_adjusted < 0.05) ')
dev.off()



png(file= "../10_differentialAccessibility/totalPeakCount_p_0.05_log.png", width = 500, height = 400)
ggplot(countsDA, aes(fill=Direction, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  scale_y_log10() +
  ggtitle('Number of Differentially Accessible Peaks Per Condition (p_adjusted < 0.05) ')
dev.off()


countsDA2 <- data.frame(fluni, myr, hc, fen, hq, nerol, fiveF)
rownames(countsDA2) <- c("Upregulated", "Downregulated")
countsDA2 <- reshape2::melt(as.matrix(countsDA2)) 
colnames(countsDA2) <- c("Direction","Treatment","PeakCount")

png(file= "../10_differentialAccessibility/totalPeakCount_p_0.05_noBG.png", width = 500, height = 400)
ggplot(countsDA2, aes(fill=Direction, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  # scale_y_log10() +
  ggtitle('Number of Differentially Accessible Peaks Per Condition (p_adjusted < 0.05) ')
dev.off()