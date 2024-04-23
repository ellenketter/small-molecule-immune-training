library(tidyr)
library(readr)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/")
OUT_DIR <- "background_vs_DA_peaks/"

# Write formatted file for all peaks, regardless of DA status as background peaks.

peakmat <- read.delim(paste0(OUT_DIR,"results_bg.txt"))
peakmat <- data.frame(rownames(peakmat), peakmat$genes)
colnames(peakmat) <- c("Geneid","genes")
peakmat <- separate(data = peakmat, col = genes, into= c("Chr", "Start", "End"))
peakmat <- peakmat[, c(2, 3, 4, 1)]
write.table(peakmat,sep = "\t", row.names = F, quote = F, 
            "background_vs_DA_peaks/backgroundPeaks.txt")


# Write formatted files containing only significantly upregulated genes
### NB results_bg_up etc from 10-1_volcanoPlots.R. Also available in txt files in
### background_vs_DA_peaks.

significant_bg_up <- data.frame(rownames(results_bg_up), results_bg_up$genes)
colnames(significant_bg_up) <- c("Geneid","genes")
significant_bg_up <- separate(data = significant_bg_up, col = genes, into= c("Chr", "Start", "End"))
significant_bg_up <- significant_bg_up[, c(2, 3, 4, 1)]
write.table(significant_bg_up, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/bg_genePeaksUp_p_0.05.txt")

significant_fluni_up <- data.frame(rownames(results_fluni_up), results_fluni_up$genes)
colnames(significant_fluni_up) <- c("Geneid","genes")
significant_fluni_up <- separate(data = significant_fluni_up, col = genes, into= c("Chr", "Start", "End"))
significant_fluni_up <- significant_fluni_up[, c(2, 3, 4, 1)]
write.table(significant_fluni_up, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/fluni_genePeaksUp_p_0.05.txt")

significant_myr_up <- data.frame(rownames(results_myr_up), results_myr_up$genes)
colnames(significant_myr_up) <- c("Geneid","genes")
significant_myr_up <- separate(data = significant_myr_up, col = genes, into= c("Chr", "Start", "End"))
significant_myr_up <- significant_myr_up[, c(2, 3, 4, 1)]
write.table(significant_myr_up, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/myr_genePeaksUp_p_0.05.txt")

significant_hc_up <- data.frame(rownames(results_hc_up), results_hc_up$genes)
colnames(significant_hc_up) <- c("Geneid","genes")
significant_hc_up <- separate(data = significant_hc_up, col = genes, into= c("Chr", "Start", "End"))
significant_hc_up <- significant_hc_up[, c(2, 3, 4, 1)]
write.table(significant_hc_up, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/hc_genePeaksUp_p_0.05.txt")

significant_fen_up <- data.frame(rownames(results_fen_up), results_fen_up$genes)
colnames(significant_fen_up) <- c("Geneid","genes")
significant_fen_up <- separate(data = significant_fen_up, col = genes, into= c("Chr", "Start", "End"))
significant_fen_up <- significant_fen_up[, c(2, 3, 4, 1)]
write.table(significant_fen_up, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/fen_genePeaksUp_p_0.05.txt")


# Write formatted files containing only significantly downregulated genes
### NB results_bg_up etc from 10-1_volcanoPlots.R. Also available in txt files in
### background_vs_DA_peaks.

significant_bg_down <- data.frame(rownames(results_bg_down), results_bg_down$genes)
colnames(significant_bg_down) <- c("Geneid","genes")
significant_bg_down <- separate(data = significant_bg_down, col = genes, into= c("Chr", "Start", "End"))
significant_bg_down <- significant_bg_down[, c(2, 3, 4, 1)]
write.table(significant_bg_down, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/bg_genePeaksDown_p_0.05.txt")

significant_fluni_down <- data.frame(rownames(results_fluni_down), results_fluni_down$genes)
colnames(significant_fluni_down) <- c("Geneid","genes")
significant_fluni_down <- separate(data = significant_fluni_down, col = genes, into= c("Chr", "Start", "End"))
significant_fluni_down <- significant_fluni_down[, c(2, 3, 4, 1)]
write.table(significant_fluni_down, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/fluni_genePeaksDown_p_0.05.txt")

significant_myr_down <- data.frame(rownames(results_myr_down), results_myr_down$genes)
colnames(significant_myr_down) <- c("Geneid","genes")
significant_myr_down <- separate(data = significant_myr_down, col = genes, into= c("Chr", "Start", "End"))
significant_myr_down <- significant_myr_down[, c(2, 3, 4, 1)]
write.table(significant_myr_down, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/myr_genePeaksDown_p_0.05.txt")

significant_hc_down <- data.frame(rownames(results_hc_down), results_hc_down$genes)
colnames(significant_hc_down) <- c("Geneid","genes")
significant_hc_down <- separate(data = significant_hc_down, col = genes, into= c("Chr", "Start", "End"))
significant_hc_down <- significant_hc_down[, c(2, 3, 4, 1)]
write.table(significant_hc_down, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/hc_genePeaksDown_p_0.05.txt")

significant_fen_down <- data.frame(rownames(results_fen_down), results_fen_down$genes)
colnames(significant_fen_down) <- c("Geneid","genes")
significant_fen_down <- separate(data = significant_fen_down, col = genes, into= c("Chr", "Start", "End"))
significant_fen_down <- significant_fen_down[, c(2, 3, 4, 1)]
write.table(significant_fen_down, sep = "\t", row.names = F, quote = F,
            "background_vs_DA_peaks/fen_genePeaksDown_p_0.05.txt")
