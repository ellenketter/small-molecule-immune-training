# Specific library
library(UpSetR)
library(dplyr)
# library(ComplexHeatmap)

# Dataset
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/")
OUT_DIR <- "../10_differentialAccessibility/upsetPlots/"


# Call results from text files ---------------------------------------------

results_5f <- read.table("results_5f.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_bg <- read.table("results_bg.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_fen <- read.table("results_fen.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_fluni <- read.table("results_fluni.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_hc <- read.table("results_hc.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_hq <- read.table("results_hq.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_myr <-read.table("results_myr.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
results_nerol <- read.table("results_nerol.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 

# Reorder results ----------------------------------
results_bg <- results_bg[order(as.numeric(row.names(results_bg))), ]
results_5f <- results_5f[order(as.numeric(row.names(results_5f))), ]
results_fen <- results_fen[order(as.numeric(row.names(results_fen))), ]
results_fluni <- results_fluni[order(as.numeric(row.names(results_fluni))), ] 
results_hc <- results_hc[order(as.numeric(row.names(results_hc))), ]
results_hq <- results_hq[order(as.numeric(row.names(results_hq))), ]
results_myr <- results_myr[order(as.numeric(row.names(results_myr))), ] 
results_nerol <- results_nerol[order(as.numeric(row.names(results_nerol))), ]


listInput <- list(BetaGlucan = results_bg$genes[results_bg$DA != "NO"],
                  # significant_5f = results_5f$genes[results_5f$DA != "NO"],
                  Fenoterol = results_fen$genes[results_fen$DA != "NO"],
                  Flunisolide = results_fluni$genes[results_fluni$DA != "NO"],
                  Hydrocortisone = results_hc$genes[results_hc$DA != "NO"],
                  # significant_hq = results_hq$genes[results_hq$DA != "NO"],
                  Myricetin = results_myr$genes[results_myr$DA != "NO"])
                  # significant_nerol = results_nerol$genes[results_nerol$DA != "NO"],)

tiff("upsetPlot_allDA.tiff", units="in", width=6, height=3, res=500)
upset(fromList(listInput), order.by = "freq",
      mainbar.y.label = "Shared DA Peaks", sets.x.label = "DA Peak Count")
dev.off()

tiff("upsetPlot_allDA_logScale.tiff", units="in", width=6, height=3, res=500)
upset(fromList(listInput), order.by = "freq",
      mainbar.y.label = "Shared DA Peaks: Log Scale", sets.x.label = "DA Peak Count: Log",
      scale.intersections = "log10", scale.sets = "log10")
dev.off()

listInput1 <- list(BetaGlucan = results_bg$genes[results_bg$DA == "UP"],
                  Fenoterol = results_fen$genes[results_fen$DA == "UP"],
                  Flunisolide = results_fluni$genes[results_fluni$DA == "UP"],
                  Hydrocortisone = results_hc$genes[results_hc$DA == "UP"],
                  Myricetin = results_myr$genes[results_myr$DA == "UP"])

tiff("upsetPlot_upDA.tiff", units="in", width=4.5, height=3, res=500)
upset(fromList(listInput1), order.by = "freq",
      mainbar.y.label = "Shared Up DA Peaks", sets.x.label = "Up DA Peak Count")
dev.off()

tiff("upsetPlot_upDA_logScale.tiff", units="in", width=4.5, height=3, res=500)
upset(fromList(listInput1), order.by = "freq",
      mainbar.y.label = "Shared Up DA Peaks: Log Scale", sets.x.label = "Up DA Peak Count: Log",
      scale.intersections = "log10", scale.sets = "log10")
dev.off()
