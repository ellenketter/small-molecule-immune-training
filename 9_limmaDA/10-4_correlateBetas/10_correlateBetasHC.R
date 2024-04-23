##
# Correlate Betas - Based on Sarah Sun's script
# Ellen Ketter
# 8 August 2023
##

library(ggplot2)
library(tidyverse)


setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/")
OUT_DIR <- "../10_differentialAccessibility/correlateBetas/"


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

# Correlate betas ---------------------------------------------------------
treatments <- list(results_5f=results_5f,results_bg=results_bg, results_fen=results_fen,
                   results_fluni= results_fluni, 
                   # results_hc=results_hc, 
                   results_hq=results_hq, results_myr=results_myr,
                   results_nerol=results_nerol)

# confirm concordance between gene names, should print "integer(0)" for all
for (i in 1:length(treatments)){
   not <- which(results_hc$genes != treatments[i]$genes)
   print(not)
}

da_hc <- results_hc[which(results_hc$DA != "NO"),]
da_hc <- da_hc %>%
   select(genes, logFC)
da_hc$peakID <- rownames(da_hc)

columnselect <- function(X){
   X[,c("genes","logFC")]
}
mylist <- lapply(X=treatments,FUN=columnselect)
names(mylist) <- gsub("results_","",names(mylist))

allLogFC <- Reduce(function(x, y) merge(x, y, by = "genes", all = TRUE), mylist)
# allLogFC <- mylist %>% reduce(inner_join, by='genes')
allLogFC <- merge(da_hc, allLogFC, by= 'genes')
allLogFC_names <- c("genes","hc","peakID", names(mylist))
# allLogFC_names <- gsub("results_","",allLogFC_names)
colnames(allLogFC) <- allLogFC_names
allLogFC <- allLogFC %>% select(peakID, genes, everything())
allLogFC <- allLogFC[order(as.numeric(allLogFC$peakID)), ]

# Plot dataframes 

for (i in 1:length(names(mylist))){
   
   i_file_name <- paste0("beta_correlation_hc_DA_",names(mylist)[i],"_logFC.tiff")

   p <- allLogFC %>%  
      ggplot(aes(x=.data[["hc"]], y= .data[[names(mylist)[i]]])) + 
      geom_point(shape=21, color="black", size=0.1) + 
      theme_classic() + 
      labs(x="Hydrocortisone LogFC", y=paste0(names(mylist)[i], " LogFC"), fill="")+
      geom_hline(yintercept=0) + 
      geom_vline(xintercept=0) + 
      geom_smooth(method = "lm", se=FALSE) +
      theme(axis.text=element_text(size=15.5),
            axis.title=element_text(size=15.5)) + 
      theme(legend.text = element_text(size=15.5))
 
 ggsave(filename = i_file_name, plot = p, width=4.5, height=3)

}
