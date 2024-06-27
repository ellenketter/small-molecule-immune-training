# make a heatmap
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/diffFootprinting/")
molecule <- c("5F", "BG", "Fen", "Fluni", 
              "HC", "HQ", "Myr", "Nerol")

extractStats<- function(x){
  x <- read.delim(paste0(x,"/differential_statistics.txt"))}

list <- lapply(molecule, extractStats)
names(list) <- molecule








dstatNames <- colnames(footStats) 
# protectionScoreNames <- statNames[3:28]
tagCountNames <- statNames[29:54]
moleculeN<- gsub("TC_","",tagCountNames)

# protectionScores <- footStats[3:28]
# tagCounts <- footStats[29:54]

library(dplyr)
library(ggplot2)
library(ggrepel)

footStats <- read.delim("Downloads/diffFootprinting/differential_statistics.txt")

outputs <- colnames(footStats)
moleculeN <- gsub("TC_","",outputs[29:54])
molecule <- c("5F", "BG", "Fen", "Fluni", 
              "HC", "HQ", "Myr", "Nerol")
# calculate activity scores per replicate 
for(i in 3:28){
  footStats <- cbind(footStats, rowSums(footStats[c(i,i+26)]))}
colnames(footStats)[55:80] <- moleculeN
rownames(footStats) <- footStats[,1]
footStats <- subset(footStats, Num > 1000) # removed 5 rows


`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant


subsets <- function(x){
  test <- select(footStats, starts_with(x))
  testMeans<- rowMeans(test)
  pbs <- select(footStats, starts_with("PBS"))
  pbsMeans <- rowMeans(pbs)
  foldChange <- testMeans/pbsMeans
  m <- mean(foldChange)
  s <- sd(foldChange)
  z <- (foldChange-m)/s
  p <- pnorm(z)
  # p_adj <- p.adjust(p)  # not standard to use adj p value according to tutorial
  
  # fcL2 <- log2(foldChange)
  # pL10 <- log10(p)
  x <- cbind(footStats$Motif,foldChange, p)
}
subFootStats<- lapply(molecule, subsets)
names(subFootStats) <- molecule
colnames(subFootStats$BG)<- c("Motif", "foldChange","p")



# old strategy
for(i in 3:28){
  footStats <- cbind(footStats, rowSums(footStats[c(i,i+26)]))}
colnames(footStats)[55:80] <- moleculeN
rownames(footStats) <- footStats[,1]

footStats <- subset(footStats, Num > 1000) # removed 5 rows
df <- subset(footStats, select = moleculeN)

df$Var <- apply(df, 1, sd)

df <- df %>%
  top_n(100, wt = Var) %>%
  select(-Var)

df <- t(scale(t(df)))

options(repr.plot.width = 18, repr.plot.height = 18)

p <- Heatmap(as.matrix(df),
             name = "TF Activity",
             cluster_columns = TRUE,
             cluster_rows = TRUE,
             show_row_names = TRUE,
             rect_gp = gpar(col = "black", lwd = 0.5)
)

pdf("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/diffFootprinting/heatmapReplicates.pdf", 
    width = 8,
    height = 14)
p
dev.off()
