# make a heatmap
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(RColorBrewer))
suppressMessages(library(tidyverse))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))

footStats <- read.delim("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/diffFootprinting/2D/differential_statistics.txt")
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/diffFootprinting/2D/")


library(ggrepel)

ggplot(data=footStats, aes(x=log2(Protection_Score_Nerol/Protection_Score_PBS), y=-log10(P_values), label = Motif)) +
  geom_point(size=.5) +
  theme_minimal() +
  geom_text_repel()













statNames <- colnames(footStats) 
# protectionScoreNames <- statNames[3:28]
tagCountNames <- statNames[29:54]
moleculeN<- gsub("TC_","",tagCountNames)

# protectionScores <- footStats[3:28]
# tagCounts <- footStats[29:54]



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
