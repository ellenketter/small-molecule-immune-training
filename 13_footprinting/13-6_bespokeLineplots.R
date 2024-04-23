library(dplyr)
library(ggplot2)


setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/diffFootprinting/multi/Lineplots/")
OUT_DIR <- "bespokeLineplots/"
# dir.create(OUT_DIR)
mashInput <- readRDS("/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/16_footprinting/figure4_TFs.rds")

allFiles <- list.files(pattern = "txt")

master <- list()
for (i in 1:length(mashInput)){
  df <- read.delim(allFiles[grep(paste0("\\.", mashInput[i],".txt"), allFiles)])
  master[[i]] <- df
         }
names(master) <- mashInput

# saveRDS(master, "lineplotDF_SignifcantTFs.rds")
master <- readRDS("lineplotDF_SignifcantTFs.rds")

library(ggplot2)
library(dplyr)
library(tidyr)
# library(hrbrthemes)
library(viridis)
library(reshape2)
library(gridExtra)
library(sjplot)
# Keep only 3 names
test <- master
test <- melt(test)
nreps <- nrow(test)/200
nt <- rep(-100:99,nreps)
test <- cbind(test, nt)


head(test)
test <- split(test, test$L1)

p <- lapply(test, function(m) {ggplot(m, aes(x = nt, y = value, group=variable, color=variable)) + 
  geom_line() +
  scale_color_manual(values = c("#11A290","#B3B7B9", "#7C543E", "#2BC4E4", "#F9987A","#7BD1C1", "#485695", 
                                "thistle4", "black")) +
  ggtitle(paste0('count of ',unique(m$L1))) +
  # theme_ipsum() +
  ylab("Normalized Number of Reads")+ 
  xlab("Coordinates From Motif Center")})

pdf("test.pdf", width = 30, height = 40)
plot_grid(plotlist = p, nrow = 10, ncol = 4)
dev.off()
# 1. change molecule names to full 
# 2. make single legend
# 3. install hrbr, fix background to white
# 4. make second page with table of jaspar motifs and n of motif matches

