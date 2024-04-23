library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/otherModels/dataDriven_commonBaseline/")


res_list <- readRDS("mashMeans_dd.rds")

stims <- levels(res_list$treatment)

res_list <- res_list %>% separate(peak, c('chromosome', 'start', 'end'))
res_list <- cbind(res_list, paste0(res_list$chromosome,".", res_list$start,".", res_list$end))
colnames(res_list)[8] <- "peak"

subsetHomer<- function(x){subset(res_list, treatment == x & significant == 'UP') %>%
    select(c("chromosome","start", "end", "peak"))}
homer_res_list <- lapply(stims, subsetHomer)
names(homer_res_list) <- stims

OUT_DIR = ("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/15_transcriptionFactorMotifs/dataDriven_commonBaseline")
# dir.create(OUT_DIR)
setwd(OUT_DIR)

write.table(homer_res_list$`bg-pbs`, paste0("annotateInput_UP_bg.txt"), 
              sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`5f-pbs`, paste0("annotateInput_UP_5f.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`fen-pbs`, paste0("annotateInput_UP_fen.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`fluni-pbs`, paste0("annotateInput_UP_fluni.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`hc-pbs`, paste0("annotateInput_UP_hc.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`hq-pbs`, paste0("annotateInput_UP_hq.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`myr-pbs`, paste0("annotateInput_UP_myr.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
write.table(homer_res_list$`nerol-pbs`, paste0("annotateInput_UP_nerol.txt"), 
            sep = '\t', quote = FALSE, row.names=FALSE)
