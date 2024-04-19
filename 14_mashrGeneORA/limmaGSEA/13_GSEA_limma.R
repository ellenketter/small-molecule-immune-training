library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/")
hallmarkDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/12_geneSetEnrichment/"
# annotateInputDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/"


hallmarks_symbols <- fgsea::gmtPathways(gmt.file = paste0(hallmarkDir,"mh.all.v2023.1.Mm.symbols.gmt"))

stims <- c("bg", "5f", "fen", "fluni", "hc", "hq", "myr", "nerol")
openFiles<- function (x){read.table(paste0("background_vs_DA_peaks/results_", x, ".txt")) %>% 
    subset(select = c('genes', 't'))}
limmaResults <- lapply(stims, openFiles)
names(limmaResults) <- stims

limmaResults <- reshape2::melt(limmaResults)
colnames(limmaResults) <- c("peak", "stat", "tStat", "treatment")
limmaResults <- subset(limmaResults, select = c("peak", "tStat", "treatment"))

# re-ran homer with concatenated peak info as peakid
# annotateInput <- read.delim("../../10_differentialAccessibility/background_vs_DA_peaks/backgroundPeaks.txt")
# annotateInput$Geneid <- paste0(annotateInput$Chr,".",annotateInput$Start,".",annotateInput$End)

res_list <- limmaResults %>% separate(peak, c('chromosome', 'start', 'end'))
res_list <- cbind(res_list, paste0(res_list$chromosome,".", res_list$start,".", res_list$end))
colnames(res_list)[6] <- "peak"



annotate <- read.delim("../13_mashr/annotateOutput.txt")
annotate$PeakID <- paste0(annotate$Chr,".",annotate$Start,".",annotate$End)

annotate_names <- colnames(annotate)
annotate_names[1] <- "peak"
colnames(annotate) <- annotate_names
# check that peakids match peak dimensions
# note that peak start is +1bp from annotate$PeakID!!!! 

# remove extra columns
annotate <- annotate[c("peak","Gene.Name")]

# remove row 59925, which contains colname "geneid"
annotate <-  annotate[-which(annotate$peak == "Geneid"),]

# merge Gene.Name to each df in a list
# if peak is in column, assign gene.name 
mash_res_list <- merge(res_list, annotate)
mash_res_list <- arrange(mash_res_list, desc(tStat))
mash_res_list <- split(mash_res_list, f = mash_res_list$treatment)

mash_res_list$bg <- mash_res_list$bg[!duplicated(mash_res_list$bg$Gene.Name),]
mash_res_list$'5f' <- mash_res_list$'5f'[!duplicated(mash_res_list$'5f'$Gene.Name),]
mash_res_list$fen <- mash_res_list$fen[!duplicated(mash_res_list$fen$Gene.Name),]
mash_res_list$fluni <- mash_res_list$fluni[!duplicated(mash_res_list$fluni$Gene.Name),]
mash_res_list$hc <- mash_res_list$hc[!duplicated(mash_res_list$hc$Gene.Name),]
mash_res_list$hq <- mash_res_list$hq[!duplicated(mash_res_list$hq$Gene.Name),]
mash_res_list$myr <- mash_res_list$myr[!duplicated(mash_res_list$myr$Gene.Name),]
mash_res_list$nerol <- mash_res_list$nerol[!duplicated(mash_res_list$nerol$Gene.Name),]

# saveRDS(mash_res_list, "tStat_annotated.rds")

ranks <- c()
ranks$bg <- pull(mash_res_list$bg,tStat)
names(ranks$bg) <- mash_res_list$bg$Gene.Name

ranks$'5f' <- pull(mash_res_list$'5f',tStat)
names(ranks$'5f') <- mash_res_list$'5f'$Gene.Name

ranks$fen <- pull(mash_res_list$fen,tStat)
names(ranks$fen) <- mash_res_list$fen$Gene.Name

ranks$fluni <- pull(mash_res_list$fluni,tStat)
names(ranks$fluni) <- mash_res_list$fluni$Gene.Name

ranks$hc <- pull(mash_res_list$hc,tStat)
names(ranks$hc) <- mash_res_list$hc$Gene.Name

ranks$hq <- pull(mash_res_list$hq,tStat)
names(ranks$hq) <- mash_res_list$hq$Gene.Name

ranks$myr <- pull(mash_res_list$myr,tStat)
names(ranks$myr) <- mash_res_list$myr$Gene.Name

ranks$nerol <- pull(mash_res_list$nerol,tStat)
names(ranks$nerol) <- mash_res_list$nerol$Gene.Name
 
# Run gene set enrichment 
fgseaRes <- c()
fgseaRes$bg <- fgsea(pathways = hallmarks_symbols, 
                  stats    = ranks$bg,
                  minSize  = 15,
                  maxSize  = 500)
fgseaRes$'5f' <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$'5f',
                     minSize  = 15,
                     maxSize  = 500)
fgseaRes$fen <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$fen,
                     minSize  = 15,
                     maxSize  = 500)
fgseaRes$fluni <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$fluni,
                     minSize  = 15,
                     maxSize  = 500)
fgseaRes$hc <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$hc,
                     minSize  = 15,
                     maxSize  = 500)

fgseaRes$hq <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$hq,
                     minSize  = 15,
                     maxSize  = 500)
fgseaRes$myr <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$myr,
                     minSize  = 15,
                     maxSize  = 500)
fgseaRes$nerol <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$nerol,
                     minSize  = 15,
                     maxSize  = 500)

saveRDS(fgseaRes, "fgseaResult.rds")



ggplot(fgseaRes$bg, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Beta Glucan")

ggplot(fgseaRes$'5f', aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: 5F")

ggplot(fgseaRes$fen, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Fenoterol")

ggplot(fgseaRes$fluni, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Flunisolide")

ggplot(fgseaRes$hc, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Hydrocortisone")

ggplot(fgseaRes$myr, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Myricetin")

ggplot(fgseaRes$nerol, aes(reorder(pathway, NES), NES)) +
  geom_col(aes(fill=padj<0.05)) +
  coord_flip() +
  labs(x="Pathway", y="Normalized Enrichment Score",
       title="Hallmark GSEA: Nerol")
