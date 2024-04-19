library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/otherModels/dataDriven_commonBaseline/")
hallmarkDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/12_geneSetEnrichment/"
# annotateInputDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/"


# hallmarks_symbols <- fgsea::gmtPathways(gmt.file = paste0(hallmarkDir,"h.all.v7.2.symbols.gmt"))

hallmarks_symbols <- gmtPathways(system.file(
  "extdata", "mouse.reactome.gmt", package="fgsea"))



res_list <- readRDS("mashMeans_dd.rds")

res_list$lfsr[which(res_list$lfsr == 0)] <- min(res_list$lfsr[res_list$lfsr > 0]) # prevent error with infinite values
res_list$product <- res_list$posteriorMean*-log10(res_list$lfsr)

stims <- levels(res_list$treatment)

# re-ran homer with concatenated peak info as peakid
# annotateInput <- read.delim("../../10_differentialAccessibility/background_vs_DA_peaks/backgroundPeaks.txt")
# annotateInput$Geneid <- paste0(annotateInput$Chr,".",annotateInput$Start,".",annotateInput$End)

res_list <- res_list %>% separate(peak, c('chromosome', 'start', 'end'))
res_list$peak <- paste0(res_list$chromosome,".", res_list$start,".", res_list$end)

annotate <- read.delim("../../annotateOutput.txt")
annotate$PeakID <- paste0(annotate$Chr,".",annotate$Start,".",annotate$End)

annotate_names <- colnames(annotate)
annotate_names[1] <- "peak"
colnames(annotate) <- annotate_names
# check that peakids match peak dimensions
# note that peak start is +1bp from annotate$PeakID!!!! 

# remove extra columns
annotate <- annotate[c("peak","Entrez.ID")]

# remove row 59925, which contains colname "geneid"
# annotate <-  annotate[-which(annotate$peak == "Entrez.ID"),]

# merge Gene.Name to each df in a list
# if peak is in column, assign gene.name 
mash_res_list <- merge(res_list, annotate)
mash_res_list <- arrange(mash_res_list, desc(product))
mash_res_list <- split(mash_res_list, f = mash_res_list$treatment)

mash_res_list$bg <- mash_res_list$bg[!duplicated(mash_res_list$bg$Entrez.ID),]
mash_res_list$'5f' <- mash_res_list$'5f'[!duplicated(mash_res_list$'5f'$Entrez.ID),]
mash_res_list$fen <- mash_res_list$fen[!duplicated(mash_res_list$fen$Entrez.ID),]
mash_res_list$fluni <- mash_res_list$fluni[!duplicated(mash_res_list$fluni$Entrez.ID),]
mash_res_list$hc <- mash_res_list$hc[!duplicated(mash_res_list$hc$Entrez.ID),]
mash_res_list$hq <- mash_res_list$hq[!duplicated(mash_res_list$hq$Entrez.ID),]
mash_res_list$myr <- mash_res_list$myr[!duplicated(mash_res_list$myr$Entrez.ID),]
mash_res_list$nerol <- mash_res_list$nerol[!duplicated(mash_res_list$nerol$Entrez.ID),]

saveRDS(mash_res_list, "mash_res_list_reactome.rds")
# mash_res_list <- readRDS("mash_res_list_reactome.rds")

ranks <- c()
ranks$bg <- pull(mash_res_list$bg,product)
names(ranks$bg) <- mash_res_list$bg$Entrez.ID

ranks$'5f' <- pull(mash_res_list$'5f',product)
names(ranks$'5f') <- mash_res_list$'5f'$Entrez.ID

ranks$fen <- pull(mash_res_list$fen,product)
names(ranks$fen) <- mash_res_list$fen$Entrez.ID

ranks$fluni <- pull(mash_res_list$fluni,product)
names(ranks$fluni) <- mash_res_list$fluni$Entrez.ID

ranks$hc <- pull(mash_res_list$hc,product)
names(ranks$hc) <- mash_res_list$hc$Entrez.ID

ranks$hq <- pull(mash_res_list$hq,product)
names(ranks$hq) <- mash_res_list$hq$Entrez.ID

ranks$myr <- pull(mash_res_list$myr,product)
names(ranks$myr) <- mash_res_list$myr$Entrez.ID

ranks$nerol <- pull(mash_res_list$nerol,product)
names(ranks$nerol) <- mash_res_list$nerol$Entrez.ID
 
# Run gene set enrichment 
fgseaRes <- c()
fgseaRes$bg <- fgsea(pathways = hallmarks_symbols, 
                  stats    = ranks$bg,
                  minSize  = 15,
                  maxSize  = 500,
                  nPermSimple = 10000)
fgseaRes$'5f' <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$'5f',
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)
fgseaRes$fen <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$fen,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)
fgseaRes$fluni <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$fluni,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)
fgseaRes$hc <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$hc,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)

fgseaRes$hq <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$hq,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)
fgseaRes$myr <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$myr,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)
fgseaRes$nerol <- fgsea(pathways = hallmarks_symbols, 
                     stats    = ranks$nerol,
                     minSize  = 15,
                     maxSize  = 500,
                     nPermSimple = 10000)

saveRDS(fgseaRes, "fgseaResult_reactome.rds")

# Filter results to include top several from each treatment

bgTop <- head(arrange(fgseaRes$bg, by = padj), 10)
X5fTop <- head(arrange(fgseaRes$`5f`, by = padj), 10)
fenTop <- head(arrange(fgseaRes$fen, by = padj), 10)
fluniTop <- head(arrange(fgseaRes$fluni, by = padj), 10)
hcTop <- head(arrange(fgseaRes$hc, by = padj), 10)
hqTop <- head(arrange(fgseaRes$hq, by = padj), 10)
myrTop <- head(arrange(fgseaRes$myr, by = padj), 10)
nerolTop <- head(arrange(fgseaRes$nerol, by = padj), 10)

sigPathways <- unique(c(bgTop$pathway, X5fTop$pathway, fenTop$pathway, 
                      fluniTop$pathway, hcTop$pathway, hqTop$pathway, 
                      myrTop$pathway, nerolTop$pathway))

fgseaRes$bg <- fgseaRes$bg[fgseaRes$bg$pathway %in% sigPathways]
fgseaRes$`5f` <- fgseaRes$`5f`[fgseaRes$`5f`$pathway %in% sigPathways]
fgseaRes$fen <- fgseaRes$fen[fgseaRes$fen$pathway %in% sigPathways]
fgseaRes$fluni <- fgseaRes$fluni[fgseaRes$fluni$pathway %in% sigPathways]
fgseaRes$hc <- fgseaRes$hc[fgseaRes$hc$pathway %in% sigPathways]
fgseaRes$hq <- fgseaRes$hq[fgseaRes$hq$pathway %in% sigPathways]
fgseaRes$myr <- fgseaRes$myr[fgseaRes$myr$pathway %in% sigPathways]
fgseaRes$nerol <- fgseaRes$nerol[fgseaRes$nerol$pathway %in% sigPathways]

saveRDS(fgseaRes, "fgseaResult_reactome_top10perTreatment.rds")


# QC ---------------


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

dev.off()
