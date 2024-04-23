# code from vignette:
# https://pnnl-comp-mass-spec.github.io/proteomics-data-analysis-tutorial/ora.html
# pathways overview
# https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp

library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)



setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")

OUT_DIR <- "fgseaTwoList/treatmentOverlaps/withPositiveControl/"
# dir.create(OUT_DIR)

m <- readRDS("mash_res_list.rds")

universe <- m$bg$Gene.Name

geneNames <- list.files("upsetPlots/" ,pattern = "\\.csv$")

intersections <- c()  
importIntersections <- function(x){
  new <- read.csv(paste0("upsetPlots/", x))
  intersections <- new$Gene.Name
}

gcUnique <- lapply (geneNames, importIntersections)

finalNames <- gsub(".csv", "", geneNames)
finalNames <- gsub("geneNames_","",finalNames)
names(gcUnique) <- finalNames

# make composite for steroids
hcFluni_hc_fluni <- c(gcUnique$hcFluni_5, gcUnique$fluniOnly_10, gcUnique$hcOnly_12)
gcUnique$hcFluni_hc_fluni <- hcFluni_hc_fluni


# reoroder gcUnique
gcUnique <- gcUnique[c("bgAlone_1", "no_hcFluni_2","myrBg_3" , "all_4", 
                       "hcFluni_5", "hcFluniBg_6", "fluniBg_7","no_hcFluniBg_8", 
                       "noHC_9", "fluniOnly_10", "myrOnly_11", "hcOnly_12", "noBg_13",
                       "hcFluni_hc_fluni", "allUpPeaks",
                       "noHcFluni5f_14", "fenBg_15", "hcBg_16")]
                               
       

library(msigdbr)
msigdbr::msigdbr_collections() # available collections
BP_db <- msigdbr(species = "Mus musculus", 
                 category = "C2", subcategory = "CP:REACTOME")
head(BP_db)


# Convert to a list of gene sets
BP_conv <- unique(BP_db[, c("gene_symbol", "gs_exact_source")])
BP_list <- split(x = BP_conv$gene_symbol, f = BP_conv$gs_exact_source)
# First ~6 IDs of first 3 terms
lapply(head(BP_list, 3), head)

fgsea_ora <- lapply(seq_along(gcUnique), function(i) {
  fora(pathways = BP_list, 
       genes = gcUnique[[i]], # genes in cluster i
       universe = universe, # all genes
       minSize = 15, 
       maxSize = 500) %>% 
    mutate(cluster = names(gcUnique)[i]) # add cluster column
}) %>% 
  data.table::rbindlist() %>% # combine tables
  # filter(padj < 0.2) %>% 
  arrange(cluster, padj) %>% 
  # Add additional columns from BP_db
  left_join(distinct(BP_db, gs_subcat, gs_exact_source, 
                     gs_name, gs_description),
            by = c("pathway" = "gs_exact_source")) %>% 
  # Reformat descriptions
  mutate(gs_name = sub("^GOBP_", "", gs_name),
         gs_name = gsub("_", " ", gs_name))

# First 6 rows
head(fgsea_ora)
fgsea_ora$percentOverlap <- (fgsea_ora$overlap/fgsea_ora$size)*100
fgsea_ora <- split(fgsea_ora, f = fgsea_ora$cluster)

saveRDS(fgsea_ora,"fgseaTwoList/treatmentOverlaps/withPositiveControl/treatmentOverlaps1.rds")

bgAlone_1 <- head(arrange(fgsea_ora$bgAlone_1, by = padj), 5)
no_hcFluni_2 <- head(arrange(fgsea_ora$no_hcFluni_2, by = padj), 5)
myrBg_3 <- head(arrange(fgsea_ora$myrBg_3, by = padj), 5)
all_4 <- head(arrange(fgsea_ora$all_4, by = padj), 5)
hcFluni_5 <- head(arrange(fgsea_ora$hcFluni_5, by = padj), 5)
hcFluniBg_6 <- head(arrange(fgsea_ora$hcFluniBg_6, by = padj), 5)
fluniBg_7 <- head(arrange(fgsea_ora$fluniBg_7, by = padj), 5)
no_hcFluniBg_8<- head(arrange(fgsea_ora$no_hcFluniBg_8, by = padj), 5)
noHC_9 <- head(arrange(fgsea_ora$noHC_9, by = padj), 5)
fluniOnly_10 <- head(arrange(fgsea_ora$fluniOnly_10, by = padj), 5)
myrOnly_11 <- head(arrange(fgsea_ora$myrOnly_11, by = padj), 5)
hcOnly_12 <- head(arrange(fgsea_ora$hcOnly_12, by = padj), 5)
noBg_13 <- head(arrange(fgsea_ora$noBg_13, by = padj), 5)
hcFluni_hc_fluni <- head(arrange(fgsea_ora$hcFluni_hc_fluni, by = padj), 5)
allBgPeaksUp <- head(arrange(fgsea_ora$allBgPeaksUp, by = padj), 5)

noHcFluni5f_14 <- head(arrange(fgsea_ora$noHcFluni5f, by = padj), 5)
fenBg_15 <- head(arrange(fgsea_ora$fenBg_15, by = padj), 5)
hcBg_16 <- head(arrange(fgsea_ora$hcBg_16, by = padj), 5)

sigPathways <- unique(c(bgAlone_1$gs_description, no_hcFluni_2$gs_description,
                 myrBg_3$gs_description, all_4$gs_description, 
                 hcFluni_5$gs_description, hcFluniBg_6$gs_description, 
                 fluniBg_7$gs_description,no_hcFluniBg_8$gs_description, 
                 noHC_9$gs_description, fluniOnly_10$gs_description, 
                 myrOnly_11$gs_description, hcOnly_12$gs_description, 
                 noBg_13$gs_description, hcFluni_hc_fluni$gs_description,allBgPeaksUp$gs_description))
                 # noHcFluni5f_14$gs_description, 
                 # fenBg_15$gs_description, hcBg_16$gs_description))



fgsea_ora$bgAlone_1 <- fgsea_ora$bgAlone_1[fgsea_ora$bgAlone_1$gs_description %in%
                                             sigPathways]
fgsea_ora$no_hcFluni_2 <- fgsea_ora$no_hcFluni_2[fgsea_ora$no_hcFluni_2$gs_description %in%
                                             sigPathways]
fgsea_ora$myrBg_3 <- fgsea_ora$myrBg_3[fgsea_ora$myrBg_3$gs_description %in%
                                             sigPathways]
fgsea_ora$all_4 <- fgsea_ora$all_4[fgsea_ora$all_4$gs_description %in%
                                             sigPathways]
fgsea_ora$hcFluni_5 <- fgsea_ora$hcFluni_5[fgsea_ora$hcFluni_5$gs_description %in%
                                             sigPathways]
fgsea_ora$hcFluniBg_6 <- fgsea_ora$hcFluniBg_6[fgsea_ora$hcFluniBg_6$gs_description %in%
                                             sigPathways]
fgsea_ora$fluniBg_7 <- fgsea_ora$fluniBg_7[fgsea_ora$fluniBg_7$gs_description %in%
                                             sigPathways]
fgsea_ora$no_hcFluniBg_8 <- fgsea_ora$no_hcFluniBg_8[fgsea_ora$no_hcFluniBg_8$gs_description %in%
                                             sigPathways]
fgsea_ora$noHC_9 <- fgsea_ora$noHC_9[fgsea_ora$noHC_9$gs_description %in%
                                             sigPathways]
fgsea_ora$fluniOnly_10 <- fgsea_ora$fluniOnly_10[fgsea_ora$fluniOnly_10$gs_description %in% # 
                                             sigPathways]
fgsea_ora$myrOnly_11 <- fgsea_ora$myrOnly_11[fgsea_ora$myrOnly_11$gs_description %in%
                                             sigPathways]
fgsea_ora$hcOnly_12 <- fgsea_ora$hcOnly_12[fgsea_ora$hcOnly_12$gs_description %in% # 
                                             sigPathways]
fgsea_ora$noBg_13 <- fgsea_ora$noBg_13[fgsea_ora$noBg_13$gs_description %in%
                                             sigPathways]
fgsea_ora$hcFluni_hc_fluni<- fgsea_ora$hcFluni_hc_fluni[fgsea_ora$hcFluni_hc_fluni$gs_description %in%
                                         sigPathways]

fgsea_ora$noHcFluni5f_14 <- fgsea_ora$noHcFluni5f_14[fgsea_ora$noHcFluni5f_14$gs_description %in%
                                                sigPathways]
fgsea_ora$fenBg_15 <- fgsea_ora$fenBg_15[fgsea_ora$fenBg_15$gs_description %in%
                                          sigPathways]
fgsea_ora$hcBg_16 <- fgsea_ora$hcBg_16[fgsea_ora$hcBg_16$gs_description %in%
                                         sigPathways]
allBgPeaksUp <- fgsea_ora$allBgPeaksUp[fgsea_ora$allBgPeaksUp$gs_description %in%
                                         sigPathways]


saveRDS(fgsea_ora, paste0(OUT_DIR, "fgsea_ora_reactome_top5perTreatment_allBgUp.rds"))

