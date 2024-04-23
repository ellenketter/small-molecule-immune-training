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

OUT_DIR <- "fgseaTwoList/"
#dir.create(OUT_DIR)

m <- readRDS("mash_res_list.rds")

universe <- m$bg$Gene.Name
  

gcUnique <- c()
gcUnique$bg <- m$bg$Gene.Name[m$bg$significant == "UP"]
gcUnique$`5f` <- m$`5f`$Gene.Name[m$`5f`$significant == "UP"]
gcUnique$fen <- m$fen$Gene.Name[m$fen$significant == "UP"]
gcUnique$fluni <- m$fluni$Gene.Name[m$fluni$significant == "UP"]
gcUnique$hc <- m$hc$Gene.Name[m$hc$significant == "UP"]
gcUnique$hq <- m$hq$Gene.Name[m$hq$significant == "UP"]
gcUnique$myr <- m$myr$Gene.Name[m$myr$significant == "UP"]
gcUnique$nerol <- m$nerol$Gene.Name[m$nerol$significant == "UP"]




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

saveRDS(fgsea_ora,"fgseaTwoList/CP_reactome.rds")

bgTop <- head(arrange(fgsea_ora$bg, by = padj), 10)
X5fTop <- head(arrange(fgsea_ora$`5f`, by = padj), 10)
fenTop <- head(arrange(fgsea_ora$fen, by = padj), 10)
fluniTop <- head(arrange(fgsea_ora$fluni, by = padj), 10)
hcTop <- head(arrange(fgsea_ora$hc, by = padj), 10)
hqTop <- head(arrange(fgsea_ora$hq, by = padj), 10)
myrTop <- head(arrange(fgsea_ora$myr, by = padj), 10)
nerolTop <- head(arrange(fgsea_ora$nerol, by = padj), 10)

sigPathways <- unique(c(bgTop$gs_description, X5fTop$gs_description, fenTop$gs_description, 
                        fluniTop$gs_description, hcTop$gs_description, 
                        hqTop$gs_description, 
                        myrTop$gs_description, nerolTop$gs_description))

fgsea_ora$bg <- fgsea_ora$bg[fgsea_ora$bg$gs_description %in% sigPathways]
fgsea_ora$`5f` <- fgsea_ora$`5f`[fgsea_ora$`5f`$gs_description %in% sigPathways]
fgsea_ora$fen <- fgsea_ora$fen[fgsea_ora$fen$gs_description %in% sigPathways]
fgsea_ora$fluni <- fgsea_ora$fluni[fgsea_ora$fluni$gs_description %in% sigPathways]
fgsea_ora$hc <- fgsea_ora$hc[fgsea_ora$hc$gs_description %in% sigPathways]
fgsea_ora$hq <- fgsea_ora$hq[fgsea_ora$hq$gs_description %in% sigPathways]
fgsea_ora$myr <- fgsea_ora$myr[fgsea_ora$myr$gs_description %in% sigPathways]
fgsea_ora$nerol <- fgsea_ora$nerol[fgsea_ora$nerol$gs_description %in% sigPathways]

saveRDS(fgsea_ora, "fgsea_ora_reactome_top10perTreatment_new.rds")
