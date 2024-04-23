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

OUT_DIR <- "fgseaTwoList/treatmentOverlaps/vennDiagram/"
# dir.create(OUT_DIR)

m <- readRDS("mash_res_list.rds")

universe <- m$bg$Gene.Name

geneNames <- list.files("upsetPlots/vennDiagram/" ,pattern = "\\.csv$")

intersections <- c()  
importIntersections <- function(x){
  new <- read.csv(paste0("upsetPlots/vennDiagram/", x))
  intersections <- new$Gene.Name
}

gcUnique <- lapply (geneNames, importIntersections)

finalNames <- gsub(".csv", "", geneNames)
finalNames <- gsub("geneNames_","",finalNames)
names(gcUnique) <- finalNames

# make composite for steroids
steroidAll <- c(gcUnique$bgSteroids_3, gcUnique$steroidsOnly_4, 
                      gcUnique$steroidsOtherDrugs_7, gcUnique$all_5)
gcUnique$steroidAll <- steroidAll

# make composite for otherDrugs
otherDrugsAll <- c(gcUnique$bgOtherDrugs_2, gcUnique$otherDrugsOnly_6, 
                gcUnique$steroidsOtherDrugs_7, gcUnique$all_5)
gcUnique$otherDrugsAll <- otherDrugsAll
                               
       

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

saveRDS(fgsea_ora,"fgseaTwoList/treatmentOverlaps/vennDiagram/treatmentOverlaps1.rds")

top5 <- lapply(fgsea_ora, function(x){head(arrange(x, by=padj), 5)})
matrixTop5 <- bind_rows(top5, .id = "overlap")

sigPathways <- unique(matrixTop5a$gs_description)

union <- lapply(fgsea_ora, function(x){x[x$gs_description %in% sigPathways]})
# remove rows where no condition is <0.1


saveRDS(union, paste0(OUT_DIR, "fgsea_ora_reactome_top5perTreatment_venn_0.1.rds"))
