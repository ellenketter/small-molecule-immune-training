library(fgsea)
library(data.table)
library(ggplot2)
library(dplyr)
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
hallmarkDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/12_geneSetEnrichment/"
# annotateInputDir <- "/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/"


hallmarks_symbols <- fgsea::gmtPathways(gmt.file = paste0(hallmarkDir,"mh.all.v2023.1.Mm.symbols.gmt"))



res_list <- readRDS("mashMeans_dd.rds")
res_list$lfsr[which(res_list$lfsr == 0)] <- 1e-17 # prevent error with infinite values
res_list$product <- res_list$posteriorMean*-log10(res_list$lfsr)

stims <- unique(res_list$treatment)
mashPriorsAndPosteriors$peak == outliers[i]
# re-run homer with concatenated peak info as peakid
# annotateInput <- read.delim("../../10_differentialAccessibility/background_vs_DA_peaks/backgroundPeaks.txt")
# annotateInput$Geneid <- paste0(annotateInput$Chr,".",annotateInput$Start,".",annotateInput$End)
# write.table(annotateInput, "../annotateInput.txt", sep = '\t', quote = FALSE, row.names=FALSE)

annotate <- read.delim("../annotateOutput.txt")
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
mash_res_list <- arrange(mash_res_list, desc(product))
mash_res_list <-  mash_res_list[-which(mash_res_list$treatment == "bg"),]

mash_res_list <- split(mash_res_list, f = mash_res_list$treatment)

saveRDS(mash_res_list, "mash_res_list.rds")
 
# daa is a list of limma diff expression results, 
# cts are the names of daa for which you want to run fgsea 
# name - name of the column in the limma results that points to the symbol gene name 
# value - column that will be used to rank the genes for gsea
# p_name - if there's a need to add pvalue then will perform as.numeric(ict_gene_list[,value]) * -log10(abs(ict_gene_list[,p_name]))
# GO - FALSE (default) will run fgsea for the hallmark pathways - GO=TRUE will run gsea for gene ontologies you need org.Hs.eg.db activated and library(org.Hs.eg.db) 
## EXAMPLE
# GSEA_run(daa=limma_res_list, cts=c("BG","stim2")) 

GSEA_run <- function(daa, cts=NULL, GO=FALSE, value="posteriorMean", p_name="lfsr", name="Gene.Name"){
  
  if(!is.null(cts)){
    if(!all(names(daa) %in% cts)){
      stop("some cell types are not present in the DAA summary stats")
    }
    daa <- daa[cts]
  }
  res_list <- future.apply::future_lapply(cts, function(x){
    
    ict_res <- as.data.frame(daa[[x]])
    if(is.null(p_name)){
      ict_gene_v <- as.numeric(ict_res[,value])
    } else {
      ict_gene_v <- as.numeric(ict_res[,value]) * -log10(abs(ict_res[,p_name]))
    }
    names(ict_gene_v) <- ict_res[,name] 
    
    ict_gene_v <- ict_gene_v[!is.na(ict_gene_v)]
    ict_gene_v <- sort(ict_gene_v, decreasing = T)
    
    if(GO){
      # dont run GO unless there are no significant hits in hallmarks. 
      ict_gsea <- gseGO(geneList = ict_gene_v, 
                        keyType = "SYMBOL",
                        OrgDb        = org.Hs.eg.db,
                        ont          = "BP",
                        pvalueCutoff = 1,
                        maxGSSize = 600, 
                        minGSSize = 80)
      
    } else{
      
      fgseaRes <- fgsea::fgsea(pathways = hallmarks_symbols, 
                               stats= ict_gene_v,
                               maxSize  = 500,
                               nproc=1)
      
      ict_gsea <- data.frame(
        pathways = fgseaRes$pathway,
        pval = fgseaRes$pval,
        NES = fgseaRes$NES,
        padj = fgseaRes$padj
      )
    }
    
    ict_gsea
  })
  names(res_list) <- cts
  return(res_list)
}


mash_res_list <- list(fiveF=mash_res_list$'5f', 
                      fen=mash_res_list$fen, fluni=mash_res_list$fluni, 
                      hc=mash_res_list$hc, hq=mash_res_list$hq, 
                      myr=mash_res_list$myr, nerol=mash_res_list$nerol)
# test first just element 1 and 2 from the list 
test <- GSEA_run(daa=mash_res_list, cts = names(mash_res_list))

saveRDS(test, "fgseaResult.rds")

