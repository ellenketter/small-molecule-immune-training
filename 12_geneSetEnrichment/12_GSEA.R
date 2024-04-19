library(fgsea)
library(data.table)
library(ggplot2)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/12_geneSetEnrichment/")
hallmarks_symbols <- fgsea::gmtPathways(gmt.file = "h.all.v7.2.symbols.gmt")

limma_res_files <- list.files(path = "../../demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/",
                              pattern = "^results",
                              full.names = T)

names(limma_res_files) <- gsub(pattern = "^results_",replacement = "",
                               gsub(pattern =".txt$",replacement="",
                                    basename(limma_res_files)))


stims <- names(limma_res_files)
res_list <- lapply(limma_res_files, function(x){
  read.delim(x)
 })

res_list$`5f`$PeakID <- rownames(res_list$`5f`)
res_list$bg$PeakID <- rownames(res_list$bg)
res_list$fen$PeakID <- rownames(res_list$fen)
res_list$fluni$PeakID <- rownames(res_list$fluni)
res_list$hc$PeakID <- rownames(res_list$hc)
res_list$hq$PeakID <- rownames(res_list$hq)
res_list$myr$PeakID <- rownames(res_list$myr)
res_list$nerol$PeakID <- rownames(res_list$nerol)

annotate <- read.delim("../../demultiplexed/FastQ/11_annotatePeaks/old/background_annotate_out.txt")
annotate_names <- colnames(annotate)
annotate_names[1] <- "PeakID"
colnames(annotate) <- annotate_names
annotate <- annotate[c("PeakID","Gene.Name")]
annotate <- annotate %>% arrange(as.numeric(PeakID))




# merge Gene.Name to each df in a list

mergeGeneNames <- function(x){merge(annotate,x, by="PeakID") %>%
    arrange(as.numeric(PeakID))}
limma_res_list <- lapply(res_list, mergeGeneNames)
# saveRDS(limma_res_list, "limma_res_list.rds")
 
# daa is a list of limma diff expression results, 
# cts are the names of daa for which you want to run fgsea 
# name - name of the column in the limma results that points to the symbol gene name 
# value - column that will be used to rank the genes for gsea
# p_name - if there's a need to add pvalue then will perform as.numeric(ict_gene_list[,value]) * -log10(abs(ict_gene_list[,p_name]))
# GO - FALSE (default) will run fgsea for the hallmark pathways - GO=TRUE will run gsea for gene ontologies you need org.Hs.eg.db activated and library(org.Hs.eg.db) 
## EXAMPLE
# GSEA_run(daa=limma_res_list, cts=c("BG","stim2")) 

GSEA_run <- function(daa, cts=NULL, GO=FALSE, value="logFC", p_name="adj.P.Val", name="Gene.Name"){
  
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


# limma_res_list <- list(bg=bg, fiveF=fiveF, fen=fen, fluni=fluni, hc=hc,
#                        hq=hq, myr=myr, nerol=nerol)
# test first just element 1 and 2 from the list 
test <- GSEA_run(daa=limma_res_list, cts = names(limma_res_list), name = "")

saveRDS(test, "fgseaResult.rds")

