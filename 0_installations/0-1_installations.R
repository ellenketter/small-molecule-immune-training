# Install R packages ###################
# Requisite for bulk ATAC-seq analysis in Knight et al. PNAS. (2024)
# R version 4.1.0 
# on RCC linux cluster - module load gcc gsl
# 
# Ellen Ketter
# April 2024
########################################

install.packages(c("tidyverse", "stringr", "ggplot2", "hrbrthemes", "viridis",
                   "pheatmap", "wesanderson", "statmod", "goseq", "RColorBrewer",
                   "ggfortify", "devtools", "FactoMineR", "factoextra", 
                   "cluster", "ggrepel", "data.table", "VennDiagram", "UpSetR",
                   "readr"))


# install BiocManager
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
library(BiocManager)

# Install BiocManager packages
BiocManager::install(c("ChIPseeker", "ATACseqQC", "ChIPpeakAnno", "MotifDb", 
                       "GenomicAlignments","BSgenome.Mmusculus.UCSC.mm10", 
                       "TxDb.Mmusculus.UCSC.mm10.knownGene","Rsamtools", 
                       "DESeq2", "goseq", "edgeR"))


