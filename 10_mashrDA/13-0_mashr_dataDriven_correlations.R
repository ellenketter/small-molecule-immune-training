library(ashr)
library(mashr)
library(gplots)
library(viridis)
library(ggplot2)
library(corrplot)
library(rmeta)
library(dplyr)
library(reshape2)
library("heatmaply")
set.seed(2020)

# could also use DA genes from limma as 'strong' input

current = "/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/"
setwd(current)
folder = "dataDriven_withCorrelations/"

# Make input matrix of Bhat (effect size = fit$coefficients) and Shat 
# (standard error = fit$stdev.unscaled * fit$sigma), see p. 76
# https://bioconductor.riken.jp/packages/3.1/bioc/manuals/limma/man/limma.pdf

betaMatrix <- readRDS( "betaMatrix.rds")

## Run mash to define canonical covariance matrix ---------------
data = mash_set_data(betaMatrix$Bhat, betaMatrix$Shat)

# estimate with and without correlations
V.simple = estimate_null_correlation_simple(data)
data.Vsimple = mash_update_data(data, V=V.simple)

m.1by1 = mash_1by1(data.Vsimple) # takes a minute 
strong = get_significant_results(m.1by1,0.05)

U.pca = cov_pca(data.Vsimple,5,subset=strong)
U.ed = cov_ed(data.Vsimple, U.pca, subset=strong) # apply extreme deconvolution, ~5min

m = mash(data.Vsimple, U.ed)
print(get_loglik(m),digits = 10) # fit the model

saveRDS(m, paste0(folder,"mashResults_dd.rds"))
m <- readRDS(paste0(folder,"mashResults_dd.rds"))
## Processing -------------
localFalseSignRates <- melt(get_lfsr(m))
colnames(localFalseSignRates) <- c("peak","treatment","lfsr")
posteriorMeans <- melt(get_pm(m))
colnames(posteriorMeans) <- c("peak","treatment","posteriorMean")
mashMeans <- inner_join(posteriorMeans, localFalseSignRates, 
                        by = c("peak","treatment"))
mashMeans$significant <- "NO"
mashMeans$significant[mashMeans$posteriorMean > 0 & mashMeans$lfsr < 0.05] <- "UP"
mashMeans$significant[mashMeans$posteriorMean < 0 & mashMeans$lfsr < 0.05] <- "DOWN"

saveRDS(mashMeans, paste0(folder,"mashMeans_dd.rds"))

## Get all significant results ----------
significantMASH <- get_significant_results(m)


## Export significant results list ----------
treatments <- colnames(betaMatrix$Bhat)

getSigResults <- function(x){get_significant_results(m, conditions=x)}
mashSigResults <- lapply(treatments, getSigResults)

names(mashSigResults) <- treatments

saveRDS(mashSigResults, paste0(folder,"mashSigResults_dd.rds"))

## Visualize pairwise sharing -------------------
sharedMatrix <- get_pairwise_sharing(m)


heatmaply(sharedMatrix, 
          #dendrogram = "row",
          xlab = "", ylab = "", 
          main = "",
          scale = "column",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          hide_colorbar = TRUE,
          branches_lwd = 0.1,
          fontsize_row = 10, fontsize_col = 10,
          labCol = colnames(sharedMatrix),
          labRow = rownames(sharedMatrix),
          heatmap_layers = theme(axis.line=element_blank()),
          file = "dataDriven_withCorrelations/covarianceMatrixHeatmap_dd.svg")


heatmaply(sharedMatrix, file = paste0(folder,"covarianceMatrixHeatmap_dd.html"))

# final version for figure 
fullNames <- c("5-Fluoroindole 2-Carboxylic Acid", "Beta Glucan", "Fenoterol", 
               "Flunisolide", "Hydrocortisone", "Hydroquinone", "Myricetin", "Nerol")
  
heatmaply(sharedMatrix, 
          #dendrogram = "row",
          xlab = "", ylab = "", 
          main = "",
          margins = c(60,100,40,20),
          grid_color = "white",
          grid_width = 0.00001,
          titleX = FALSE,
          branches_lwd = 0.1,
          fontsize_row = 10, fontsize_col = 10,
          labCol = fullNames,
          labRow = fullNames,
          heatmap_layers = theme(axis.line=element_blank()),
          file = "dataDriven_withCorrelations/covarianceMatrixHeatmap_dd_new.pdf")

## PCA -----------------------------------------

pca = prcomp(U.pca$tPCA)
summary(pca)

group <- colnames(U.pca$tPCA)

pca_values <- pca$x
df <- data.frame(
  pc1 <- as.numeric(pca_values[,1]),
  pc2 <- as.numeric(pca_values[,2]),
  g <- group
)
row.names(df) <- row.names(pca_values)

tiff(paste0(folder,"PCA_mash_dataDriven.tiff"), units="in", width=4, height=3, res=350)
ggplot(df,aes(x=pc1,y=pc2)) + 
  theme_classic()+geom_point(aes(color=g), size = 2, stroke = 1)+ 
  labs(x="PC1 (75.65%)", y="PC2 (20.21%)") # +
dev.off()

