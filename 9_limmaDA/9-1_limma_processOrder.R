library("ggplot2")
library ("DESeq2")
library("pheatmap")
library("wesanderson")
library(statmod)
library(goseq)
library("RColorBrewer")
library(edgeR)
library(ggfortify)
library(devtools)
library("FactoMineR")
library("factoextra")
library(tidyverse)
library(cluster)
library(factoextra)
library(ggrepel)
library(data.table)
library(VennDiagram)


setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks/")
OUT_DIR <- "../10_differentialAccessibility/"

# import raw counts -------------------------------------------------------
rawdata <- read.table("countMatrix1.txt", check.names=FALSE, stringsAsFactors=FALSE, header=TRUE) 
# outliers <- c("Fluni3", "HC3")
pbs <- rawdata[c(22:27)]
rawdata_reorder <- rawdata[-c(10,13,22,23,24,25,26,27)]
# removes Fluni3, HC3
rawdata <- cbind(rawdata_reorder,pbs)
# adds back all 6 pbs



# Create DGEList object ---------------------------------------------------
y <- DGEList(counts=rawdata[,2:length(colnames(rawdata))], genes=rawdata[,1]) 
y <- calcNormFactors(y) #converts observed library sizes into effective library sizes


# Meta Data ---------------------------------------------------------------

meta_data <- read.table(file="../../../analysis/10_differentialAccessibility/meta_data_processOrder.txt", header=TRUE)
row.names(meta_data) <- meta_data$sample_name
meta_data <- subset(meta_data, 
                    meta_data$sample_name %in% colnames(rawdata[2:length(colnames(rawdata))]))
meta_data <- meta_data[,2:length(colnames(meta_data))]
meta_data$molecule <- as.factor(meta_data$molecule)
meta_data$molecule <- relevel(meta_data$molecule, ref = "pbs")
meta_data$processOrder <- as.numeric(meta_data$processOrder)

# Confirm meta data and count matrix sample index 

all(rownames(meta_data) == colnames(rawdata[2:length(rawdata)]))


# Create a design matrix --------------------------------------------------
design <- model.matrix(~1 + molecule + processOrder, data=meta_data)
colnames(design) <- c("intercept","5f","bg", "fen", "fluni",
                     "hc", "hq", "myr", "nerol","processOrder")


# Filter out low-expressed genes ------------------------------------------

v <- voom(y, design, plot=TRUE) #voom to get the cpm values
cpm <- v$E
A <- rowMedians(cpm) 

y2 <- DGEList(counts=rawdata[,2:length(colnames(rawdata))], genes=rawdata[,1]) 
#filter out low-expressed genes from the original count  matrix
y_filtered <- y2[A>1,]
y_filtered <- calcNormFactors(y_filtered)


# Voom on filtered count matrix -------------------------------------------

v2 <- voom(y_filtered, design, plot=TRUE)
cpm2 <- v2$E



# Fit linear model no contrasts -------------------------------------------

fit <- lmFit(v2, design)
fit <- eBayes(fit)
summary(decideTests(fit))
saveRDS(fit, file = paste0(OUT_DIR, "fit_allPBS_processOrder.rds"))

expression <- v2$E  
write.table(expression, file=paste0("../10_differentialAccessibility/expressionCPMs_allPBS_processOrder.txt"), quote = FALSE, sep = ",") 

# PCA ---------------------------------------------------------------------

pca = prcomp(t(expression))
summary(pca)

group <- meta_data$molecule

pca_values <- pca$x
df <- data.frame(
  pc1 <- as.numeric(pca_values[,1]),
  pc2 <- as.numeric(pca_values[,2]),
  molecule <- meta_data$molecule,
  g <- group
)
row.names(df) <- row.names(pca_values)

# myColorScale <- c('wt_pup'='#B7E4F9', 'ko_pup'='#24325F', 'wt_ad'='#FAE48B', 'ko_ad' = '#E89242')

tiff(paste0(OUT_DIR,"PCA_allPBS_processOrder.tiff"), units="in", width=4, height=3, res=350)
ggplot(df,aes(x=pc1,y=pc2)) + 
  theme_classic()+geom_point(aes(color=g), size = 2, stroke = 1)+ 
  labs(x="PC1 (17.89%)", y="PC2 (14.26%)") # +
  # geom_label(
  #  label=rownames(meta_data), 
  #  nudge_x = 0.25, nudge_y = 0.25)
dev.off()

######################
# End of script
######################
# Test for outliers driving intercept variation --------------------------
# Raul advised on this approach to DA peak validation
# Find 10 most differentially accessible peaks in each category, including intercept
results_intercept <- topTable(fit, coef = "intercept", n = Inf, sort = "p")
results_5f <- topTable(fit, coef = "5f", n = Inf, sort = "p")
results_bg <- topTable(fit, coef = "bg", n = Inf, sort = "p")
results_fen <- topTable(fit, coef = "fen", n = Inf, sort = "p")
results_fluni <- topTable(fit, coef = "fluni", n = Inf, sort = "p")
results_hc <- topTable(fit, coef = "hc", n = Inf, sort = "p")
results_hq <- topTable(fit, coef = "hq", n = Inf, sort = "p")
results_myr <- topTable(fit, coef = "myr", n = Inf, sort = "p")
results_nerol <- topTable(fit, coef = "nerol", n = Inf, sort = "p")

top10_bg <- head(results_bg[order(results_bg$adj.P.Val),],10)
top10_5f <- head(results_bg[order(results_bg$adj.P.Val),],10)
top10_intercept <- head(results_intercept[order(results_intercept$adj.P.Val),],10)

results_bg$DA <- "NO"
results_bg$DA[results_bg$logFC > 1 & results_bg$adj.P.Val < 0.01] <- "UP"
results_bg$DA[results_bg$logFC < -1 & results_bg$adj.P.Val < 0.01] <- "DOWN"

results_5f$DA <- "NO"
results_5f$DA[results_5f$logFC > 1 & results_5f$adj.P.Val < 0.01] <- "UP"
results_5f$DA[results_5f$logFC < -1 & results_5f$adj.P.Val < 0.01] <- "DOWN"

results_fen$DA <- "NO"
results_fen$DA[results_fen$logFC > 1 & results_fen$adj.P.Val < 0.01] <- "UP"
results_fen$DA[results_fen$logFC < -1 & results_fen$adj.P.Val < 0.01] <- "DOWN"

results_fluni$DA <- "NO"
results_fluni$DA[results_fluni$logFC > 1 & results_fluni$adj.P.Val < 0.01] <- "UP"
results_fluni$DA[results_fluni$logFC < -1 & results_fluni$adj.P.Val < 0.01] <- "DOWN"

results_hc$DA <- "NO"
results_hc$DA[results_hc$logFC > 1 & results_hc$adj.P.Val < 0.01] <- "UP"
results_hc$DA[results_hc$logFC < -1 & results_hc$adj.P.Val < 0.01] <- "DOWN"

results_hq$DA <- "NO"
results_hq$DA[results_hq$logFC > 1 & results_hq$adj.P.Val < 0.01] <- "UP"
results_hq$DA[results_hq$logFC < -1 & results_hq$adj.P.Val < 0.01] <- "DOWN"

results_myr$DA <- "NO"
results_myr$DA[results_myr$logFC > 1 & results_myr$adj.P.Val < 0.01] <- "UP"
results_myr$DA[results_myr$logFC < -1 & results_myr$adj.P.Val < 0.01] <- "DOWN"

results_nerol$DA <- "NO"
results_nerol$DA[results_nerol$logFC > 1 & results_nerol$adj.P.Val < 0.01] <- "UP"
results_nerol$DA[results_nerol$logFC < -1 & results_nerol$adj.P.Val < 0.01] <- "DOWN"

head(results_bg[order(results_bg$adj.P.Val) & results_bg$DA == 'NO', ])

# Plot the expression of that peak for each replicate
df <- data.frame()
for (i in 1:length(rownames(top10_bg))){
  peakNumber <- rownames(top10_bg[i,])
  peakE<- cpm2[peakNumber,]
  df<-rbind(df,peakE)
}
colnames(df) <- colnames(cpm2)
rownames(df) <- rownames(top10_bg)

reshape2::melt(as.matrix(df)) %>% 
  mutate(condition = gsub(pattern = "-.*", "",Var2)) %>% 
  ggplot(aes(y=value, x=condition)) + 
    geom_point(size=2) + 
    geom_text(aes(label=Var2)) + 
    facet_wrap(~Var1, scales = 'free', nrow = 3) + 
    theme_minimal()

df_intercept <- data.frame()
for (i in 1:length(rownames(top10_intercept))){
  peakNumber <- rownames(top10_intercept[i,])
  peakE<- cpm2[peakNumber,]
  df_intercept<-rbind(df_intercept,peakE)
}
colnames(df_intercept) <- colnames(cpm2)
rownames(df_intercept) <- rownames(top10_intercept)

reshape2::melt(as.matrix(df_intercept)) %>% 
  mutate(condition = gsub(pattern = "-.*", "",Var2)) %>% 
  ggplot(aes(y=value, x=condition)) + 
  geom_point(size=2) + 
  geom_text(aes(label=Var2)) + 
  facet_wrap(~Var1, scales = 'free', nrow = 3) + 
  theme_minimal()

# Volcano Plots -------------------------------------
png(file= "../10_differentialAccessibility/results_bg.png", width = 600, height = 400)
ggplot(data = results_bg, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Beta Glucan') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "../10_differentialAccessibility/results_5f.png", width = 600, height = 400)
ggplot(data = results_5f, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: 5F') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_fen.png", width = 600, height = 400)
ggplot(data = results_fen, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Fenoterol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) 
dev.off()

png(file= "../10_differentialAccessibility/results_fluni.png", width = 600, height = 400)
ggplot(data = results_fluni, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Flunisolide') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_hc.png", width = 600, height = 400)
ggplot(data = results_hc, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydrocortisone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_hq.png", width = 600, height = 400)
ggplot(data = results_hq, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Hydroquinone') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_myr.png", width = 600, height = 400)
ggplot(data = results_myr, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Myricetin') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

png(file= "../10_differentialAccessibility/results_nerol.png", width = 600, height = 400)
ggplot(data = results_nerol, aes(x = logFC, y = -log10(adj.P.Val), col = DA)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.01), col = "gray", linetype = 'dashed') + 
  geom_point(size=0.1) +
  ggtitle('Differentially Accessible Peaks: Nerol') +
  scale_color_manual(values = c("#009999", "grey", "red"), 
                     # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated"))
dev.off()

# with help from this plot tutorial https://biostatsquid.com/volcano-plots-r-tutorial/

# Barplot of DA Peaks: Up vs Down -----------------------------------------

bg <- c(sum(results_bg$DA == "UP"),sum(results_bg$DA == "DOWN"))
fiveF <- c(sum(results_5f$DA == "UP"),sum(results_5f$DA == "DOWN"))
fen <- c(sum(results_fen$DA == "UP"),sum(results_fen$DA == "DOWN"))
fluni <- c(sum(results_fluni$DA == "UP"),sum(results_fluni$DA == "DOWN"))
hc <- c(sum(results_hc$DA == "UP"),sum(results_hc$DA == "DOWN"))
hq <- c(sum(results_hq$DA == "UP"),sum(results_hq$DA == "DOWN"))
myr <- c(sum(results_myr$DA == "UP"),sum(results_myr$DA == "DOWN"))
nerol <- c(sum(results_nerol$DA == "UP"),sum(results_nerol$DA == "DOWN"))

countsDA <- data.frame(myr, fen, hc, fluni,bg, hq,fiveF, nerol)
rownames(countsDA) <- c("Upregulated", "Downregulated")
countsDA <- reshape2::melt(as.matrix(countsDA)) 
colnames(countsDA) <- c("Direction","Treatment","PeakCount")

png(file= "../10_differentialAccessibility/totalPeakCount_p_0.01.png", width = 500, height = 400)
ggplot(countsDA, aes(fill=Direction, y=PeakCount, x=Treatment)) + 
  geom_bar(position="stack", stat="identity") + 
  ggtitle('Number of Differentially Accessible Peaks Per Condition (p<0.01) ')
dev.off()
