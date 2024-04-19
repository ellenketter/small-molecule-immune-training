library(ggplot2)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/fgseaTwoList/")

results <- readRDS("humanPhenotypeOntology.rds")

sig_pathways <- vector()

bg_results <- head(subset(results$bg, results$bg$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(bg_results$gs_name))

`5f_results` <- head(subset(results$`5f`, results$`5f`$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(`5f_results`$gs_name))

fen_results <- head(subset(results$fen, results$fen$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(fen_results$gs_name))

fluni_results <- head(subset(results$fluni, results$fluni$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(fluni_results$gs_name))

hc_results <- head(subset(results$hc, results$hc$overlap >=15),15)
sig_pathways <- append(sig_pathways, as.character(hc_results$gs_name))

hq_results <- head(subset(results$hq, results$hq$overlap >=15),15)
sig_pathways <- append(sig_pathways, as.character(hq_results$gs_name))

myr_results <- head(subset(results$myr, results$myr$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(myr_results$gs_name))

nerol_results <- head(subset(results$nerol, results$nerol$overlap >=15), 15)
sig_pathways <- append(sig_pathways, as.character(nerol_results$gs_name))

sig_pathways <- unique(sig_pathways)
results_list <- list(`5f_results`, bg_results, fen_results, fluni_results, hc_results,
                     hq_results, myr_results, nerol_results) 

pathway <- vector()
cluster_name <- vector()
padj <- vector()
molecules_renamed <- names(results)


for(i in 1:length(results_list)){
  cluster <- molecules_renamed[i]
  result <- results_list[[i]]
  subset <- which(result$gs_name %in% sig_pathways)
  pathway <- append(pathway, as.character(result$gs_name[subset]))
  padj <- append(padj, result$padj[subset])
  cluster_name <- append(cluster_name, rep(cluster, length(subset)))
}

# Plot padj < 0.05 -----------------------------------------------------------
group <- vector(length=length(padj))
group[which(padj <= 0.05)] <- "padj <= 0.05"
group[which(padj > 0.05)] <- "padj > 0.05"

df <- data.frame(
  clust <- cluster_name,
  path <- pathway,
  q <- padj,
  g <- group
)
colnames(df) <- c("clust", "path", "padj", "g")
df$g <- factor(df$g, c("padj > 0.05", "padj <= 0.05"))
df$clust <- factor(df$clust, molecules_renamed)
# sig_pathways <- unique(sig_pathways)
df$path <- factor(df$path, c(sig_pathways))

tiff(paste0("enriched_reactome_pathways_padj0.05.tiff"), units="in", width=8, height=8, res=400)
ggplot(df, aes(x = clust, y = path)) +
  geom_point(aes(size = -log10(padj), fill = g), alpha = 0.75, shape = 21) + theme_classic() + 
  scale_fill_manual(values=c("white", "#D82632")) +
  scale_size_continuous(limits = c(0, 10)) +
  labs( x= "", y = "", size = "Significance: -log10 padj value", fill = "Significance padj <= 0.05") + 
  theme(axis.text.y = element_text(size=6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

