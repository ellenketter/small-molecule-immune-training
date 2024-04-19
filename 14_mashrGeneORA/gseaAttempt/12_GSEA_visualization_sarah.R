library(ggplot2)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/homerAnnotation/")
OUT_DIR <- "../homerOverenrichment/"
# dir.create(OUT_DIR)


# Reactome pathways -------------------------------------------------------

molecules <- c("bg","fen", "myr","hc","fluni","5f","hq","nerol")
sig_pathways <- vector()
for(i in 1:length(molecules)){
  name <- molecules[i]
  
  results <- read.delim(file=paste0(name, "_up_GO/reactome.txt"), 
                        header=TRUE)
  results$pval <- 10^(results$logP)
  results$q <- p.adjust(results$pval, method = "BH")
  results_sub1 <- subset(results, results$q < 0.05)[1:20,]
  # instead of just subsetting the 1st 15 by q, select by enrichment
  results_sub2 <- subset(results_sub1, results_sub1$Target.Genes.in.Term >= 15)
  sig_pathways <- append(sig_pathways, as.character(results_sub2$Term))
  
  assign(paste0(name, "_results"), results)
  
}
sig_pathways <- unique(sig_pathways)
results_list <- list(bg_results, fen_results, hc_results,fluni_results,
                     myr_results,`5f_results`,hq_results, nerol_results) 

pathway <- vector()
pval <- vector()
enrichment <- vector()
cluster_name <- vector()
qval <- vector()

molecules_renamed <- molecules

for(i in 1:length(results_list)){
  cluster <- molecules_renamed[i]
  result <- results_list[[i]]
  subset <- which(result$Term %in% sig_pathways)
  pathway <- append(pathway, as.character(result$Term[subset]))
  pval <- append(pval, result$pval[subset])
  qval <- append(qval, result$q[subset])
  cluster_name <- append(cluster_name, rep(cluster, length(subset)))
}


# plot p<0.05  ------------------------------------------------------------


group <- vector(length=length(pval))
group[which(pval <= 0.05)] <- "sig"
group[which(pval > 0.05)] <- "not_sig"

df <- data.frame(
  clust <- cluster_name,
  path <- pathway,
  p <- pval,
  g <- group
)
colnames(df) <- c("clust", "path", "p", "g")
df$clust <- factor(df$clust, molecules)
sig_pathways <- unique(sig_pathways)
df$path <- factor(df$path, c(sig_pathways))

tiff(paste0(OUT_DIR,"enriched_reactome_pathways_p0.05.tiff"), units="in", width=10, height=6, res=400)
ggplot(df, aes(x = clust, y = path)) +
  geom_point(aes(size = -log10(pval), fill = g), alpha = 0.75, shape = 21) + theme_light() + scale_fill_manual(values=c("grey", "red"))+
  scale_size_continuous(limits = c(0, 10)) +
  labs( x= "", y = "Pathway", size = "Significance", fill = "Significance p < 0.05") + theme(axis.text.y = element_text(size=6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1))
dev.off()


# Plot q < 0.05 -----------------------------------------------------------
group <- vector(length=length(qval))
group[which(qval <= 0.05)] <- "q <= 0.05"
group[which(qval > 0.05)] <- "q > 0.05"

df <- data.frame(
  clust <- cluster_name,
  path <- pathway,
  q <- qval,
  g <- group
)
colnames(df) <- c("clust", "path", "q", "g")
df$g <- factor(df$g, c("q > 0.05", "q <= 0.05"))
df$clust <- factor(df$clust, molecules)
# sig_pathways <- unique(sig_pathways)
df$path <- factor(df$path, c(sig_pathways))

tiff(paste0(OUT_DIR,"enriched_reactome_pathways_q0.05.tiff"), units="in", width=8, height=8, res=400)
ggplot(df, aes(x = clust, y = path)) +
  geom_point(aes(size = -log10(pval), fill = g), alpha = 0.75, shape = 21) + theme_classic() + scale_fill_manual(values=c("white", "#D82632"))+
  scale_size_continuous(limits = c(0, 10)) +
  labs( x= "", y = "", size = "Significance: -log10 q value", fill = "Significance q <= 0.05") + theme(axis.text.y = element_text(size=6))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
dev.off()

