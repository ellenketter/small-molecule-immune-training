library(readr)
library(harmony)
library(stringr)
library(pbapply)
library(parallel)
library(tidyr)
library(dplyr)
library(colorspace)
library("ggplot2")

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/15_transcriptionFactorMotifs/")
OUT_DIR <- "15_homer_motif_plots/"
# dir.create(OUT_DIR)


clusters <- c("bg", "5f","fen","fluni", "hc", "hq","myr", "nerol")

results_list <- list()
sig_motifs <- vector()
for(i in 1:length(clusters)){
  name <- clusters[i]
  
  results <- read.delim(file=paste0("tfMotifOutput_UP_", name, "/knownResults.txt"), header=TRUE)
  sig_motifs <- append(sig_motifs, subset(results, results$q.value..Benjamini. < 0.05)$Motif.Name)

  assign(paste0(name, "_results"), results)
  
  results_list[[i]] <- results
  names(results_list)[i] <- paste0(name, "_results")

  print(i)
  
}
sig_motifs <- unique(sig_motifs)



motif <- vector()
pval <- vector()
qval <- vector()
enrichment <- vector()
cluster_name <- vector()

for(i in 1:length(results_list)){
  cluster <- clusters[i]
  result <- results_list[[i]]
  subset <- which(result$Motif.Name %in% sig_motifs)
  motif <- append(motif, result$Motif.Name[subset])
  pval <- append(pval, result$P.value[subset])
  qval <- append(qval, result$q.value..Benjamini.[subset])
  target_perc <- as.numeric(sub("%", "", result$X..of.Target.Sequences.with.Motif[subset]))
  background_perc <- as.numeric(sub("%", "", result$X..of.Background.Sequences.with.Motif[subset])) 
  enrichment <- append(enrichment, target_perc/background_perc)
  cluster_name <- append(cluster_name, rep(cluster, length(subset)))
}

group <- vector(length=length(qval))
group[which(qval <= 0.05)] <- "sig"
group[which(qval > 0.05)] <- "not_sig"

motif <- gsub("\\(.*","", motif)


df <- data.frame(
  clust <- cluster_name,
  m <- motif,
  e <- enrichment,
  q <- qval ,
  g <- group)
colnames(df) <- c("clust", "m", "e", "q", "g")
uMotif<- unique(motif)
df$m <- factor(df$m, c(uMotif))

o <- split(df, df$clust)
oo <- data.frame(o$bg$q, o$`5f`$q, o$fen$q, o$fluni$q, o$hc$q, o$hq$q, o$myr$q, o$nerol$q)
oo$m<- o$bg$m
colnames(oo) <- c(clusters,"m")
library(tidyverse)
ooo<- oo[1:105,]

ooo<- arrange(ooo, bg)
ooo <- ooo[1:50,]

df1 <- reshape2::melt(ooo)
names(df1) <- c("m", "clust", "q")


# tiff(paste0(OUT_DIR,"enriched_motifs_q0.1_peaks0.01.tiff"), units="in", width=9, height=6, res=400)
# ggplot(df, aes(x = clust, y = m)) + 
#   geom_point(aes(size = enrichment, fill = g), alpha = 0.75, shape = 22) + theme_light() + scale_fill_manual(values=c("grey", "blue"))+
#   scale_size_continuous(limits = c(1.1, 3)) + 
#   labs( x= "", y = "", size = "Fold enrichment", fill = "Significance q.Benjamini < 0.1") + theme(axis.text.y = element_text(size=6))+
#   theme(axis.text.x = element_text(angle = 45, vjust = 0.8, hjust=1))
# dev.off()


pdf(paste0(OUT_DIR,"enriched_motifs_q0.1_final.pdf"), width=5, height=10)

df1 %>% 
  mutate(clust=factor(clust, levels=c("fluni", "hc", "myr","fen", "nerol", "5f", "hq", "bg"))) %>%
  ggplot(aes(x = clust, y = m)) + 
  geom_point(aes(fill = q, color=q<0.05, size = -q), alpha = 0.75, shape = 21) +
  scale_color_manual(values = c("lightgray", "#53565A")) +
  scale_size("q",breaks=c(-0.2,-0.15,-0.10,-0.05,0),labels=c(0.2,0.15,0.10,0.05,0)) +
  theme_classic() + 
  scale_fill_gradient(low="brown", high="white", 
                      oob=scales::squish)+
  # scale_size_continuous(limits = c(1.1, 3)) + 
  labs( x= "", y = "", size = "q", fill = "q") + theme(axis.text.y = element_text(size=9))+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size=7))
dev.off()

