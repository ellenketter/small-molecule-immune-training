library(ggplot2)
library(dplyr)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/15_transcriptionFactorMotifs/")
new <- readRDS("dotPlotData.rds")

new %>%
  #filter(pathways == “up_il17") %>%
  arrange(pAdj) %>%
  filter(pAdj > 0.05)
top_gos <- new %>%
  group_by(molecule) %>%
  arrange(pAdj)
#top_gos <- top_gos %>% filter(padj < 0.2)
#top_gos <- top_gos %>%
#  filter(grepl(grep_pattern, Description))
top_gos <- top_gos %>%
  slice_head(n=10) %>%
  ungroup() %>%
  pull(motif) %>%
  unique()
NES_mat <- new %>%
  filter(motif %in% top_gos) %>%
  select(molecule,motif,enrichment) %>%
  spread(motif,enrichment)
NES_mat <- as.data.frame(NES_mat)
rownames(NES_mat) <- NES_mat[,1]
NES_mat <- NES_mat[,-1]
NES_mat <- as.matrix(NES_mat)
NES_mat <- t(NES_mat)
#NES_mat <- NES_mat[which(apply(NES_mat, 1, function(x)sum(is.na(x))) == 0),]
NES_hclust <- hclust(d = dist(NES_mat))
top_gos_final <- rownames(NES_mat)[NES_hclust$order]
top_gos_final_label <- gsub("HALLMARK_","",top_gos_final)
top_gos_final_label <- str_to_title(top_gos_final_label)
go_bp_dat <- new %>%
  filter(motif %in% top_gos) %>%
  mutate(label = str_to_title(gsub("HALLMARK_","",motif))) %>%
  mutate(label = factor(label, levels = c(top_gos_final_label))) %>%
  ungroup()
go_bp_plot <- go_bp_dat %>%
  group_by(motif) %>%
  # filter(any(padj < 0.05)) %>%
  ungroup() %>%
  drop_na() %>%
  ggplot(aes(y=label, x=molecule))+
  geom_point(aes(size=-log10(pAdj), fill=log2(enrichment), color=pAdj<0.005), shape=21, stroke=1) +
  scale_fill_continuous_diverging(rev=F, palette= "Red-Green") +
  scale_color_manual(values = c("lightgray","#53565A")) +
  scale_size(range = c(0.75,4.5)) +
  ylab("") +
  xlab("") +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 7)) +
  # ipub_th +
  theme(panel.background = element_rect(fill = "white", color = "black"),strip.text = element_text(size=9, family = "Helvetica",angle=0, hjust = 0.5),
        axis.text.x = element_text(size=7,angle=45, hjust = 1),
        legend.position = "right",
        strip.background = element_blank(),
        legend.box="vertical")
gsea_deg_bar <- go_bp_plot # + plot_layout(heights = c(0.2,1)) #+ 
# DEG_res_bar
pre_tt_gsea_dot_plot_file <- "TF_dot_plot_pAdj-top10.pdf"
ggsave(plot = gsea_deg_bar,
       filename = pre_tt_gsea_dot_plot_file,
       height = 6, width = 5) # (edited) 
