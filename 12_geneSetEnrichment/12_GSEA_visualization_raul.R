library(ggplot2)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/12_geneSetEnrichment/")
test <- readRDS("fgseaResult.rds")

gsea_res <- bind_rows(test, .id="stims")
gsea_res %>%
  #filter(pathways == “up_il17") %>%
  arrange(padj) %>%
  filter(padj < 0.2)
top_gos <- gsea_res %>%
  group_by(stims) %>%
  arrange(padj)
#top_gos <- top_gos %>% filter(padj < 0.2)
#top_gos <- top_gos %>%
#  filter(grepl(grep_pattern, Description))
top_gos <- top_gos %>%
  #slice_head(n=10) %>%
  ungroup() %>%
  pull(pathways) %>%
  unique()
NES_mat <- gsea_res %>%
  filter(pathways %in% top_gos) %>%
  select(stims,pathways,NES) %>%
  spread(pathways,NES)
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
go_bp_dat <- gsea_res %>%
  filter(pathways %in% top_gos) %>%
  mutate(label = str_to_title(gsub("HALLMARK_","",pathways))) %>%
  mutate(label = factor(label, levels = c(top_gos_final_label))) %>%
  ungroup()
go_bp_plot <- go_bp_dat %>%
  group_by(pathways) %>%
  # filter(any(padj < 0.05)) %>%
  ungroup() %>%
  drop_na() %>%
  ggplot(aes(y=label, x=stims))+
  geom_point(aes(size=-log10(padj), fill=NES, color=padj<0.1), shape=21, stroke=1) +
  # scale_fill_continuous_diverging(rev=F, palette= "Red-Green") +
  scale_color_manual(values = c("lightgray","#53565A")) +
                                 scale_size(range = c(0.75,4.5)) +
                                   ylab("") +
                                   xlab("") +
                                   guides(fill = guide_colourbar(barwidth = 0.5, barheight = 7)) +
                                   # ipub_th +
                                   theme(strip.text = element_text(size=9, family = "Helvetica",angle=0, hjust = 0.5),
                                         axis.text.x = element_text(size=7,angle=45, hjust = 1),
                                         legend.position = "right",
                                         strip.background = element_blank(),
                                         legend.box="vertical")
                                 gsea_deg_bar <- go_bp_plot # + plot_layout(heights = c(0.2,1)) #+ 
                                 # DEG_res_bar
                                 pre_tt_gsea_dot_plot_file <- "response_gsea_dot_plot.pdf"
                                 ggsave(plot = gsea_deg_bar,
                                        filename = pre_tt_gsea_dot_plot_file,
                                        height = 7, width = 7) (edited) 
                                 