library(ggplot2)
library(tidyr)
library(stringr)
library(dplyr)
library(colorspace)
library(heatmaply)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/fgseaTwoList/treatmentOverlaps/vennDiagram//")
test <- readRDS("fgsea_ora_reactome_top5perTreatment_venn_0.1.rds")
# test1 <- apply(test$noHcFluni,2,as.character)
# write.csv(test1, "noHcFluni.csv")


gsea_res <- bind_rows(test, .id="stims")

# gsea_res <- subset(gsea_res, gsea_res$padj < 0.1) 
gsea_res <- arrange(gsea_res, gsea_res$padj) 
  
top_gos <- gsea_res %>%
  group_by(stims) %>%
  arrange(padj)
#top_gos <- top_gos %>% filter(padj < 0.2)
#top_gos <- top_gos %>%
#  filter(grepl(grep_pattern, Description))
top_gos <- top_gos %>%
  #slice_head(n=10) %>%
  ungroup() %>%
  pull(gs_description) %>%
  unique()
NES_mat <- gsea_res %>%
  filter(gs_description %in% top_gos) %>%
  select(stims,gs_description,percentOverlap) %>%
  spread(gs_description,percentOverlap)
NES_mat <- as.data.frame(NES_mat)
rownames(NES_mat) <- NES_mat[,1]
NES_mat <- NES_mat[,-1]
NES_mat <- as.matrix(NES_mat)
NES_mat <- t(NES_mat)
#NES_mat <- NES_mat[which(apply(NES_mat, 1, function(x)sum(is.na(x))) == 0),]

NES_mat[is.na(NES_mat)] <- 0
NES_hclust <- hclust(d = dist(NES_mat))
top_gos_final <- rownames(NES_mat)[NES_hclust$order]
top_gos_final_label <- gsub("HALLMARK_","",top_gos_final)
top_gos_final_label <- str_to_title(top_gos_final_label)
go_bp_dat <- gsea_res %>%
  filter(gs_description %in% top_gos) %>%
  mutate(label = str_to_title(gsub("HALLMARK_","",gs_description))) %>%
  mutate(label = factor(label, levels = c(top_gos_final_label))) %>%
  ungroup()
go_bp_plot <- go_bp_dat %>%
  group_by(gs_description) %>%
  # filter(any(padj < 0.05)) %>%
  ungroup() %>%
  drop_na() %>%
  ggplot(aes(y=label, x=stims))+
  geom_point(aes(size=-log10(padj), fill=percentOverlap, color=padj<0.1), shape=21, stroke=1) +
  scale_fill_continuous_diverging(rev=F, palette= "Red-Green") +
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
                                 gsea_deg_bar <-go_bp_plot # + plot_layout(heights = c(0.2,1))
                                 pre_tt_gsea_dot_plot_file <- "reactome_fgsea_ora_dot_plot_top5perTreatment_padj0.1.pdf"
                                 ggsave(plot = gsea_deg_bar,
                                        filename = pre_tt_gsea_dot_plot_file,
                                        height = 6, width = 8)

                                 
heatmaply(cor(t(NES_mat)), file = "./reactome_gsea_overlap_pathwayHeatmap_dd_0.1.html")
