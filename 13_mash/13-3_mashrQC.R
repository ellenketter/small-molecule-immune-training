library(dplyr)
library(ggplot2)
library(grid)
library(gridExtra)
library(reshape2)

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "correlatingBetas/"
# dir.create(OUT_DIR)

mashInput <- readRDS("../../10_differentialAccessibility/fit_allPBS_processOrder.rds")
mashInputSD <- data.frame(mashInput$stdev.unscaled, mashInput$genes)
mashInputSD <- select(mashInputSD, -c("intercept", "processOrder"))
mashInputSD<- melt(mashInputSD)
names(mashInputSD) <- c("peak", "treatment", "PriorSD")

mashInputBetas <- data.frame(mashInput$coefficients, mashInput$genes)
mashInputBetas <- select(mashInputBetas, -c("intercept", "processOrder"))
mashInputBetas<- melt(mashInputBetas)
names(mashInputBetas) <- c("peak", "treatment", "PriorBeta")


mashOutput <- readRDS("mashResults_dd.rds") 
mashOutputSD <- mashOutput$result$PosteriorSD
mashOutputSD <- melt(mashOutputSD)
names(mashOutputSD) <- c("peak", "treatment", "PosteriorSD")

mashOutputBetas <- mashOutput$result$PosteriorMean
mashOutputBetas <- melt(mashOutputBetas)
names(mashOutputBetas) <- c("peak", "treatment", "PosteriorBeta")

mashInput <- inner_join(mashInputBetas,mashInputSD)
levels(mashInput$treatment)[levels(mashInput$treatment)=="X5f"] <- "5f"
mashOutput <- inner_join(mashOutputBetas,mashOutputSD)
nrow(mashInput) == nrow(mashOutput)

mashPriorsAndPosteriors <- inner_join(mashInput,mashOutput)
mashPriorsAndPosteriors <- split(mashPriorsAndPosteriors, f = mashPriorsAndPosteriors$treatment)

# Plot correlation between priors and posteriors
for (i in 1:length(names(mashPriorsAndPosteriors))){
  
  i_file_name <- paste0("correlatingBetas/priorPosteriorBetaCorrelation_",names(mashPriorsAndPosteriors)[i],".tiff")
  
  p<- mashPriorsAndPosteriors[[names(mashPriorsAndPosteriors)[i]]] %>%  
    ggplot(aes(x=PriorBeta, 
               y=PosteriorBeta)) + 
    geom_point(shape=21, color="black", size=0.1) + 
    ggtitle(label = paste0("Correlation of Betas: ", names(mashPriorsAndPosteriors)[i])) +
    theme_classic() + 
    labs(x="Prior Beta (LogFC)", y="Posterior Beta", fill="")+
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0) + 
    geom_smooth(method = "lm", se=FALSE) +
    theme(axis.text=element_text(size=15.5),
          axis.title=element_text(size=15.5)) + 
    theme(legend.text = element_text(size=15.5))
  
  ggsave(filename = i_file_name, plot = p, width=4.5, height=3)
  
}

# Per peak, show mean and SD of prior vs posterior beta
# Picked peaks which are DA in all conditions, according to MASH
allIntersect <- readRDS("../canonicalCovariance/allTreatmentIntersectPeaks.rds")
mashPriorsAndPosteriors <- inner_join(mashInput,mashOutput)

ppPlot <- list()
for(i in 1:100){
prior <- ggplot(mashPriorsAndPosteriors[mashPriorsAndPosteriors$peak == allIntersect[i],], 
                aes(x= PriorBeta, y= treatment, color = treatment)) +
  geom_dotplot(binaxis='y', binwidth = 0.2, stackdir='center') +
  geom_errorbar(aes(xmax = PriorBeta + PriorSD, xmin = PriorBeta - PriorSD),
                position = "dodge") +
  ggtitle(label = allIntersect[i]) +
  theme(legend.position = "none") + 
  xlim(-1,1) 
posterior <- ggplot(mashPriorsAndPosteriors[mashPriorsAndPosteriors$peak == allIntersect[i],], 
                    aes(x= PosteriorBeta, y= treatment, color = treatment)) +
  geom_dotplot(binaxis='y',binwidth = 0.2, stackdir='center') +
  geom_errorbar(aes(xmax = PosteriorBeta + PosteriorSD, xmin = PosteriorBeta - PosteriorSD),
                position = "dodge") + 
  theme(legend.position = "none") + 
  xlim(-1,1)
ppPlot[[i]] <- grid.arrange(prior, posterior)
print(i)
}
dev.off()

pdf(paste0("ppPlot.pdf"))
marrangeGrob(ppPlot, nrow = 2, ncol = 2)
dev.off()

mashInput$coefficients[mashInput$genes == "chr17.84956870.84958199"]
topTable(mashInput)[2,]


# Compare DD + common baseline posterior betas vs DD + corr posteriors ------

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/dataDriven_withCorrelations/")
OUT_DIR <- "correlatingBetas/"
# dir.create(OUT_DIR)

mashOutput <- readRDS("mashResults_dd.rds") 
mashOutputSD <- mashOutput$result$PosteriorSD
mashOutputSD <- melt(mashOutputSD)
names(mashOutputSD) <- c("peak", "treatment", "PosteriorSD")

mashOutputBetas <- mashOutput$result$PosteriorMean
mashOutputBetas <- melt(mashOutputBetas)
names(mashOutputBetas) <- c("peak", "treatment", "PosteriorBeta")

setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/otherModels/dataDriven_commonBaseline/")

mashOutput_cb <- readRDS("mashResults_dd.rds")
colnames(mashOutput_cb$result$PosteriorMean) <- gsub("-pbs","",colnames(mashOutput_cb$result$PosteriorMean))
colnames(mashOutput_cb$result$PosteriorSD) <- gsub("-pbs","",colnames(mashOutput_cb$result$PosteriorSD))

mashOutput_cbSD <- mashOutput_cb$result$PosteriorSD
mashOutput_cbSD <- melt(mashOutput_cbSD)
names(mashOutput_cbSD) <- c("peak", "treatment", "cbPosteriorSD")

mashOutput_cbBetas <- mashOutput_cb$result$PosteriorMean
mashOutput_cbBetas <- melt(mashOutput_cbBetas)
names(mashOutput_cbBetas) <- c("peak", "treatment", "cbPosteriorBeta")


mashOutput_cb <- inner_join(mashOutput_cbBetas,mashOutput_cbSD)

mashOutput <- inner_join(mashOutputBetas,mashOutputSD)
levels(mashOutput$treatment)[levels(mashOutput$treatment)=="X5f"] <- "5f"
nrow(mashOutput_cb) == nrow(mashOutput)

mashPriorsAndPosteriors <- inner_join(mashOutput_cb,mashOutput)
mashPriorsAndPosteriors <- split(mashPriorsAndPosteriors, f = mashPriorsAndPosteriors$treatment)

# Plot correlation between priors and posteriors
for (i in 1:length(names(mashPriorsAndPosteriors))){
  
  i_file_name <- paste0("correlatingBetas/corrVScbPosteriorBetaCorrelation_",names(mashPriorsAndPosteriors)[i],".tiff")
  
  p<- mashPriorsAndPosteriors[[names(mashPriorsAndPosteriors)[i]]] %>%  
    ggplot(aes(x=PosteriorBeta, 
               y=cbPosteriorBeta)) + 
    geom_point(shape=21, color="black", size=0.1) + 
    ggtitle(label = paste0("Correlation of Betas: ", names(mashPriorsAndPosteriors)[i])) +
    theme_classic() + 
    labs(x="DD + Corr Posterior Beta", y=" DD + CB Posterior Beta", fill="")+
    geom_hline(yintercept=0) + 
    geom_vline(xintercept=0) + 
    geom_smooth(method = "lm", se=FALSE) +
    theme(axis.text=element_text(size=15.5),
          axis.title=element_text(size=15.5)) + 
    theme(legend.text = element_text(size=15.5))
  
  ggsave(filename = i_file_name, plot = p, width=4.5, height=3)
  
}




