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



# Per peak, show mean and SD of prior vs posterior beta
# Picked peaks which are DA in all conditions, according to MASH
peaksScarb1 <- readRDS("../13_mashr/peaksScarb1.rds")

allIntersect <- paste(peaksScarb1$Chr,peaksScarb1$Start,peaksScarb1$End, sep = ".")
# final digit off in start: e.g.125332936 should be 125332935
allIntersect <- c("chr5.125339059.125341372")

# c("chr5.125332935.125333471","chr5.125348186.125348767","chr5.125339059.125341372",
#                 "chr5.125272636.125272895","chr5.125357960.125358737","chr5.125271661.125271954",
#                 "chr5.125273602.125274389","chr5.125329216.125329640","chr5.125351684.125352221",
#                 "chr5.125305232.125305684","chr5.125291846.125292521","chr5.125313191.125314505",
#                 "chr5.125278440.125278717")


mashPriorsAndPosteriors <- inner_join(mashInput,mashOutput)

ppPlot <- list()
i=1
# for(i in 1:13){

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
# }
dev.off()

svg(paste0("ppPlot_Scarb1.svg"))
marrangeGrob(ppPlot, nrow = 2, ncol = 2)
dev.off()

mashInput$coefficients[mashInput$genes == "chr17.84956870.84958199"]
topTable(mashInput)[2,]




chr3 144718854 144719722

