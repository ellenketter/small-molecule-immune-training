 footStats <- read.delim("Downloads/diffFootprinting/differential_statistics.txt")

outputs <- colnames(footStats)
moleculeN <- gsub("TC_","",outputs[29:54])
# calculate activity scores per replicate 
for(i in 3:28){
  footStats <- cbind(footStats, rowSums(footStats[c(i,i+26)]))}
colnames(footStats)[55:80] <- moleculeN
rownames(footStats) <- footStats[,1]
footStats <- subset(footStats, Num > 1000) # removed 5 rows

`5F1` <- footStats[c(55:57,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

BG <- footStats[c(58:59,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

Fen <- footStats[c(60:62,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

Fluni <- footStats[c(63:64,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

HC <- footStats[c(65:66,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

LFC_5F1 <- log2(rowMeans(`5F1`[1:3])/rowMeans(`5F1`[4:9]))
z <- 
  
  
  HQ <- footStats[c(67:68,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

Myr <- footStats[c(69:71,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

Nerol <- footStats[c(72:74,75:80)]
`5F1`$pValT <- sapply(1:nrow(`5F1`), function(i) t.test(as.numeric(as.character(unlist(`5F1`[i,1:3]))), 
                                                        as.numeric(as.character(unlist(`5F1`[i,4:9]))))[c("p.value")])
`5F1`$pAdjT<- p.adjust(`5F1`$pValT) # none significant

tTestRow<- function(x) {}
footStats$pVal <- sapply(1:nrow(footStats), function(i) t.test(as.numeric(as.character(unlist(footStats[i,55:57]))), 
                                                               as.numeric(as.character(unlist(footStats[i,75:80]))))[c("p.value")])

