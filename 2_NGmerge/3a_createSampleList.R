setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/")
samples <- list.files()
samples <- samples[-c(1,56,59,60)]
samples <- gsub("_R1_001.fastq.gz","",samples)
samples <- gsub("_R2_001.fastq.gz","",samples)
samples <- unique(samples)
length(samples)
samples <- paste(as.character(samples), sep="' '", collapse="','")
