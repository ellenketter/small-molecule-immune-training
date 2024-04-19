# library(BiocManager)
# BiocManager::install(c("ATACseqQC", "ChIPpeakAnno", "MotifDb", "GenomicAlignments"))
# BiocManager::install(c("BSgenome.Mmusculus.UCSC.mm10", "TxDb.Mmusculus.UCSC.mm10.knownGene"))
### # BiocManager::install(c("phastCons100way.UCSC.mm10"))????

library(ATACseqQC)
## input the bamFile from the ATACseqQC package 
setwd("/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/5_removeDuplicates/merged/")
source(system.file("extdata", "IGVSnapshot.R", package = "ATACseqQC"))

## bamfile tags to be read in
possibleTag <- list("integer"=c("AM", "AS", "CM", "CP", "FI", "H0", "H1", "H2", 
                                "HI", "IH", "MQ", "NH", "NM", "OP", "PQ", "SM",
                                "TC", "UQ"), 
                    "character"=c("BC", "BQ", "BZ", "CB", "CC", "CO", "CQ", "CR",
                                  "CS", "CT", "CY", "E2", "FS", "LB", "MC", "MD",
                                  "MI", "OA", "OC", "OQ", "OX", "PG", "PT", "PU",
                                  "Q2", "QT", "QX", "R2", "RG", "RX", "SA", "TS",
                                  "U2"))

# samtools reheader -c 'perl -pe "s/^(@SQ.*)(\tSN:)(\d+|X|Y|MT)(\s|\$)/\$1chr\$2\$3/"' 5F1_REP1.mLb.clN.sorted.bam > 5F1.sorted.bam
library(Rsamtools)
bamTop100 <- scanBam(BamFile("Fen1.chr1.bam", yieldSize = 100),
                     param = ScanBamParam(tag=unlist(possibleTag)))[[1]]$tag
tags <- names(bamTop100)[lengths(bamTop100)>0]
tags
## files will be output into outPath
outPath <- "splited"
# dir.create(outPath)
## shift the coordinates of 5'ends of alignments in the bam file
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
seqinformation <- seqinfo(TxDb.Mmusculus.UCSC.mm10.knownGene)
which <- as(seqinformation, "GRanges")
gal <- readBamFile("Fen1.chr1.bam", tag=tags, which=which, asMates=TRUE, bigFile=TRUE)
shiftedBamfile <- file.path(outPath, "shifted.bam")
gal1 <- shiftGAlignmentsList(gal, outbam=shiftedBamfile)
tsse <- TSSEscore(gal1, txs)
tsse$TSSEscore
