# small-molecule-immune-training
Code implemented to analyze bulk ATAC-sequencing of BMDMs treated with novel small molecules from Knight et al. PNAS. (2024)

Visualizations of this analysis can be found in Figure 4, Supplementary Table 3, and Supplementary Figure 4.

Raw sequencing data is deposited in Geo. 


# 1 - Initial Installations
All ATAC-seq analysis was performed on the high performance research computing cluster Midway3 at the University of Chicago. To preprocess the raw sequencing, we began with command line installations found in `0-installations`. 

For subsequent analysis, we primarily applied the R programming language using version 4.1.0. Occasional implementation of python scripts or pip based installations used python-anaconda-2022.05-el8-x86_64. 

# 2 - Preprocessing
`1_mergeFastQ`
Merge FastQ files from the two sequencing runs

`2_NGmerge`
Merge paired-end reads and remove sequencing adapters using <i>NGmerge</i>

`3_alignment`
Align to mouse genome mm10 with <i>bowtie2</i>

`4_markdup`
Remove PCR duplicates with <i>Samtools markdup</i>

`5_atacSeqQC`
Perform basic ATACseq quality control, eg. determine percent mitochondrial DNA and transcription start site enrichment with R package <i>ATACseqQC</i>

`6_removeChrM`
Remove reads aligned to mitochondrial sequences with <i>Samtools</i>

`7_peakCalling`
Call peaks from background with <i>macs2</i>

`8_consensusPeaks`
Remove blacklisted regions and call consensus peaks with beta glucan #3 as reference file (greatest number of peaks) using <i>bedtools intersect -wa -a</i>


# 3 - Differential Accessibility Modeling


## limma

## mashr


# 4 - Gene Overenrichment Analysis


# 5 - Differential Regulation of Transcription Factors

## Motif Enrichment

## Footprinting


# Scripts used to generate each figure panel

## Figure Four

### B

### C

### D

### E

## Supplementary Table One

## Supplementary Table Two

## Supplementary Table Three

## Supplementary Figure Four

### A

### B

### C

### D

### E

### F
