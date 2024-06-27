# small-molecule-immune-training
Code implemented to analyze bulk ATAC-sequencing of BMDMs treated with novel small molecule inducers of trained immunity from Knight et al. <i>PNAS</i>. (2024)

Visualizations of this analysis can be found in Figure 4, Supplementary Table 3, and Supplementary Figure 4.

Raw sequencing data is deposited in [GEO]{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE270608}


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

Perform basic ATACseq quality control, eg. percent mitochondrial DNA and transcription start site enrichment with R package <i>ATACseqQC</i>

`6_removeChrM`

Remove reads aligned to mitochondrial sequences with <i>Samtools</i>

`7_peakCalling`

Call peaks from background with <i>macs2</i>

`8_consensusPeaks`

Remove blacklisted regions and define consensus peaks with reference file as replicate with greatest number of peaks using <i>bedtools intersect -wa -a</i>


# 3 - Differential Accessibility Modeling
Differential accessibility between each stimulation condition and the untrained control (PBS) was first modeled with limma and then secondly with mashr to increase sensitivity, as described in the methods of Knight et al (2024)

## limma
`9_limmaDA`

## mashr
`10_mashrDA`

# 4 - Gene Overenrichment Analysis
`11_mashrGeneORA`

Using the <i>ora</i> function in the <i>fgsea</i> R package, compare two lists of genes: one associated with peaks that were differentially accessible by mashr and one list of all genes associated with any consensus peaks identified.

# 5 - Differential Regulation of Transcription Factors

## Motif Enrichment
`12_TFmotifEnrichment`



## Footprinting
`13_footprinting`

Using HINT-ATAC from Regulatory Genomics Toolbox, determine transcription factor activity and differential regulation between trained and untrained conditions.


# 6- Scripts used to generate figure panels

## Figure Four

### B - 
Barplot summary of significant differentially accessible peaks determined by mashr (pAdj < 0.1)

`10_mashrDA/10-2_mashrVolcanoPlots.R`

### C - 
Heatmap of shared differentially accessible peaks by mashr (pAdj <0.1)

`10_mashrDA/10-0_mashr_dataDriven_correlations.R`

### D - 
Dot plot of overrepresentation analysis for genes near differentially accessible peaks in reactome pathways (pAdj < 0.1)

`11_mashrGeneORA/11-2-fgsea_ora_visualization.R`

### E - 
Dot plot of footprinting analysis for differentially regulated transcription factors (pAdj < 0.1)

`13_footprinting/13-9_dotPlot.R`

## Supplementary Table One -
Adjusted p values and mean betas from mashr analysis for all consensus peaks per treatment

`10_mashrDA/10-0_mashr_dataDriven_correlations.R`


## Supplementary Table Two -
Adjusted p values for two-list over-representation analysis of significantly upregulated genes present in REACTOME database pathways versus a background of all genes proximal to peaks measured in this ATACseq experiment, per treatment condition.

`11_mashrGeneORA/11-1_fgsea_ora.R`


## Supplementary Table Three -
Fold changes, adjusted p values, and JASPAR motif logos associated with all significant differentially regulated transcription factors identified by footprinting

`13_footprinting/13-10_significantMotifs.R`


## Supplementary Figure Four

### A - 
Volcano plots of significant differentially accessible peaks determined by mashr (pAdj > 0.1)

`10_mashrDA/10-2_mashrVolcanoPlots.R`

### B - 
Principal component analysis of mashr-computed betas

`10_mashrDA/10-0_mashr_dataDriven_correlations.R`

### C - 
Pathway over-enrichment of overlaps among differentially acessible peaks by mashr between treatment condtions 

`11_mashrGeneORA/11-3_oraOverlaps/11-3-3_overlaps_ora_visualization.R`

### D - 
Summary of differentially regulated transcription factors (pAdj < 0.1)

`13_footprinting/13-5_footprintingResults.R`

### E - 
Lineplots from HINT-ATAC footprinting output with custom themed color for treatment conditions

`13_footprinting/13-6_bespokeLineplots.R`

### F - 
Dot plot of transcription factor motif enrichment analysis by HOMER

`12_TFmotifEnrichment/12-3_dotPlot.R`

