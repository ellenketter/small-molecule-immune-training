#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=diff
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=47
#SBATCH --mem=100GB

# exclude samples identified as outliers in limma analysis: HC3, Fluni3
module load python
cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/motifMatching

rgt-hint differential --organism=mm10 --bc --nc 47 \
--mpbs-files=./5F_mpbs.bed,./BG_mpbs.bed,./Fen_mpbs.bed,./Fluni_mpbs.bed,./HC_mpbs.bed,./HQ_mpbs.bed,./Myr_mpbs.bed,./Nerol_mpbs.bed,./PBS_mpbs.bed \
--reads-files=../5F.bam,../BG.bam,../Fen.bam,../Fluni.bam,../HC.bam,../HQ.bam,../Myr.bam,../Nerol.bam,../PBS.bam \
--conditions=5F,BG,Fen,Fluni,HC,HQ,Myr,Nerol,PBS \
--output-location=../diffFootprinting/multi/