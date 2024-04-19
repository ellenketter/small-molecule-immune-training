#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=diff
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=48
#SBATCH --mem=100GB

# exclude samples identified as outliers in limma analysis: HC3, Fluni3
module load python
cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/motifMatching

rgt-hint differential --organism=mm10 --bc --nc 48 \
--mpbs-files=./5F1_mpbs.bed,./5F2_mpbs.bed,./5F3_mpbs.bed,./BG1_mpbs.bed,./BG3_mpbs.bed,./Fen1_mpbs.bed,./Fen2_mpbs.bed,./Fen3_mpbs.bed,./Fluni1_mpbs.bed,./Fluni2_mpbs.bed,./HC1_mpbs.bed,./HC2_mpbs.bed,./HQ1_mpbs.bed,./HQ2_mpbs.bed,./Myr1_mpbs.bed,./Myr2_mpbs.bed,./Myr3_mpbs.bed,./Nerol1_mpbs.bed,./Nerol2_mpbs.bed,./Nerol3_mpbs.bed,./PBS1_mpbs.bed,./PBS2_mpbs.bed,./PBS3_mpbs.bed,./PBSa1_mpbs.bed,./PBSa2_mpbs.bed,./PBSa3_mpbs.bed \
--reads-files=../../7_removeChrM/merged/5F1.bam,../../7_removeChrM/merged/5F2.bam,../../7_removeChrM/merged/5F3.bam,../../7_removeChrM/merged/BG1.bam,../../7_removeChrM/merged/BG3.bam,../../7_removeChrM/merged/Fen1.bam,../../7_removeChrM/merged/Fen2.bam,../../7_removeChrM/merged/Fen3.bam,../../7_removeChrM/merged/Fluni1.bam,../../7_removeChrM/merged/Fluni2.bam,../../7_removeChrM/merged/HC1.bam,../../7_removeChrM/merged/HC2.bam,../../7_removeChrM/merged/HQ1.bam,../../7_removeChrM/merged/HQ2.bam,../../7_removeChrM/merged/Myr1.bam,../../7_removeChrM/merged/Myr2.bam,../../7_removeChrM/merged/Myr3.bam,../../7_removeChrM/merged/Nerol1.bam,../../7_removeChrM/merged/Nerol2.bam,../../7_removeChrM/merged/Nerol3.bam,../../7_removeChrM/merged/PBS1.bam,../../7_removeChrM/merged/PBS2.bam,../../7_removeChrM/merged/PBS3.bam,../../7_removeChrM/merged/PBSa1.bam,../../7_removeChrM/merged/PBSa2.bam,../../7_removeChrM/merged/PBSa3.bam \
--conditions=5F1,5F2,5F3,BG1,BG3,Fen1,Fen2,Fen3,Fluni1,Fluni2,HC1,HC2,HQ1,HQ2,Myr1,Myr2,Myr3,Nerol1,Nerol2,Nerol3,PBS1,PBS2,PBS3,PBSa1,PBSa2,PBSa3 \
--output-location=../diffFootprinting
