#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=counts
#SBATCH --mem=200GB

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/bedtools2/bin/
OUT_DIR=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks
peakCalls=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks/noBlacklist

#generate consensus peaks
bedtools intersect -wa -a $peakCalls/BG3_noBlacklist.bed -b $peakCalls/5F1_noBlacklist.bed $peakCalls/5F2_noBlacklist.bed $peakCalls/5F3_noBlacklist.bed $peakCalls/BG1_noBlacklist.bed $peakCalls/Fen1_noBlacklist.bed $peakCalls/Fen2_noBlacklist.bed $peakCalls/Fen3_noBlacklist.bed $peakCalls/Fluni1_noBlacklist.bed $peakCalls/Fluni2_noBlacklist.bed $peakCalls/Fluni3_noBlacklist.bed $peakCalls/HC1_noBlacklist.bed $peakCalls/HC2_noBlacklist.bed $peakCalls/HC3_noBlacklist.bed $peakCalls/HQ1_noBlacklist.bed $peakCalls/HQ2_noBlacklist.bed $peakCalls/Myr1_noBlacklist.bed $peakCalls/Myr2_noBlacklist.bed $peakCalls/Myr3_noBlacklist.bed $peakCalls/Nerol1_noBlacklist.bed $peakCalls/Nerol2_noBlacklist.bed $peakCalls/Nerol3_noBlacklist.bed $peakCalls/PBS1_noBlacklist.bed $peakCalls/PBS2_noBlacklist.bed $peakCalls/PBS3_noBlacklist.bed $peakCalls/PBSa1_noBlacklist.bed $peakCalls/PBSa2_noBlacklist.bed $peakCalls/PBSa3_noBlacklist.bed -sorted > $OUT_DIR/consensusPeaks.bed