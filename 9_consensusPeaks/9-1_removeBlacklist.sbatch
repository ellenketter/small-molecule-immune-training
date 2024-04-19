#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=rmBl

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/bedtools2/bin/
blacklist=/project/lbarreiro/USERS/ellen/KnightMolecules/analysis/9_consensusPeaks/mm10-blacklist.v2.sorted.bed
OUT_DIR=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks
peakCalls=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/8_peakCalling/merged

#remove blacklisted regions (only report peaks in a with no overlap in b)
bedtools intersect -a $peakCalls/$1_peaks.narrowPeak -b $blacklist -v -sorted > $OUT_DIR/$1_noBlacklist.bed