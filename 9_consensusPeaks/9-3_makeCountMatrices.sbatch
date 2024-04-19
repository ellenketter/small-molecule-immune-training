#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=counts
#SBATCH --mem=300GB

#fen2 required 400GB on bigmem
export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/bedtools2/bin/
OUT_DIR=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/9_consensusPeaks
bamFiles=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged/

#make count tables for each .bam file
echo $1
bedtools intersect -a $OUT_DIR/consensusPeaks.bed -b $bamFiles/$1.bam -c > $OUT_DIR/$1_countsRaw.txt
uniq $OUT_DIR/$1_countsRaw.txt $OUT_DIR/$1_noDup_counts.txt
