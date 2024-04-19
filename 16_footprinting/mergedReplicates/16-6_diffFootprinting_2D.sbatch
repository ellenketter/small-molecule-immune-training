#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=diff
#SBATCH --time=15:00:00
#SBATCH --cpus-per-task=47
#SBATCH --mem=100GB

# excluded samples identified as outliers in limma analysis: HC3, Fluni3
module load python
cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/motifMatching
mkdir ../diffFootprinting/$1

rgt-hint differential --organism=mm10 --bc --nc 47 \
--mpbs-files=./$1_mpbs.bed,./PBS_mpbs.bed \
--reads-files=../$1.bam,../PBS.bam \
--conditions=$1,PBS \
--output-location=../diffFootprinting/$1
