#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=foot
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB

module load python

rgt-hint footprinting --atac-seq --paired-end --organism=mm10 \
--output-location=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting \
--output-prefix=$1 \
/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged/$1.bam \
/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/8_peakCalling/merged/$1_peaks.narrowPeak

