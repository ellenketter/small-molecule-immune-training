#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=caslake
#SBATCH --job-name=merge
#SBATCH --time=36:00:00

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged/
module load samtools

samtools merge --threads 20 -o ../../16_footprinting/mergedBams/PBS.bam PBS1.bam PBS2.bam PBS3.bam PBSa1.bam PBSa2.bam PBSa3.bam
