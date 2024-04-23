#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=caslake
#SBATCH --job-name=merge
#SBATCH --time=36:00:00


cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged/
module load samtools
for twos in Fluni HC HQ
do
samtools merge --threads 20 -o ../../16_footprinting/mergedBams/${twos}.bam ${twos}1.bam ${twos}2.bam
done
