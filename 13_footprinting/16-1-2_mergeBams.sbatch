#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=caslake
#SBATCH --job-name=merge
#SBATCH --time=36:00:00


cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged/
module load samtools
for molecule in 5F Fen Myr Nerol
do
samtools merge --threads 20 -o ../../16_footprinting/mergedBams/${molecule}.bam ${molecule}1.bam ${molecule}2.bam ${molecule}3.bam 
done 
