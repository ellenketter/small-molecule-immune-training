#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=index
#SBATCH --time=5:00:00


cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/6_removeChrM

module load samtools
echo $1
samtools index -b $1.sorted.bam
