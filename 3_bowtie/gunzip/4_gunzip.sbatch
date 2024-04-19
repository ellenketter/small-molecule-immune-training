#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=bigmem
#SBATCH --job-name=gunzip2



cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/2_NGmergeOutput
# zcat AEK-RK-1s-Fen2_S22_1.fastq.gz 
zcat AEK-RK-1s-Fen2_S22_2.fastq.gz 
