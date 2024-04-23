#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=macs
#SBATCH --mem=100GB

module load python/anaconda-2020.11
cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/7_removeChrM/merged

echo $1

macs2 callpeak -f BAMPE -t $1.sorted.bam -g mm -n ../../8_peakCalling/merged/$1 -B --keep-dup all