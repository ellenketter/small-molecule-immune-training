#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=macs
#SBATCH --mem=100GB

module load python/anaconda-2022.05
cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/

echo $1

macs2 callpeak -f BAMPE -t $1.bam -g mm -n ./$1 -B --keep-dup all
