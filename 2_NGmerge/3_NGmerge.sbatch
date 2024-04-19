#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=ngm
#SBATCH --time=8:00:00
#SBATCH --ntasks=10

module load gsl/2.7 zlib/1.2

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/NGmerge

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/1_originalFiles/merged/

NGmerge -y -a -v -n 10 -1 $1_R1.fastq.gz  -2 $1_R2.fastq.gz  -o ../../2_NGmergeOutput/merged/$1
# -n is the number of threads

#################
# Run quality control analysis per replicate
# 
# Ellen Ketter, Luis Barreiro Lab
# March 2023
#################
