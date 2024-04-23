#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=bow
#SBATCH --time=10:00:00
#SBATCH --ntasks=22
#SBATCH --mem=100GB

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/bowtie2-2.5.1

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/mm10

module load samtools
echo $1
bowtie2 -x GRCm38 -1 ../2_NGmergeOutput/$1_1.fastq -2 ../2_NGmergeOutput/$1_2.fastq -S ../2_NGmergeOutput/$1.sam -p 20  --very-sensitive
samtools view -u -S -b ../2_NGmergeOutput/$1.sam > ../2_NGmergeOutput/$1.bam
samtools sort ../2_NGmergeOutput/$1.bam -o  ../2_NGmergeOutput/$1.sorted.bam

#################
# Align to mouse index sequence GRCm38 (pre-built) per replicate
# 
# Ellen Ketter, Luis Barreiro Lab
# March 2023
#################
