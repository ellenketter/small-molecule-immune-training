#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=removeM
#SBATCH --ntasks=5 
#SBATCH --mem=100GB


echo $1

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/5_removeDuplicates/merged/

module load samtools
samtools sort -@ 5 $1.bam -o $1.sorted.bam
samtools index -@ 5 $1.sorted.bam
samtools view -h -b -q 10 -@ 5 $1.sorted.bam chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chrX chrY > ../../7_removeChrM/merged/$1.bam  

samtools sort -@ 5 ../../7_removeChrM/merged/$1.bam -o ../../7_removeChrM/merged/$1.sorted.bam 
samtools coverage ../../7_removeChrM/merged/$1.sorted.bam -o ../../7_removeChrM/merged/$1_coverage.out