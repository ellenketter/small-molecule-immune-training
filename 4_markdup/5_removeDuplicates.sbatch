#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=markdup
#SBATCH --ntasks=10
#SBATCH --time=40:00:00
#SBATCH --mem=100GB


echo $1

#cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/4_bams/merged/

module load samtools
#samtools sort -@ 10 -n $1.sorted.bam -o  $1.namesorted.bam
#samtools fixmate -rm --threads=10 $1.namesorted.bam ../../5_removeDuplicates/merged/$1_fixmate.bam

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/5_removeDuplicates/merged/

#samtools sort -@ 10 -o $1_fixmate.sorted.bam $1_fixmate.bam
samtools markdup -rs -d 2500 -f $1_markdupStats.txt  --threads=10 $1_fixmate.sorted.bam $1.bam