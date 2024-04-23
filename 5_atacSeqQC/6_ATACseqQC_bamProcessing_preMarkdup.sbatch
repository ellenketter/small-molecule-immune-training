#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=qcPrep
#SBATCH --time=8:00:00


echo $1

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/4_bams/merged/

module load samtools
#samtools index $1.sorted.bam
#samtools view -h -b -q 10 $1.sorted.bam chr1 > ../../6_ATACseqQC/$1.preMarkdup.chr1.bam
samtools index ../../6_ATACseqQC/$1.preMarkdup.chr1.bam 