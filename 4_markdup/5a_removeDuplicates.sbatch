#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=caslake
#SBATCH --job-name=markdup
#SBATCH --ntasks=10

cd /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/6_removeDuplicates/

module load samtools
echo AEK-RK-1s-PBSa3_S28
samtools sort -@ 10 -o AEK-RK-1s-PBSa3_S28_fixmate.sorted.bam AEK-RK-1s-PBSa3_S28_fixmate.bam
samtools markdup -rs -d 2500  --threads=10 AEK-RK-1s-PBSa3_S28_fixmate.sorted.bam AEK-RK-1s-PBSa3_S28.bam