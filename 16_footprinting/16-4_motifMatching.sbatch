#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=match
#SBATCH --time=36:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=10GB

module load python
for celltype in 5F BG Fen Fluni HC HQ Myr Nerol PBS;
do
rgt-motifanalysis matching --organism=mm10 \
--output-location=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/motifMatching \
--input-files /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/mergedBams/footprints/${celltype}.bed
done

# background made of random regions of size 10 times the biggest genomic region in input.
