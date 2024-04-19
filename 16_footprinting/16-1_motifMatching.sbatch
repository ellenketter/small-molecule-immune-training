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
for celltype in 5F1 5F2 5F3 BG1 BG3 Fen1 Fen2 Fen3 Fluni1 Fluni2 Fluni3 HC1 \
HC2 HC3 HQ1 HQ2 Myr1 Myr2 Myr3 Nerol1 Nerol2 Nerol3 PBS1 PBS2 PBS3 PBSa1 PBSa2 PBSa3;
do
rgt-motifanalysis matching --organism=mm10 \
--output-location=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/motifMatching \
--input-files /project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/16_footprinting/${celltype}.bed
done

# background made of random regions of size 10 times the biggest genomic region in input.
