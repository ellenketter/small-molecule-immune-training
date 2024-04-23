#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=motifs
#SBATCH --time=3:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=40GB

module load samtools


PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/homer/bin
INPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr/

OUTPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/15_transcriptionFactorMotifs/dataDriven_commonBaseline/

/project/lbarreiro/USERS/ellen/programs/homer/bin/findMotifsGenome.pl $OUTPUT_PATH/$1 mm10 \
$OUTPUT_PATH/$3/ -size 500 -mask -p 3 -bg $INPUT_PATH/$2