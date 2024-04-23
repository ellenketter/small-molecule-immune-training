#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=lbarreiro
#SBATCH --job-name=homer
#SBATCH --mem=40GB

module load samtools

PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/homer/bin
INPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks

OUTPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/13_mashr

/project/lbarreiro/USERS/ellen/programs/homer/bin/annotatePeaks.pl $OUTPUT_PATH/annotateInput.txt mm10 -go $OUTPUT_PATH/annotateOutput_GO/ > $OUTPUT_PATH/annotateOutput.txt
