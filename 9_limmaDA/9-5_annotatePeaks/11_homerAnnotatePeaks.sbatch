#!/bin/bash

#SBATCH --account=pi-lbarreiro
#SBATCH --partition=caslake
#SBATCH --job-name=homer
#SBATCH --mem=40GB

module load samtools

PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/homer/bin
INPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/10_differentialAccessibility/background_vs_DA_peaks/

OUTPUT_PATH=/project/lbarreiro/USERS/ellen/KnightMolecules/demultiplexed/FastQ/11_annotatePeaks

/project/lbarreiro/USERS/ellen/programs/homer/bin/annotatePeaks.pl $INPUT_PATH/$1_genePeaksUp_p_0.05.txt mm10 -go $OUTPUT_PATH/$1_up_GO/ > $OUTPUT_PATH/$1_up_annotate_out.txt
/project/lbarreiro/USERS/ellen/programs/homer/bin/annotatePeaks.pl $INPUT_PATH/$1_genePeaksDown_p_0.05.txt mm10 -go $OUTPUT_PATH/$1_down_GO/ > $OUTPUT_PATH/$1_down_annotate_out.txt