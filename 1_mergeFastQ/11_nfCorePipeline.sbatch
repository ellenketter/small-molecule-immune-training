#!/bin/bash


#SBATCH --partition=lbarreiro
#SBATCH --account=pi-lbarreiro
#SBATCH --job-name=nfcore
#SBATCH --time=24:00:00
#SBATCH --ntasks=10
#SBATCH --mem=100GB

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/
export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/FastQC/
module load java
cd /project/lbarreiro/USERS/ellen/KnightMolecules/analysis/
nextflow run nf-core/atacseq --input samplesheet.csv --outdir 11_mergeFastQ/ --genome GRCm38 --read_length 50 -r 2.0