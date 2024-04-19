#!/bin/bash
#install bedtools
module load python/anaconda-2020.11
cd /project/lbarreiro/USERS/ellen/programs
wget https://github.com/arq5x/bedtools2/releases/download/v2.29.1/bedtools-2.29.1.tar.gz
tar -zxvf bedtools-2.29.1.tar.gz
cd bedtools2/
make

#set directories
cd /project/lbarreiro/USERS/ellen/KnightMolecules/analysis/9_consensusPeaks/
wget https://github.com/Boyle-Lab/Blacklist/blob/master/lists/mm10-blacklist.v2.bed.gz