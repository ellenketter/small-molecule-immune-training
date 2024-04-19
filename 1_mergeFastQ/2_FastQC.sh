#################
# Run quality control analysis per replicate
# 
# Ellen Ketter, Luis Barreiro Lab
# March 2023
#################
module load java/15.0.2


cd /project/lbarreiro/USERS/ellen/programs/FastQC/
./fastqc -t 6 ../../KnightMolecules/demultiplexed/FastQ/1_originalFiles/SecondRun/*.gz
# -t parameter is number of threads/cores

mkdir ../KnightMolecules/demultiplexed/QualityControl/
mkdir ../KnightMolecules/demultiplexed/QualityControl/FastQC/
cd ../KnightMolecules/demultiplexed/QualityControl/FastQC/

mv *zip ../FastQC
mv ../FastQ/*html ./FastQC/