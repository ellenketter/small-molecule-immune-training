#################
# Installation of packages necessary for bulk ATAC seq analysis
# 
# Ellen Ketter, Luis Barreiro Lab
# March 2023
#################


module load java/15.0.2
cd /project/lbarreiro/USERS/ellen/programs/

wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.12.1.zip
unzip fastqc_v0.12.1.zip 

export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/FastQC

# NGmerge #############################
module load zlib/1.2 gsl/2.7
git clone https://github.com/harvardinformatics/NGmerge
cd NGmerge/
make
export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/NGmerge

# Bowtie 2 ################################
module load zlib/1.2 gsl/2.7
cd /project/lbarreiro/USERS/ellen/programs/
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.5.1/bowtie2-2.5.1-source.zip
unzip bowtie2-2.5.1-source.zip 
cd bowtie2-2.5.1-source
make
export PATH=$PATH:/project/lbarreiro/USERS/ellen/programs/bowtie2-2.5.1

# Mus musculus ENSEMBL genome #####################
cd /project/lbarreiro/USERS/ellen/programs/
wget http://igenomes.illumina.com.s3-website-us-east-1.amazonaws.com/Mus_musculus/Ensembl/GRCm38/Mus_musculus_Ensembl_GRCm38.tar.gz

# install macs2 ##########################
module load python
pip install macs2