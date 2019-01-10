#!/bin/bash
#$ -S /bin/bash
#$ -cwd

#SBATCH --job-name=genomeIndex
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=cpu_long
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/genome/err_out/%x_%j.err 
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/genome/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

#module load bwa
#bwa index /ifs/data/proteomics/projects/Sunny/genome/hg38.fa

#module load bowtie2
#genome=$1
#bowtie2-build /ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/${genome}.fa ${genome}

# 1- .fa
# 2- .gtf
# 3- index directory

genome=$1
index_dir=`pwd`/$1_star

mkdir ${index_dir}
module load star/2.6.0a

STAR --runMode genomeGenerate --runThreadN 8 --genomeDir ${index_dir} --genomeFastaFiles ${genome}.fa --sjdbGTFfile ${genome}.gtf --outFileNamePrefix ${genome} --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfeatureExon CDS --genomeSAindexNbases 10

