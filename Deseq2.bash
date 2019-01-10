#!/bin/bash

#SBATCH --job-name=DESeq2
##SBATCH --nodes=1
##SBATCH --cpus-per-task=1
#SBATCH --mem=100GB
##SBATCH --gres=gpu:1
##SBATCH --partition=gpu4_medium
#SBATCH --partition=fn_long
#SBATCH --error=/gpfs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/rnaseq/err_out/%x_%j.err
#SBATCH --output=/gpfs/data/proteomics/projects/Sunny/jingchuan/2018-07-24_syni/rnaseq/err_out/%x_%j.out
##SBATCH --dependency=afterany:job_id

module load r/3.5.0

# Total of 6 inputs
# 1. CountTable (including path)
# 2. PhenoTable (including path)
# 3. Strain name of the Control
# 4. pvalue cutoff
# 5. log2FC cutoff
# 6. output directory

countTable=$1
phenoTable=$2
ctr=$3
pvalue=$4
log2fc=$5
outdir=$6

mkdir ${outdir}
Rscript /gpfs/data/proteomics/projects/Sunny/R/Deseq2_bash.R ${countTable} ${phenoTable} ${ctr} ${pvalue} ${log2fc} ${outdir}
