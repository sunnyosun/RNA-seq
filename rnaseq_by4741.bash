#!/bin/bash       
#$ -S /bin/bash    
#$ -cwd       

module load bowtie2/2.3.1
module load tophat/2.1.1
module load cufflinks/2.2.1
module load samtools/1.3


# input 1 is the basename of fastq file
# input 2 is the directory contains fastq files


# genome basename of the .fa file                                                                                                                                                    
genome=/ifs/data/proteomics/projects/Sunny/genome/BY4741_genome

# annotation file -basename of the gtf file                                                                                                                                          
gtf=/ifs/data/proteomics/projects/Sunny/genome/BY4741_genome.gff

dir=$2

fq1=$1_R1_001.fastq.gz
fq2=$1_R2_001.fastq.gz
sortbam=$1.sorted.bam
#mkdir ${dir}/$1
outdir=$1

# tophat alignment -> -N 2 allowing 2 mismatches in the final output reads
tophat -N 2 --GTF ${gtf} --num-threads 8 --output-dir ${dir}/${outdir} ${genome} ${dir}/${fq1} ${dir}/${fq2}

# sort bam file by name -n
# samtools sort -n /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/accepted_hits.bam > /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/${sortbam}

# convert bam to sam and output read count using htseq       
# samtools view /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/accepted_hits.bam | htseq-count - ${gtf} > /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/$1_count.txt

# bai file       
# samtools index /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/accepted_hits.bam

# htseq-count       
# htseq-count --format=bam --stranded=no ${outdir}/${sortbam}.bam ${gtf} > ${outdir}/$1_count.txt
