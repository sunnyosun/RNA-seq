#!/bin/bash                                                                                                                                                     
#$ -S /bin/bash                                                                                                                                                 
#$ -cwd

module load java/1.8
module load jre/1.8
module load bowtie2/2.3.1
module load tophat/2.1.1
module load cufflinks/2.2.1

# genome basename of the .fa file
genome=/ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/l1hs_bt2/L1HS

# annotation file -basename of the gtf file
gtf=/ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/l1hs_bt2/l1.gtf

star_index_l1hs=/ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/l1hs_index
dir=/ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1
dir2=/ifs/data/proteomics/projects/Sunny/rnaseq/maprrcon
pathtoreg=/ifs/data/proteomics/projects/Sunny/RepeatMaskerL1HsInfo
reg=ucsc.rmsk.hg38.L1HS.bed
threshmis=$2

sortbam=$1.sorted.bam
l1hsbam=$sortbam.L1HsOnly.bam

#L1HS feature extraction
cd $dir2
java -classpath sam-1.89.jar:. ExtractSubsetOfRNAseqBAMBasedOnBedFile $dir $pathtoreg $sortbam $reg $threshmis

cd $dir
# convert bam to fastq
module unload samtools/0.1.19
module load samtools/1.3
samtools fastq $dir/$l1hsbam > $dir/$l1hsbam.fastq

l1hsdir=$dir/l1hs
mkdir $l1hsdir

cd $l1hsdir
# tophat alignment -> -N 2 allowing 2 mismatches in the final output reads
tophat -N 2 --GTF ${gtf} --num-threads 8 --output-dir ${l1hsdir} ${genome} ${l1hsbam}.fastq

# sort bam file by name -n
samtools sort -n ${l1hsdir}/accepted_hits.bam ${l1hsdir}/$1_l1hs.sort

# convert bam to sam and output read count using htseq
# samtools view /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/accepted_hits.bam | htseq-count - ${gtf} > /ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/$1/$1_count.txt

# bai file
samtools index ${l1hsdir}/$1_l1hs.sort.bam

# htseq-count
htseq-count --format=bam --stranded=no ${l1hsdir}/$1_l1hs.sort.bam ${gtf} > ${l1hsdir}/$1_l1hs_count.txt

