#!/bin/bash                                                                                                                                                                                     
#$ -S /bin/bash                                                                                                                                                                                 
#$ -cwd                                                                                                                                                                                         

module load star/2.4.5a

# reads file                                                                                                                                                                                    
fq_r1=$1_R1_001.fastq.gz.cleaned.fastq.gz
fq_r2=$1_R2_001.fastq.gz.cleaned.fastq.gz

# folder of fastq files                                                                                                                                                                         
dir=/ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq
genomedir=/ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/star_index_2.4.5a
outputdir=/ifs/data/proteomics/projects/Sunny/rnaseq/clean_fastq/star_out/$1

# go to the output directory                                                                                                                                                                    
mkdir $outputdir
cd $outputdir

# Splice junction detection                                                                                                                                                                     
STAR \
--genomeDir ${genomedir} \
--readFilesIn ${dir}/${fq_r1} ${dir}/${fq_r2} \
--runThreadN 8 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMtype None \
--outSAMmode None \
--readFilesCommand zcat \

# INTERMEDIATE INDEX GENERATION                                                                                                                                                                 
STAR \
--runMode genomeGenerate \
--genomeDir ${outputdir} \
--genomeFastaFiles ${genomedir}/hg38_with_l1.fa \
--sjdbOverhang 100 \
--runThreadN 16 \
--sjdbFileChrStartEnd ${outputdir}/SJ.out.tab

# Alignment 2ND PASS                                                                                                                                                                                    
STAR \
--genomeDir $outputdir \
--readFilesIn ${dir}/${fq_r1} ${dir}/${fq_r2} \
--runThreadN 8 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 20 \
--outFilterMismatchNmax 10 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 0 \
--readFilesCommand zcat \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI NM MD AS XS \
--outSAMunmapped Within \
--outSAMtype BAM Unsorted \
--outSAMheaderHD @HD VN:1.4



module unload samtools/1.2.1
module load samtools/1.3

#sort bam file                     
samtools sort -n ${outputdir}/Aligned.out.bam > ${outputdir}/Aligned.out.bam.sort.bam

# htseq
                                 
#samtools view -F 4 ${outputdir}/Aligned.out.bam.sort.bam | \
#htseq-count \
#-i gene_id \
#-s no \
#-m intersection-nonempty \
#- /ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/star_index_2.4.5a/genes_with_l1.gtf \
#> ${outputdir}/$1_count.txt

# samtools filtering out the aligned reads                                                       
#samtools view -F 4 ${outputdir}/Aligned.sortedByCoord.out.bam | \
#python /ifs/home/xs338/HTSeq-0.6.1/HTSeq/scripts/count.py \
#-m intersection-nonempty \
#-i gene_id \
#-s no \
#-r pos \
#- /ifs/data/proteomics/projects/Sunny/rnaseq/rnaseq_genome/star_index/genes_with_l1.gtf