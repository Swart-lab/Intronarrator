#!/bin/bash

#Example script for mapping reads from multiple fastq files - replace/alter as appropriate

FDIR= #replace with directory with FASTQ files
ASM= #replace with assembly name
THREADS = 48 #number of threads for hisat2
FILE_COUNT = 25 #number of files to process

hisat2-build -f $ASM.fa $ASM.idx

for (( i=2; i<=$FILE_COUNT; i++ ))
do

  hisat2 --threads $THREADS -q -x $ASM.idx --min-intronlen 9 --max-intronlen 500 \
  -1 $FDIR/${i}_f.fq.gz \
  -2 $FDIR/${i}_r.fq.gz \
  -S ${ASM}_${i}.sam

  samtools view -@ $THREADS -bS ${ASM}_${i}.sam | samtools sort -o ${ASM}_${i}.sorted.bam
  samtools view -@ $THREADS -h -F 4 -b ${ASM}_${i}.sorted.bam > ${ASM}_${i}.mapped_only.bam
  samtools index ${ASM}_${i}.mapped_only.bam

  rm ${ASM}_${i}.sam
done
