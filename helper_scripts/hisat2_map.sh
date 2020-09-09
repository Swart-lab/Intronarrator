#!/bin/bash

#with unfiltered reads

FDIR=/ebio/abt2_projects/ag-swart-blepharisma
ASM=assembly.fixed.minus_mt

hisat2-build -f $ASM.fa $ASM.idx

FILE_COUNT=25
for (( i=2; i<=$FILE_COUNT; i++ ))
do

  hisat2 --threads 48 -q -x $ASM.idx --min-intronlen 9 --max-intronlen 500 \
  -1 $FDIR/raw/20190819_BGI_RNAseq_data/F19FTSEUET0126_BLEdkcE/clean_fastq/${i}_f.fq.gz \
  -2 $FDIR/raw/20190819_BGI_RNAseq_data/F19FTSEUET0126_BLEdkcE/clean_fastq/${i}_r.fq.gz \
  -S ${ASM}_${i}.sam

  samtools view -@ 48 -bS ${ASM}_${i}.sam | samtools sort -o ${ASM}_${i}.sorted.bam
  samtools view -@ 48 -h -F 4 -b ${ASM}_${i}.sorted.bam > ${ASM}_${i}.mapped_only.bam
  samtools index ${ASM}_${i}.mapped_only.bam


  rm ${ASM}_${i}.sam
done
