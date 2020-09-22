#!/bin/sh
#Modified script from Istvan Alberts biostars example

#$1 is the original bam file
#$2 is the prefix for split bam files

#example command:
# ~/ag-swart-blepharisma/data/alignments/intron_searches/samtools_extract_stranded_bams.sh merged.bam merged

set -ue

PROCS=16

# Forward strand.
#
# 1. alignments of the second in pair if they map to the forward strand
# 2. alignments of the first in pair if they map to the reverse  strand
#
samtools view -@ $PROCS -b -f 128 -F 16 $1 > $2.fwd1.bam
samtools index -@ $PROCS $2.fwd1.bam

samtools view -@ $PROCS -b -f 80 $1 > $2.fwd2.bam
samtools index -@ $PROCS $2.fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -@ $PROCS -f $2.fwd.bam $2.fwd1.bam $2.fwd2.bam
samtools index -@ $PROCS $2.fwd.bam

# Reverse strand
#
# 1. alignments of the second in pair if they map to the reverse strand
# 2. alignments of the first in pair if they map to the forward strand
#
samtools view -@ $PROCS -b -f 144 $1 > $2.rev1.bam
samtools index -@ $PROCS $2.rev1.bam

samtools view -@ $PROCS -b -f 64 -F 16 $1 > $2.rev2.bam
samtools index -@ $PROCS $2.rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -@ $PROCS -f $2.rev.bam $2.rev1.bam $2.rev2.bam
samtools index -@ $PROCS $2.rev.bam

# Cleanup

rm $2.fwd1.bam $2.fwd2.bam $2.rev1.bam $2.rev2.bam
rm $2.fwd1.bam.bai $2.fwd2.bam.bai $2.rev1.bam.bai $2.rev2.bam.bai
