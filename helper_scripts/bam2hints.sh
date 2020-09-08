#!/bin/bash
#exonpart hints
BAM2WIG=/ebio/abt2_projects/ag-swart-blepharisma/data/gene_prediction/augustus/augustus_mod_gencode2_ps_bjap_v2/bin/bam2wig 
WIG2HINTS=/ebio/abt2_projects/ag-swart-blepharisma/build/Augustus/scripts/wig2hints.pl 
BAM_PREFIX=assembly.fixed.minus_mt.mapped_only
FWD_BAM=$BAM_PREFIX.fwd.bam
REV_BAM=$BAM_PREFIX.rev.bam

#$BAM2WIG $FWD_BAM > $BAM_PREFIX.fwd.wig

$WIG2HINTS --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --pri=2 \
--strand=+ --radius=4.5 < $BAM_PREFIX.fwd.wig > $BAM_PREFIX.fwd.eph.gff

#$BAM2WIG $REV_BAM > $BAM_PREFIX.rev.wig

$WIG2HINTS --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --pri=2 \
--strand=- --radius=4.5 < $BAM_PREFIX.rev.wig > $BAM_PREFIX.rev.eph.gff

