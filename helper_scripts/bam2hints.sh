#!/bin/bash
#exonpart hints
AUGUSTUS_PATH=$1
BAM2WIG=$AUGUSTUS_PATH/bin/bam2wig
WIG2HINTS=$AUGUSTUS_PATH/scripts/wig2hints.pl
BAM_PREFIX=$2
HINT_PREFIX=$3

FWD_BAM=$BAM_PREFIX.fwd.bam
REV_BAM=$BAM_PREFIX.rev.bam

$BAM2WIG $FWD_BAM > $BAM_PREFIX.fwd.wig

$WIG2HINTS --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --pri=2 \
--strand=+ --radius=4.5 < $BAM_PREFIX.fwd.wig > $HINT_PREFIX.fwd.ep.gff

$BAM2WIG $REV_BAM > $BAM_PREFIX.rev.wig

$WIG2HINTS --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --width=10 --margin=10 --minthresh=2 --minscore=4 --prune=0.1 --src=W --type=ep --pri=2 \
--strand=- --radius=4.5 < $BAM_PREFIX.rev.wig > $HINT_PREFIX.rev.ep.gff

