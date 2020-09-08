#!/bin/bash

#add directory containing intronnarator scripts
INTRONARRATOR_PATH=/ebio/ag-swart/home/eswart/ag-swart-blepharisma/data/alignments/intron_searches/intronarrator
export PATH=$PATH:$INTRONARRATOR_PATH:$INTRONARRATOR_PATH/ncRNA_searches

# --- Intron splicing parameters:
MIN_INTRON_RATIO=0.2 #default 0.3 - at least 30% of mapped reads with introns
MIN_INTRONS=10 #default 10 - at least 10 reads with introns spliced out
MAX_INTRON_LEN=40 #default 40 - only find introns less than this length (bp)

ALL_PROCS=32 #number of processes for multiprocessing/threading
INFERNAL_PROCS=$ALL_PROCS #number of processes for Infernal searches of ncRNAs
PYSAM_PROCS=$ALL_PROCS #number of processes for intron splicing counting
BAM=assembly.fixed.minus_mt.mapped_only.bam
#/tmp/scratch/map/rna-seq/merged.intron_lines.sorted.bam

ASM=assembly.fixed.minus_mt #Genome assembly prefix without ".fa" extension

# --- AUGUSTUS binary, config and exon part hints paths/config
AUGUSTUS_CONFIG_PATH=/ebio/ag-swart/home/eswart/ag-swart-blepharisma/data/gene_prediction/augustus/augustus_mod_gencode2_ps_bjap_v2/config
AUGUSTUS_SCRIPTS=/ebio/abt2_projects/ag-swart-blepharisma/build/Augustus/scripts/
export AUGUSTUS_CONFIG_PATH
AUGUSTUS_BIN=/ebio/abt2_projects/ag-swart-blepharisma/build/Augustus/bin/augustus
#_shorty_2
AUGUSTUS_EXTRINSIC=$AUGUSTUS_CONFIG_PATH/species/bjap_v26/bjap_v26_extrinsic.M.RM.E.W.cfg

#The hints below for AUGUSTUS can be created with the scripts in "helper_scripts"
EP_hints=assembly.fixed.minus_mt.mapped_only.eph #AUGUSTUS hints file without ".gff" extension
SPLIT_FILE_SIZE=2000000 #Size of chunks to try to split genome into for parallel AUGUSTUS gene prediction
                        #This will automatically set the number of processes executed to the number of chunks
                        #So, a large enough number should be chosen to avoid forking too many processes
SPLIT_GENOME=split_genome
SPLIT_GFF=split_gff

JUNCS=$ASM.juncs

RFAM_DB=/ebio/ag-swart/home/eswart/ag-swart-blepharisma/db/Rfam/Rfam.cm

# --- Filename setup - It should not be necessary to change anything in this
# section.
REALTRONS_PREFIX=all.realtrons
REALTRONS_PREFIX2=$REALTRONS_PREFIX.${MIN_INTRON_RATIO}
REALTRONS_FN=$REALTRONS_PREFIX2.gff
REALTRONS_NO_ALT=all.realtrons.${MIN_INTRON_RATIO}.noalt.gff
ASM_M=$ASM.${MIN_INTRON_RATIO}
ASM_ncRNA_masked=$ASM_M.minus_introns.ncRNA_hard_masked
EP_hints_minus=$EP_hints.$MIN_INTRON_RATIO

# --- Directory setup - It should not be necessary to change anything in this
mkdir -p $SPLIT_GENOME
mkdir -p $SPLIT_GFF
mkdir -p tmp
mkdir -p tmp_juncs

# --- Code execution below ---

pysam_extract_introns_parallel.py --processes $PYSAM_PROCS --bam $BAM --genome $ASM.fa --outfile $JUNCS
sort -k 1,1 -k 2,2n $JUNCS > $JUNCS.sorted

pysam_count_reads.py --processes $PYSAM_PROCS --bam $BAM --intron_juncs $JUNCS.sorted --genome $ASM.fa --outfile $ASM.intron_counts.txt 
count both spliced reads and total reads properly 

realtrons.py $ASM.intron_counts.txt 0.01 $MIN_INTRONS $MAX_INTRON_LEN > $REALTRONS_PREFIX.0.01.gff
antisense_realtrons.py $ASM.intron_counts.txt 0.01 $MIN_INTRONS $MAX_INTRON_LEN > $REALTRONS_PREFIX.0.01.antisense.gff

realtrons.py $ASM.intron_counts.txt $MIN_INTRON_RATIO $MIN_INTRONS $MAX_INTRON_LEN > $REALTRONS_FN
antisense_realtrons.py $ASM.intron_counts.txt 0.01 $MIN_INTRONS $MAX_INTRON_LEN > $REALTRONS_PREFIX2.antisense.gff

realtrons_filter_alt.py $REALTRONS_FN > $REALTRONS_NO_ALT
#Selects only the most frequently spliced introns that share a boundary
realtrons_filter_alt.py $REALTRONS_PREFIX2.antisense.gff > $REALTRONS_NO_ALT.antisense

mask_introns.py $ASM.fa $REALTRONS_NO_ALT >  $ASM_M.intron_masked.fa
remove_introns.py $ASM_M.intron_masked.fa > $ASM_M.minus_introns.fa

convert_ep_hint_coords.py $ASM_M.intron_masked.fa $EP_hints.gff 1> $EP_hints_minus.gff 2> $EP_hints.gff.err

#ncRNA search with Infernal, followed by filtering to remove undesirable domains and weak matches:
cmsearch --cpu $INFERNAL_PROCS --tblout $ASM_M.minus_introns.cmsearch $RFAM_DB $ASM.minus_introns.fa > $ASM.minus_introns.cmsearch.aln
infernal_filter.py $ASM_M.minus_introns.cmsearch > $ASM_M.minus_introns.filtered.cmsearch
infernal_hard_masker.py $ASM_M.minus_introns.filtered.cmsearch $ASM_M.minus_introns.fa > $ASM_ncRNA_masked.fa

##Split genome for parallel processing with AUGUSTUS
$AUGUSTUS_SCRIPTS/splitMfasta.pl --minsize=$SPLIT_FILE_SIZE --outputpath=$SPLIT_GENOME $ASM_ncRNA_masked.fa
SPLIT_FILE_COUNT=`ls $SPLIT_GENOME | wc -l`

echo "$SPLIT_FILE_COUNT AUGUSTUS processes running in the background"
for (( i=1; i<=$SPLIT_FILE_COUNT; i++ ))
do
  $AUGUSTUS_BIN \
  --species=bjap_v26 --translation_table=4 \
  --genemodel=intronless \
  --hintsfile=$EP_hints_minus.gff \
  --extrinsicCfgFile=$AUGUSTUS_EXTRINSIC \
  --softmasking=1 \
  $SPLIT_GENOME/${ASM_ncRNA_masked}.split.$i.fa \
  > $SPLIT_GFF/${ASM_ncRNA_masked}.$i.gff &
  pids[${i}]=$!
done

FAIL=0
#Wait on execution of all the AUGUSTUS processes before continuing
for pid in ${pids[*]}; do
    wait $pid || let "FAIL+=1"
done

if [ "$FAIL" == "0" ];
then
  echo "AUGUSTUS jobs all completed"
else
  echo "($FAIL) AUGUSTUS jobs failed!"
fi

merge_and_fix_AUGUSTUS_gene_names.py "$SPLIT_GFF/*" > $ASM_ncRNA_masked.gff 

add_introns_back_to_AUGUSTUS_gff.py $ASM_M.intron_masked.fa $ASM_ncRNA_masked.gff $REALTRONS_NO_ALT > $ASM_M.final.gff 

##Output CDS and protein sequences
gff_extract_and_join_CDS.py $ASM.fa $ASM_ncRNA_masked.gff >  $ASM_M.final.CDS.fa
translate_CDS.py $ASM_M.final.CDS.fa > $ASM_M.final.CDS.pep

gff_extract_and_join_CDS.py $ASM_M.minus_introns.fa $ASM_ncRNA_masked.gff > $ASM_M.minus_introns.CDS.fa
translate_CDS.py $ASM_M.minus_introns.CDS.fa > $ASM_M.minus_introns.CDS.pep

##Find ncRNAs (including tRNAs) and output GFFs for them
cmsearch --cpu $INFERNAL_PROCS --tblout $ASM.cmsearch $RFAM_DB $ASM.fa > $ASM.cmsearch.aln
infernal_filter.py $ASM.cmsearch > $ASM.filtered.cmsearch
infernal_to_GFF.py $ASM.filtered.cmsearch> $ASM.filtered.gff
tRNAscan-SE --thread $INFERNAL_PROCS -E -o $ASM.tRNAscan-SE -f $ASM.tRNAscan-SE.aln $ASM.fa
tRNAscan-SE_to_GFF.py $ASM.tRNAscan-SE > $ASM.tRNAscan-SE.gff


# vim:sts=2:ts=2:sw=2:tw=300
