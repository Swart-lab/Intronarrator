#!/usr/bin/env python
import sys, pysam
from re import compile
from Bio import SeqIO
import multiprocessing
import shutil
from glob import glob
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--processes', action='store', type=int, default=8)
parser.add_argument('--bam', action='store', type=str, required=True)
parser.add_argument('--genome', action='store', type=str, required=True)
parser.add_argument('--outfile', action='store', type=str, required=True)

tmp_dir = 'tmp_juncs'

args = parser.parse_args()

intron_p = compile("[0-9]+M([0-9][0-9])N[0-9]+M")
ref_recs = SeqIO.parse(open(args.genome), "fasta")
contigs = [rec.name for rec in ref_recs]

def process_contig(contig):
  print("%s processing" % contig)
  ok_cigar_ops = set([0, 3])

  samfile = pysam.AlignmentFile(args.bam, "rb" )

  intron_d = {}
  for read in samfile.fetch(contig):
    if not read.is_secondary and read.mapping_quality >= 10 and read.cigarstring.count("N") == 1:
      #do not count secondary alignments, check that read orientation matches intron and that mapq >= 10

      cig_match = intron_p.search(read.cigarstring)
      cigar_ops = set([seg[0] for seg in read.cigartuples])
      
      if cigar_ops == ok_cigar_ops and cig_match != None and read.cigartuples[0][1] >= 6 and read.cigartuples[2][1] >= 6:

        if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse):
          strand = "+"
        elif (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):
          strand = "-"
        #only count reads with 6 or more bases on either side matching
        intron_start = read.reference_start + read.cigartuples[0][1] + 1
        intron_end = read.reference_start + read.cigartuples[0][1] + read.cigartuples[1][1]

        intron_d.setdefault((intron_start, intron_end, strand), 0)
        intron_d[(intron_start, intron_end, strand)] += 1


  samfile.close()

  intron_keys = list(intron_d.keys())
  intron_keys.sort()


  lines = []
  for intron in intron_keys:
    lines.append("%s\t%s\t%s\t.\t%s\t%s\n" % (contig, intron[0], intron[1], intron[2], intron_d[intron]))


  f = open('%s/%s.txt' % (tmp_dir, contig), 'w+')
  f.writelines(lines)
  f.close()

  print("%s processed" % contig)


def find_introns():

  pool = multiprocessing.Pool(processes=args.processes) 
  for contig in contigs:
    #print("starting processing %s" % contig)
    pool.apply_async(process_contig, args=(contig, ))

  pool.close()
  pool.join()

  fh = open(args.outfile, 'w+')
  fh.write("#contig	start	end	score	strand	count\n")
  fh.close()

  with open(args.outfile, 'wb') as wfd:
    for f in glob("%s/*" % tmp_dir):
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)

if __name__ == "__main__":
  find_introns()


# vim:sts=2:ts=2:sw=2
