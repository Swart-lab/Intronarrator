#!/usr/bin/env python
import sys, pysam
from re import compile
from Bio import SeqIO
import multiprocessing
import shutil
from glob import glob
import argparse
import os.path

parser = argparse.ArgumentParser()
parser.add_argument('--processes', action='store', type=int, default=8)
parser.add_argument('--bam', action='store', type=str, required=True)
parser.add_argument('--intron_juncs', action='store', type=str, required=True)
parser.add_argument('--genome', action='store', type=str, required=True)
parser.add_argument('--outfile', action='store', type=str, required=True)
parser.add_argument('--tmpdir', action='store', type=str, default='tmp',
  help='Folder for temporary files')

args = parser.parse_args()

#strict_intron_p = compile("^[0-9][0-9]M([0-9][0-9])N[0-9][0-9]M$")
intron_p = compile("[0-9]+M([0-9]+)N[0-9]+M")

ref_recs = SeqIO.parse(open(args.genome), "fasta")

contigs = [rec.name for rec in ref_recs]

def process_intron(samfile, contig, start, end, strand, junc):
  intron_spliced = 0
  read_tot = 0
  fwd_strand_tot = 0
  rev_strand_tot = 0

  for read in samfile.fetch(contig, start, end):
    try:
      if strand == "+" and read.is_read1 and read.is_reverse:
        proper_orientation = True
      elif strand == "+" and read.is_read2 and not read.is_reverse:
        proper_orientation = True
      elif strand == "-" and read.is_read1 and not read.is_reverse:
        proper_orientation = True
      elif strand == "-" and read.is_read2 and read.is_reverse:
        proper_orientation = True
      else:
        proper_orientation = False

      if not read.is_secondary and read.mapping_quality >= 10 and read.cigarstring.count("N") < 2:
        #Do not count reads with more than one intron
        read_tot += 1

        if (read.is_read1 and read.is_reverse) or (read.is_read2 and not read.is_reverse):
          fwd_strand_tot += 1
        elif (read.is_read1 and not read.is_reverse) or (read.is_read2 and read.is_reverse):
          rev_strand_tot += 1

      if not read.is_secondary and proper_orientation and read.mapping_quality >= 10:
        #do not count secondary alignments, check that read orientation matches intron and that mapq >= 10

        cig_match = intron_p.search(read.cigarstring)

        if cig_match != None and len(read.cigartuples) == 3 and read.cigartuples[0][1] >= 6 and read.cigartuples[2][1] >= 6:
          #only count reads with 6 or more bases on either side matching
          #Do not count reads with more than one intron

          intron_start = read.reference_start + read.cigartuples[0][1]
          intron_end = read.reference_start + read.cigartuples[0][1] + read.cigartuples[1][1]

          if intron_start == start - 1 and intron_end == end:
            #Note: convert intron start coordinate to zero-based coordinate used by pysam
            intron_spliced += 1
          
          #Debugging code:
          #if read.query_name == 'V300020225L1C002R0581258290':
          #  intron_seq = ref_seq[intron_start:intron_end]
          #  sys.stderr.write("%s %s %s %s %s %s\n" % (read.query_name, intron_start, intron_end, start, end, intron_seq))

    except TypeError:
      pass

  return "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contig, start, end, strand, junc, intron_spliced, fwd_strand_tot, rev_strand_tot, read_tot)


def process_all_introns(juncs, contig, tmpdir):

  samfile = pysam.AlignmentFile(args.bam, "rb" )
  lines = []
  for atoms in juncs[0:100]:
    #lines.append(process_intron(samfile, atoms[0], int(atoms[1]), int(atoms[2]), atoms[5], int(atoms[4])))
    lines.append(process_intron(samfile, atoms[0], int(atoms[1]), int(atoms[2]), atoms[4], int(atoms[5])))

  print("processed: ", contig)

  samfile.close()

  f = open('%s/%s.txt' % (tmpdir, contig), 'w+')
  f.writelines(lines)
  f.close()


def read_intron_juncs(fn):

  junc_d = {}
  for line in open(fn):
    atoms = line.split() 
    junc_d.setdefault(atoms[0], []).append(atoms)

  pool = multiprocessing.Pool(processes=args.processes) 
  for contig in contigs:
    try:
      pool.apply_async(process_all_introns, args=(junc_d[contig], contig, args.tmpdir, ))
    except KeyError:
      pass

  pool.close()
  pool.join()

  fh = open(args.outfile, 'w+')
  fh.write("#contig	start	end	strand	juncs	spliced	fwd_strand_all	rev_strand_all\n")
  fh.close()

  with open(args.outfile, 'ab') as wfd:
    for f in glob("%s/*" % (args.tmpdir)):
        with open(f,'rb') as fd:
            shutil.copyfileobj(fd, wfd)

#  print("contig	start	end	strand	junc	intron_spliced	intronless	total")
#  for line in open(fn):
#    atoms = line.split() 
#    if atoms[0] == 'contig_83':
#    process_intron(atoms[0], int(atoms[1]), int(atoms[2]), atoms[5], int(atoms[4]))

if __name__ == "__main__":
  if os.path.isfile(args.outfile):
    print("output file %s already exists, not overwriting" % (args.outfile))
  else:
    if os.path.exists(args.tmpdir):
      if os.path.isdir(args.tmpdir):
        read_intron_juncs(args.intron_juncs)
      else:
        print("temp path %s is not a folder" % (args.tmpdir))
    else:
      print("temp folder %s does not exist, please mkdir" % (args.tmpdir))


# vim:sts=2:ts=2:sw=2
