import argparse
from Bio import SeqIO
import re

__author__      = "Renato Augusto Correa dos Santos"
__copyright__   = "Copyright 2017"
__license__     = "GPL v3.0"
__email__       = "renatoacsantos@gmail.com"

version=sys.argv[0] + 'v.1'
parser = argparse.ArgumentParser(description='Extract regions from FASTA and GFF3 (annotation)', add_help=True)
parser.add_argument('-gff','--gff', dest='GFF3', metavar='genome.gff3', type=str, help='Annotation file (GFF3)', required=True)
parser.add_argument('-g','--genome', dest='INFASTA', metavar='genome.fasta', type=str, help='Chromosomal sequences (FASTA)', required=True)
parser.add_argument('-o','--out', dest='OUTFASTA', metavar='out.fasta', type=str, help='Output gene upstream regions (FASTA)', required=True)

args = parser.parse_args()
inFastaOBJ = args.INFASTA
outFastaOBJ = args.OUTFASTA

for seq_record in SeqIO.parse(inFastaOBJ, "fasta"):
 seqIdentifier = seq_record.id
 seqItself = seq_record.seq
 inGFFOBJ = open(args.GFF3,"r")
 for line in inGFFOBJ:
  regex1 = re.compile(r"^#")
  if regex1.search(line):
   continue
  else:
   scaff,prediction,feattype,geneinit,geneend,frame,strand,misc,notes=line.split("\t")
   notes=notes.rstrip("\n")
   if (feattype == 'mRNA') and (scaff == seqIdentifier):
    if strand == '+':
     if 0 <= (int(geneinit)-2)-1500 < len(seqItself):
      gene_sequps=seqItself[(int(geneinit)-1)-1500:(int(geneinit)-1)]
      lenseq=len(gene_sequps)
      print(">%s %s\n%s" % (notes,lenseq,gene_sequps))
     else:
      gene_sequps=seqItself[0:(int(geneinit)-1)]
      lenseq=len(gene_sequps)
      print(">%s %s\n%s" % (notes,lenseq,gene_sequps))
    elif strand == '-':
     if 0 <= (int(geneend)+1)+1500 < len(seqItself):
      gene_sequps=seqItself[int(geneend):(int(geneend))+1500]
      gene_sequps_revcomp = gene_sequps.reverse_complement()
      lenseq=len(gene_sequps_revcomp)
      print(">%s %s\n%s" % (notes,lenseq,gene_sequps_revcomp))
     else:
      gene_sequps=seqItself[int(geneend):-1]
      gene_sequps_revcomp = gene_sequps.reverse_complement()
      lenseq=len(gene_sequps_revcomp)
      print(">%s %s\n%s" % (notes,lenseq,gene_sequps_revcomp))
    else:
     print("Something is wrong.")
     exit(1)
 gffOBJ.close()

inFastaOBJ.close()
