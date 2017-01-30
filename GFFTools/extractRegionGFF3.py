from Bio import SeqIO
import re

fastaOBJ = open("Trichoderma_reesei_RUTC30.fasta","r")

for seq_record in SeqIO.parse(fastaOBJ, "fasta"):
 seqIdentifier = seq_record.id
 seqItself = seq_record.seq
 gffOBJ = open("TrireRUTC30_1.filtered_proteins.FilteredModels1.gff3","r")
 for line in gffOBJ:
  regex1 = re.compile(r"^#")
  if regex1.search(line):
   continue
  else:
   scaff,prediction,feattype,geneinit,geneend,frame,strand,misc,notes=line.split("\t")
   notes=notes.rstrip("\n")
   if (feattype == 'mRNA') and (scaff == seqIdentifier):
    #print("%s\t%s\t%s\t%s\t%s" % (scaff,geneinit,geneend,strand,notes)).rstrip("\n")
    if strand == '+':
     #continue
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

fastaOBJ.close()
