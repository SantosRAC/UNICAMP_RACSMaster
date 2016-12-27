#!/home/renato/Software/Anaconda/anaconda2/bin/python

from __future__ import print_function
import sys

#Making sure you are running a version of python that works with this script.
if sys.version_info[0] != 2:
    print("This script requires Python version 2")
    sys.exit(1)

import os
import os.path
import re

res = {}
genesInUMAmodel={}

if not len(sys.argv) > 1:
 raise Exception("You must provide input BLAST file as argument!")
#################################
## Get variables form user input
#################################
filename=sys.argv[1]
orthofile=sys.argv[2]
infoSBML=sys.argv[3]

##################################
## Check whether input file exists
##################################

if not os.path.exists(filename):
 raise Exception("File " + filename + " does not exists. Check your input!")

if not os.path.exists(orthofile):
 raise Exception("File: " + orthofile + " does not exists. Check your input!")

if not os.path.exists(infoSBML):
 raise Exception("File: " + infoSBML + " does not exists. Check your input!")
##################################
#Process blast file Ustilago NCBI vs Broad
##################################
file=open(filename,"r")
for pline in file.readlines():
 line=pline.strip()
 (sseqid, qseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore, qlen, slen)=line.split('\t')
 sseqid=sseqid.upper()
 #print("ID (subject): %s, qlen: %d,slen: %d,percAlign (subject) %f,p ident: %f" % (sseqid, int(qlen),int(slen),float(((float(send)-float(sstart))/float(slen))),float(pident)))
 if int(qlen) == int(slen) and int(length) == int(qlen) and float(pident)>=100:
  #print("1 %s %s %s %s %s" % (qseqid, sseqid, pident, length, qlen), end="\n")
  if sseqid in res:
   res[sseqid][qseqid]=1
  else:
   res[sseqid] = {qseqid:1}
 elif (int(length) == int(qlen) or int(length) == int(slen)) and float(pident)>=100:
  #print("2 %s %s %s %s %s %s" % (qseqid, sseqid, pident, length, qlen, slen), end="")
  if sseqid in res:
   res[sseqid][qseqid]=1
  else:
   res[sseqid] = {qseqid:1} 
 elif int(qlen) == int(slen) and int(length) >= 0.9*int(qlen) and float(pident)>=100:
#  print("3 %s %s %s %s %s %s" % (qseqid, sseqid, pident, length, qlen, slen), end="")
  if sseqid in res:
   res[sseqid][qseqid]=1
  else:
   res[sseqid] = {qseqid:1}
 elif int(qlen) > int(slen) and float(((float(send)-float(sstart))/float(slen))) >= 0.95 and float(pident) >= 100:
  if sseqid in res:
   res[sseqid][qseqid]=1
  else:
   res[sseqid] = {qseqid:1} 
 elif int(slen) > int(qlen) and float(((float(qend)-float(qstart))/float(qlen))) >= 0.95 and float(pident) >= 97:
  if sseqid in res:
   res[sseqid][qseqid]=1
  else:
   res[sseqid] = {qseqid:1}
 
file.close()

#for sid in res:
# print("%s\t" % (sid),end="")
# for qid in res[sid]:
#  print("%s\t" % (qid),end="")
# print("\n",end="")

##################################
# Process Ustilaginaceae orthologues
##################################
orthologues=open(orthofile,"r")

seq2group={}
group2ustilago={}

for pline in orthologues.readlines():
 line=pline.rstrip()
 (orthoid, pre_seqs)=line.split(': ')
 seqs=pre_seqs.split(' ')
 for seq in seqs:
  (spe, seqid)=str(seq).split('|')
  if spe == 'UMA1':
   seq2group[seqid]=orthoid
   if orthoid in group2ustilago:
    if spe in group2ustilago[orthoid]:
     group2ustilago[orthoid][spe][seqid]=1
    else:
     group2ustilago[orthoid][spe]={seqid:1}
   else:
    group2ustilago[orthoid]={spe:{seqid:1}}
  if spe == 'KBR1':
   if orthoid in group2ustilago:
    if spe in group2ustilago[orthoid]:
     group2ustilago[orthoid][spe][seqid]=1
    else:
     group2ustilago[orthoid][spe]={seqid:1}
   else:
    group2ustilago[orthoid]={spe:{seqid:1}}

#for oid in group2ustilago:
# for spid in group2ustilago[oid]:
#  for seqid in group2ustilago[oid][spid]:
#   if spid=='UMA1':
#    if seqid in res:
#     for i in res[seqid].keys():
#      print("%s %s %s\n" % (oid,spid,i), end="")
#    #else:
     #print("%s %s %s\n" % (oid,spid,seqid), end="")

orthologues.close()

infoSBMLfile=open(infoSBML,"r")

ec2gene={}
geneWithouEC=[]
genesNotFoundBLASTOufile='notFoundBlastUmaydis.txt';
genesNotFoundBLASTOufileOBJ=open(genesNotFoundBLASTOufile,"w")
orthologsNotFoundOufile='orthologsNotFound.txt';
orthologsNotFoundOufileOBJ=open(orthologsNotFoundOufile,"w")

for pline in infoSBMLfile.readlines():
 label=0
 line=pline.rstrip()
 mm = re.search('GENE_ASSOCIATION:',line)
 mmec = re.search('EC Number:',line)
 geneColletion={}
 if mm:
  fields=line.split('	')
  for i in fields:
   m = re.match('(GENE_ASSOCIATION: ?)(.+)',i)
   mec = re.match('(EC Number: ?)(.+)',i)
   if m:
    #print("GENE_ASSOCIATION:",end="")
    if re.search(' or ',m.group(2)):
     listOR=m.group(2).split(' or ')
     for g in listOR:
      g=g.upper()
      g=g.replace("(","")
      g=g.replace(")","")
      if not g in genesInUMAmodel.keys():
       genesInUMAmodel[g]=1
      if g in res.keys():
       ustncbiid,value=res[g].items()[0]
       if ustncbiid in seq2group.keys():
        orthoid=seq2group[ustncbiid]
        if 'KBR1' in group2ustilago[orthoid].keys():
         kalmanozymaID,value=group2ustilago[orthoid]['KBR1'].items()[0]
         #print(" (%s: %s)" % (g, kalmanozymaID),end="")
         if g in geneColletion.keys():
          if kalmanozymaID in geneColletion[g].keys():
           geneColletion[g][kalmanozymaID]=1
          else:
           geneColletion[g]={kalmanozymaID:1}
         else:
          geneColletion[g]={kalmanozymaID:1}
        else:
         orthologsNotFoundOufileOBJ.write("%s (PEDANT), %s (NCBI), %s (orthoMCL group)\n" % (g,ustncbiid,orthoid)) 
         #print(" (%s: no_gene_K_brasiliensis)" % (g),end="")
      else:
       genesNotFoundBLASTOufileOBJ.write("%s\n" % g)
    elif re.search(' and ',m.group(2)):
     listAND=m.group(2).split(' and ')
     for g in listAND:
      g=g.upper()
      g=g.replace("(","")
      g=g.replace(")","")
      if not g in genesInUMAmodel.keys():
       genesInUMAmodel[g]=1
      if g in res.keys():
       ustncbiid,value=res[g].items()[0]
       if ustncbiid in seq2group.keys():
        orthoid=seq2group[ustncbiid]
        if 'KBR1' in group2ustilago[orthoid].keys():
         kalmanozymaID,value=group2ustilago[orthoid]['KBR1'].items()[0]
         #print(" (%s: %s)" % (g, kalmanozymaID),end="")
         if g in geneColletion.keys():
          if kalmanozymaID in geneColletion[g].keys():
           geneColletion[g][kalmanozymaID]=1
          else:
           geneColletion[g]={kalmanozymaID:1}
         else:
          geneColletion[g]={kalmanozymaID:1}
        else:
         orthologsNotFoundOufileOBJ.write("%s (PEDANT), %s (NCBI), %s (orthoMCL group)\n" % (g,ustncbiid,orthoid))
         #print(" (%s: no_gene_K_brasiliensis)" % (g),end="")
      else:
       genesNotFoundBLASTOufileOBJ.write("%s\n" % g)
    else:
     g=m.group(2)
     g=g.upper()
     g=g.replace("(","")
     g=g.replace(")","")
     if not g in genesInUMAmodel.keys():
      genesInUMAmodel[g]=1
     if g in res.keys():
       ustncbiid,value=res[g].items()[0]
       if ustncbiid in seq2group.keys():
        orthoid=seq2group[ustncbiid]
        if 'KBR1' in group2ustilago[orthoid].keys():
         kalmanozymaID,value=group2ustilago[orthoid]['KBR1'].items()[0]
         #print(" (%s: %s)" % (g, kalmanozymaID),end="")
         if g in geneColletion.keys():
          if kalmanozymaID in geneColletion[g].keys():
           geneColletion[g][kalmanozymaID]=1
          else:
           geneColletion[g]={kalmanozymaID:1}
         else:
          geneColletion[g]={kalmanozymaID:1}
        else:
          orthologsNotFoundOufileOBJ.write("%s (PEDANT), %s (NCBI), %s (orthoMCL group)\n" % (g,ustncbiid,orthoid))
         #print(" (%s: no_gene_K_brasiliensis)" % (g), end="")
     else:
       genesNotFoundBLASTOufileOBJ.write("%s\n" % g)
   elif mec:
   # print("EC Number:",end="")
    listEC=re.split(r', |/',mec.group(2))
    for g in listEC:
     ec2gene[g]=geneColletion
   else:
    True
  if mmec:
   True
  else:
   geneWithouEC.append(geneColletion)

genesNotFoundBLASTOufileOBJ.close()
orthologsNotFoundOufileOBJ.close()
infoSBMLfile.close()

geneAssociationOufile='associationUmaydisKbrasiliensis.txt';
geneAssociationOufileOBJ=open(geneAssociationOufile,"w")
geneAssociationNoECsOufile='associationUmaydisKbrasiliensis.noECs.txt';
geneAssociationNoECsOufileOBJ=open(geneAssociationNoECsOufile,"w")

for ec in ec2gene.keys():
 for geneustilago in ec2gene[ec].keys():
  for genekalmanozyma in ec2gene[ec][geneustilago].keys():
   geneAssociationOufileOBJ.write("%s %s %s\n" % (ec,geneustilago,genekalmanozyma))
   #print("%s %s %s" % (ec,geneustilago,genekalmanozyma))

for genesNoEC in geneWithouEC:
 for geneustilago in genesNoEC.keys():
  absent=True
  for ec in ec2gene.keys():
   if geneustilago in ec2gene[ec].keys():
    absent=False
  for genekalmanozyma in genesNoEC[geneustilago].keys():
   if absent:
    geneAssociationNoECsOufileOBJ.write("%s %s\n" % (geneustilago,genekalmanozyma))

geneAssociationOufileOBJ.close()
geneAssociationNoECsOufileOBJ.close()

listGeneInModel='listGenesInModelUMA.txt'
listGeneInModelOBJ=open(listGeneInModel, "w")
for gene in genesInUMAmodel.keys():
 listGeneInModelOBJ.write("%s\n" % (gene))

listGeneInModelOBJ.close()
