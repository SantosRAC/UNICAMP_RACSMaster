#!/home/renato/Software/Anaconda/anaconda2/bin/python

# Using Python 2.7.12
# Should set PATH to geckodriver

import argparse
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import re
from Bio import SeqIO
import os.path
import os
import sys

parser = argparse.ArgumentParser(description='Run Euk-mPLoc2 for a set of proteins')
parser.add_argument('-o','--out', dest='out', metavar='file.tab', type=str, help='Output with proteins and corresponding locations', required=True)
parser.add_argument('-f','--fasta', dest='fasta', metavar='proteome.fasta', type=str, help='FASTA with proteome', required=True)
parser.add_argument('-t','--time', dest='wTime', metavar='N', type=int, help='Integer with waiting time for each request (it depends on the internet connect and server responsinevess)', required=True)
parser.add_argument('-s','--sound', dest='sound', action='store_true', help='Make sound if the script finished', required=False)

listOfProteinsWithSubLoc=[]

args = parser.parse_args()
fastaOBJ = open(args.fasta,"r")
waitingTime = args.wTime
makeSound=False

# Control sound
if args.sound:
 makeSound=args.sound

if os.path.isfile(args.out):
 outOBJ = open(args.out,"r")
 for line in outOBJ:
  a,b=line.split("\t")
  #print("%s\t%s" % (a,b)).rstrip("\n")
  if not a in listOfProteinsWithSubLoc:
   listOfProteinsWithSubLoc.append(a)
 outOBJ.close()

driver = webdriver.Firefox()

for seq_record in SeqIO.parse(fastaOBJ, "fasta"):
 seqIdentifier = seq_record.id
 seqItself = seq_record.seq
 print("Checking protein %s" % (seqIdentifier))
 if seqIdentifier in listOfProteinsWithSubLoc:
  print("Protein %s is already with predicted location" % (seqIdentifier))
  continue
 # Open Firefox
 print("Analizing protein %s" % (seqIdentifier))
 # Go to iLoc-Euk page, fill form, then submit request
 driver.get("http://www.csbio.sjtu.edu.cn/bioinf/euk-multi-2/")
 time.sleep(10)
 form = driver.find_element_by_tag_name("textarea")
 form.clear()
 form.send_keys(">%s\n%s" % (seqIdentifier,seqItself))
 time.sleep(10)
 #<input type="button" id="sbmitBtn" value="Submit">
 sendButton = driver.find_element_by_name("B1")
 sendButton.send_keys(Keys.RETURN)
 assert "No results found." not in driver.page_source
 # Waiting time to process result
 time.sleep(waitingTime)
 # Try to find element on the page (driver source code)
 results = driver.find_element_by_tag_name('table')
 results_text = results.text
 print(results_text)
 regex1 = re.compile(r"evm.model.(.+)\.(\d+)\.(\d+) (.+)")
 if regex1.search(results_text):
  for res in regex1.findall(results_text):
   outOBJ = open(args.out,"a")
   outOBJ.write("%s\t%s" % (seqIdentifier,res[3:]))
   outOBJ.write("\n")
   outOBJ.close()
   listOfProteinsWithSubLoc.append(seqIdentifier)
 else:
  exit(1)
  outOBJ = open(args.out,"a")
  outOBJ.write("%s\tIMPOSSIBLE_COMPLETE_PROBLEM" % (seqIdentifier))
  outOBJ.write("\n")
  outOBJ.close()
  listOfProteinsWithSubLoc.append(seqIdentifier)
  #print("Something weird happened. Check online result for protein %s" % seqIdentifier)
  if makeSound:
   os.system("/usr/bin/canberra-gtk-play --id='bell' --loop=10") 
  #sys.exit()
 if seqIdentifier in listOfProteinsWithSubLoc:
  continue
 else:
  print("There might be something weird with the server - protein %s was not added to output file!" % (seqIdentifier))
  if makeSound:
   os.system("/usr/bin/canberra-gtk-play --id='bell' --loop=10")
  sys.exit()

driver.quit()
fastaOBJ.close()
if makeSound:
 os.system("/usr/bin/canberra-gtk-play --id='bell' --loop=10")
