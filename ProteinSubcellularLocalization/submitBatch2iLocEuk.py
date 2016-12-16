#!/home/renato/Software/Anaconda/anaconda2/bin/python

# Using Python 2.7.12
# Should set PATH to geckodriver

import argparse
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import re
from Bio import SeqIO
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

parser = argparse.ArgumentParser(description='Run iLoc-Euk for a set of proteins')
parser.add_argument('-o','--out', dest='out', metavar='file.tab', type=str, help='Output with proteins and corresponding locations', required=True)
parser.add_argument('-f','--fasta', dest='fasta', metavar='proteome.fasta', type=str, help='FASTA with proteome', required=True)
parser.add_argument('-t','--time', dest='wTime', metavar='N', type=int, help='Integer with waiting time for each request (it depends on the internet connect and server responsinevess)', required=True)

args = parser.parse_args()
fastaOBJ = open(args.fasta,"r")
outOBJ = open(args.out,"w")
waitingTime = args.wTime

for seq_record in SeqIO.parse(fastaOBJ, "fasta"):
 seqIdentifier = seq_record.id
 seqItself = seq_record.seq
 # Open Firefox
 driver = webdriver.Firefox()
 # Go to iLoc-Euk page, fill form, then submit request
 driver.get("http://www.jci-bioinfo.cn/iLoc-Euk")
 time.sleep(3)
 form = driver.find_element_by_id("seq")
 form.clear()
 form.send_keys(">%s\n%s" % (seqIdentifier,seqItself))
 time.sleep(3)
 #<input type="button" id="sbmitBtn" value="Submit">
 sendButton = driver.find_element_by_id("sbmitBtn")
 sendButton.send_keys(Keys.RETURN)
 assert "No results found." not in driver.page_source
 # Waiting time to process result
 time.sleep(waitingTime)
 # Try to find element on the page (driver source code)
 results = driver.find_element_by_id('resultDetail')
 results_text = results.text
 #print("%s\n" % results_text)
 regex = re.compile(r"Predicted Result: (\w+) \( Predicted By PSS\)")
 for res in regex.findall(results_text):
  outOBJ.write("%s\t%s" % (seqIdentifier,res))
  outOBJ.write("\n")
 driver.quit()

outOBJ.close()
