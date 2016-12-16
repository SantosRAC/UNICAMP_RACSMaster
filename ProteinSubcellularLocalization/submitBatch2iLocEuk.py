#!/home/renato/Software/Anaconda/anaconda2/bin/python

# Using Python 2.7.12
# Should set PATH to geckodriver

from selenium import webdriver
from selenium.webdriver.common.keys import Keys
import time
import re
from Bio import SeqIO
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.common.exceptions import TimeoutException

for seq_record in SeqIO.parse("/home/renato/Projects/CNPEM/Pseudozyma/Annotation/AnnotationCuration_03112016/Proteome03162016/KbrasiliensisGHG001_evm.out.all.mod.proteome.mod.fasta", "fasta"):
 seqIdentifier = seq_record.id
 seqItself = seq_record.seq
 # Open Firefox
 driver = webdriver.Firefox()
 # Go to iLoc-Euk page, fill form, then submit request
 driver.get("http://www.jci-bioinfo.cn/iLoc-Euk")
 time.sleep(3)
 assert "iLic-Euk" in driver.title
 form = driver.find_element_by_id("seq")
 form.clear()
 form.send_keys(">%s\n%s" % (seqIdentifier,seqItself))
 time.sleep(3)
 #<input type="button" id="sbmitBtn" value="Submit">
 sendButton = driver.find_element_by_id("sbmitBtn")
 sendButton.send_keys(Keys.RETURN)
 assert "No results found." not in driver.page_source
 # Waiting time to process result
 time.sleep(120)
 # Try to find element on the page (driver source code)
 results = driver.find_element_by_id('resultDetail')
 results_text = results.text
 #print("%s\n" % results_text)
 regex = re.compile(r"Predicted Result: (\w+) \( Predicted By PSS\)")
 for res in regex.findall(results_text):
  print("%s\t%s" % (seqIdentifier,res))
 driver.quit()
