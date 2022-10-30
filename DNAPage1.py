#!/usr/bin/env python3

# Import modules for CGI handling 
import cgi, cgitb 
cgitb.enable()
import collections
import os
import argparse
import sys
import re
# from Bio import SeqIO 
# from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
# from Bio import AlignIO
# from Bio.PDB import *
# from Bio.PDB.Polypeptide import PPBuilder
import subprocess
import urllib
from urllib.request import urlopen
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import itertools
import shutil
import pandas as pd
# from Bio import motifs 
# from Bio.Seq import Seq 

import random
import string

sys.path.append('/usr/local/Anaconda3/lib/')
url="https://www.rcsb.org/structure/"


# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from field
if form.getvalue('textcontent'):
   text_content = form.getvalue('textcontent')

else:
   text_content = None
   
print ("Content-type:text/html\r\n\r\n")
print ("<html>")
print ("<head>")

#HTML style details for web page
print ("<style>")
print ("ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}")
print ("li{float:left;}")
print ("li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}")
print ("li a:hover { background-color: #111111;}")
print ("table, th, td { border: 2px solid black;}")
print (".footer { position: absolute; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }")
print ("</style>")
#Style ends here

print ("<title>MSALigMap</title>")
print ("</head>")
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
print ("<div align='center'>")
print ("</div>")
print ("<body>")
print ("<ul>")
print ("<li><a href='HomePage.py'>Home</a></li>")
print ("<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>")
print ("<li><a href='ProteinDNA.py'>Protein-DNA</a></li>")
print ("<li><a href='ProteinPeptide.py'>Protein-Peptide</a></li>")
print ("<li><a href='Contact.py'>Contact</a></li>")
print ("</ul>")
print ("<div align='center'>")
print ("<h2> Protein Details! </h2>")#print lig_content


if text_content != None:
    for i in text_content:
        text_content1=text_content.replace(' ', '')
        l=text_content1.split(',')
    

   

    foldernamer= ''.join(random.choices(string.ascii_letters, k=4))
    ProtLigfolder= '/opt/lampp/htdocs/MSALigMap/tmp/' + 'ProtDNA'

    isExist = os.path.exists(ProtLigfolder)

    if not isExist:
        os.mkdir( ProtLigfolder)

    folderpath= ProtLigfolder + '/' +foldernamer 
    os.mkdir(folderpath)


    fileitem = form['filename']
    fileattached = fileitem.value

    InputFileName = folderpath + '/sequenceInputfile.fasta'
    with open(InputFileName, 'wb') as fout:
        fout.write(fileattached)
    
    pdbid_list=[]
    url_list=[]
    prot_lig_dict=collections.defaultdict(dict)

    pdbid_chain_dict={}
    
    for i in l:
        
        each= i#.split(':')
        pdbcode= each[0]
        chain=each[1] # i am removing the +"="+"[]"
        pdbid_list.append(each)
        pdbid_chain_dict[pdbcode]=chain
    

    #linking the PDB url address to the selected PDB ids
    for x in l:
        pdbcodename= x.split(':')[0]
        pdbid=url+pdbcodename#+".pdb"
        url_list.append(pdbid)

    

    
    #Listing the PDB ids and their bound ligands
    Protein_DNAComplex_dict={}
    for link,id  in zip(url_list,pdbid_list):#list(pdbid_chain_dict.keys()

        PDB_id = id #+ '_' + Chain
        Protein_DNAComplex_dict[PDB_id]={}

        page = urlopen(link)
        html = page.read().decode("utf-8")
        soup = BeautifulSoup(html, "html.parser")

        divs = soup.find(id="MacromoleculesButton_NucleicAcids")

        
        divs1 = soup.find(id="MacromoleculeTable")

        
        divs2 = soup.find(id="MacromoleculeTableDNA")

        
        # Defining of the dataframe
        df = pd.DataFrame(columns=['PDB_Code','Molecule', 'Chain', 'Length'])
        for row in divs1.tbody.find_all('tr', id="macromolecule-entityId-1-rowDescription"):    
            # Find all data for each column
            columns = row.find_all('td')
            #print (columns)
            if(columns != []):
                molecule= columns[0].text.strip()
                chain=columns[1].text.strip()
                
                Protein_DNAComplex_dict[PDB_id][chain]={}
                length= columns[2].text.strip()
                Protein_DNAComplex_dict[PDB_id][chain]['Molecule']= "Protein"
                Protein_DNAComplex_dict[PDB_id][chain]['Detail']= molecule
                Protein_DNAComplex_dict[PDB_id][chain]['Length']= length
                


                df = df.append({'PDB_Code':PDB_id,'Molecule':molecule, 'Chain':chain, 'Length':length}, ignore_index=True)
        


        
        for tr in soup.find_all('table',class_='table table-bordered table-condensed tableEntity',id=re.compile("^table_macromolecule-dna-entityId")):
  
            columns2 = tr.find_all('td')
            if(columns2 != []):
                moleculeDNA= columns2[0].text.strip()
                chainDNA=columns2[1].text.strip()
                Protein_DNAComplex_dict[PDB_id][chainDNA]={}
                lengthDNA= columns2[2].text.strip()
                
                Protein_DNAComplex_dict[PDB_id][chainDNA]['Molecule']="DNA"
                Protein_DNAComplex_dict[PDB_id][chainDNA]['Detail']= moleculeDNA
                Protein_DNAComplex_dict[PDB_id][chainDNA]['Length']= lengthDNA
                
                df = df.append({'PDB_Code':PDB_id,'Molecule':moleculeDNA, 'Chain':chainDNA, 'Length':lengthDNA}, ignore_index=True)
                
      

        # Form for selecting the Peptide chains of interest for the USERR
        print ("<form enctype='multipart/form-data' action='LigPageDNA.py' method = 'post' target = '_blank'>")
        
        print ("<table style=width:50%>")
        print ("<tr>")
        print ("<th>PDB ID</th>")
        print ("<th >Chain</th>")
        print ("<th >Molecule</th>")
        
        print ("<th >Length</th>")
        print ("<th >Select</th>")

        for k,dk in Protein_DNAComplex_dict.items():
            length_of_rowspan = len(Protein_DNAComplex_dict[k])
            print ("<tr>")
            print ("<th rowspan='%d'>"% length_of_rowspan,k, "</th>")

            for chainKey, lengthValue in dk.items():
                print ("<td align=center>", chainKey,"</td>")
                print ("<td align=center>", dk[chainKey]['Detail'],"</td>")
                print ("<td align=center>", dk[chainKey]['Length'],"</td>")
                print ("<td align=center>")
                values_add = k + ':' + chainKey 
                if dk[chainKey]['Molecule']=='Protein':
                        
                    print  ("<input type='checkbox' name='DNASelection' value='%s'/>" % (values_add))
                print ("</td>")
                print ("</tr>")
                print ("<br/>")
                print ("<br/>")
        print ("</table>")
        print ("<br/>")

        print ("<input type = 'submit'  value = 'Submit'  />")
        print ("<input type = 'reset'  value = 'Clear'  />")
        print ("<input type='hidden' name='SequencePath' value='%s'>" %(InputFileName))
        
        print ("</form>")

else:
    print ("<h3>Please enter a PDB code and try again!</h3>")
    sys.exit()














print ("</body>")
print ("</html>")