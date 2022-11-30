#!/usr/bin/env python3

# Import modules for CGI handling 
import cgi, cgitb 
cgitb.enable()
import os
import argparse
import sys
# from Bio import SeqIO 
# from Bio.Align.Applications import MafftCommandline
#from StringIO import StringIO
# from Bio import AlignIO
# from Bio.PDB import *
# from Bio.PDB.Polypeptide import PPBuilder
import subprocess
import urllib
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import itertools
import shutil
# from Bio import motifs 
# from Bio.Seq import Seq 

import random
import string

sys.path.append('/usr/local/Anaconda3/lib/')
url="https://files.rcsb.org/view/"


# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from field
if form.getvalue('textcontent'):
   text_content = form.getvalue('textcontent')

else:
   text_content = ""
   
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
print ("<h2> Select a Ligand of interest! </h2>")#print lig_content

##Sequence file
#print lig_content
#capturing the entered pdb ids into list
if any(form.getvalue('filename').decode("utf-8")) == False  :
      print ("<h1>Error!</h1>")
      print ("<h3>Please upload a sequence file and try again!.</h3>")
      sys.exit()

if text_content != None:
    try:
        for i in text_content:
                text_content1=text_content.replace(' ', '')
                l=text_content1.split(',')
       
  
   

        foldernamer= ''.join(random.choices(string.ascii_letters, k=4))
        ProtLigfolder= '/opt/lampp/htdocs/MSALigMap/tmp/' + 'ProtLig'


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

      



        ##
        f2_list=[]
        combined_list=[]
        #print "entered id's",l
        pdbid_list=[]
        url_list=[]
        prot_lig_dict={}
        count=1

        
        pdbid_chain_dict={}
        for i in l:
            
            each= i.split(':')
            
            pdbcode= each[0].upper()
            chain=each[1].upper() # i am removing the +"="+"[]"
            pdbid_list.append(each)
            
            pdbid_chain_dict[pdbcode]=chain
        

     


        #linking the PDB url address to the selected PDB ids
        for x in l:
            pdbcodename= x.split(':')[0]
            pdbid=url+pdbcodename+".pdb"
            url_list.append(pdbid)

        
        


        #Listing the PDB ids and their bound ligands
        
        for link,id  in zip(url_list,list(pdbid_chain_dict.keys())):
                
                count=count+1
                try:
                    
                    f=urllib.request.urlopen(link)
                    f=f.readlines()
                    #print (f)
                    PDBCode_Chains=[]
                    for li in f:
                        if li.startswith(b'HET '):
                            
                            li= li.decode('UTF-8').split()
                            
                            
                            PDBCode_Chains.append(li[2])
                            
                            if li[2]==pdbid_chain_dict[id]:
                                f2=li[1]
                                prot_lig_dict.setdefault(id,[]).append(f2)
                            
                    if  not pdbid_chain_dict[id] in PDBCode_Chains:
                        print ("<h3>Error! Please enter the correct chain for the PDB code and try again!</h3>")
                        sys.exit()
                        
                                
                            
                            
                except:
                    if pdbid_chain_dict=={}:

                        print("<p>Error! Check the PDB code and the chain. Try again!</p>")

        ## Checking 
        



        # Form for selecting the ligands for interest for the USER
        if not prot_lig_dict =={}:

            print ("<form enctype='multipart/form-data' action='LigPage.py' method = 'post' target = '_blank'>")
            print ("<table style=width:50%>")
            print ("<tr>")
            print ("<th>PDB ID</th>")
            print ("<th colspan=2>LIGANDS</th>")

            


            for k,dk in prot_lig_dict.items():

                count= len(dk)
                print ("<tr>")
                print ("<th rowspan='%d'>"% count,k,":",pdbid_chain_dict[k], "</th>")
                #print "</tr>"
                for x in dk:
                    value_code = k + ':' +  pdbid_chain_dict[k]
                    print ("<td align=center>", "%s" % x,"</td>")
                    print ("<td align=center>")
                    values_add = value_code + '_' + x
                    print  ("<input type='radio' name='LigSelection' value='%s'/>" % (values_add))
                    
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
    
    except:
        print ("<br/>")
        print ("<h3><b>ERROR!! Please enter PDB code with chain and try again!</b></h3>")


print ("</body>")
print ("</html>")