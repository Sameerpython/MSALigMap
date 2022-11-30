#!/usr/bin/env python3

# Import modules for CGI handling 
import cgi, cgitb 
cgitb.enable()
import collections
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
print ("<h2> Select a Protein chain and a Peptide chain of Interest! </h2>")#print lig_content

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
        ProtLigfolder= '/opt/lampp/htdocs/MSALigMap/tmp/' + 'ProtPep'

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
        prot_lig_dict=collections.defaultdict(dict)
        count=1

        
        pdbid_chain_dict={}
        for i in l:
            each= i#.split(':')
            pdbcode= each.upper()
            ##chain=each[1] # i am removing the +"="+"[]"
            pdbid_list.append(pdbcode)
            #pdbid_chain_dict[pdbcode]=chain #(not used i later part of the code)
        


        #linking the PDB url address to the selected PDB ids
        for x in l:
            pdbcodename= x.split(':')[0]
            pdbid=url+pdbcodename+".pdb"
            url_list.append(pdbid)


        #Listing the PDB ids and their bound ligands
        for link,id  in zip(url_list,pdbid_list): #list(pdbid_chain_dict.keys())):
                #prot_lig_dict[id]={}
                dict_primarykey = id #+ '_' + Chain
                count=count+1
                try:
                    f=urllib.request.urlopen(link)
                    f=f.readlines()
                    #print (f)
                    for li in f:
                        if li.startswith(b'SEQRES'):
                            
                            li= li.decode('UTF-8').split()
                            Chain=li[2]
                            prot_lig_dict[dict_primarykey][Chain]={}
                            Seq_length= int(li[3])
                            
                            if Seq_length > 100:
                                prot_lig_dict[dict_primarykey][Chain]['Length']= Seq_length
                                prot_lig_dict[dict_primarykey][Chain]['Type']= 'Protein'
          
                            if Seq_length < 100:
                                prot_lig_dict[dict_primarykey][Chain]['Length']= Seq_length
                                prot_lig_dict[dict_primarykey][Chain]['Type']= 'Peptide'
                        
                except:
                    print("<p>The entered PDB code is not identitifed in PDB database. Try again!</p>")

        



        print("<h3>Sequences that are shorter than 100 amino acid in length is named as peptide</h3>")

        # Form for selecting the Peptide chains of interest for the USERR
        

        print ("<form enctype='multipart/form-data' action='LigPagePeptide.py' method = 'post' target = '_blank'>")
        
        print ("<table style=width:50%>")
        print ("<tr>")
        print ("<th>PDB ID</th>")
        print ("<th >Chain</th>")
        print ("<th >Length</th>")
        print ("<th >Type</th>")
        print ("<th >Select</th>")

        #length_of_rowspan= len(prot_lig_dict) 


        for k,dk in prot_lig_dict.items():
            
            length_of_rowspan = len(prot_lig_dict[k])
            k= k.split('_')[0]
            

            
            print ("<tr>")
            print ("<th rowspan='%d'>"% length_of_rowspan,k, "</th>")
            #print "</tr>"
            for chainKey, lengthValue in dk.items():
                
                print ("<td align=center>", chainKey,"</td>")
                print ("<td align=center>", dk[chainKey]['Length'],"</td>")
                print ("<td align=center>", dk[chainKey]['Type'],"</td>")
                print ("<td align=center>")
                if dk[chainKey]['Type']== 'Protein':
                    values_add = k + ':' + chainKey 
                    print  ("<input type='radio' name='ProteinSelection' value='%s'/>" % (values_add))
                
                if dk[chainKey]['Type']== 'Peptide':
                    values_add = k + ':' + chainKey 
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
        print ("<p><b>Please upload your sequence file and PDB code. Try again!</b></p>")

else:
    print ("<h3> Please enter a PDB code and try again!</h3>")
    sys.exit()



print ("</body>")
print ("</html>")