#!/usr/bin/env python3

# Import modules for CGI handling 
import cgi, cgitb 
cgitb.enable()
import os
import argparse
import sys
from Bio import SeqIO 
from Bio.Align.Applications import MuscleCommandline
from Bio.Align.Applications import ClustalOmegaCommandline
#from StringIO import StringIO
from Bio import AlignIO
from Bio.PDB import *
from Bio.PDB.Polypeptide import PPBuilder
import subprocess
import urllib
from collections import OrderedDict
import requests
from bs4 import BeautifulSoup
import itertools
import shutil
from Bio import motifs 
from Bio.Seq import Seq 
sys.path.append('/Users/xhasam/anaconda2/lib')

# Create instance of FieldStorage
form = cgi.FieldStorage()

  
# Get data from field
if form.getvalue('textcontent'):
   text_content = form.getvalue('textcontent')
else:
   text_content = "Not entered"
   
#if form.getvalue('textcontent'):
#   lig_content = form.getvalue('ligatom')
#else:
#   lig_content = "Not entered"

print ("Content-type:text/html\r\n\r\n")
print ("<html>")
print ("<head>")

#HTML style details for web page
print( "<style>")
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
print ("<img src='Pic.gif' align='middle' width='80%' height='200'")
print ("</div>")
print ("<body>")
print ("<ul>")
print ("<li><a href='HomePage.py'>Home</a></li>")
print ("<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>")
print ("<li><a href='ProteinDNA.py'>Protein-DNA</a></li>")
print ("<li><a href='ProteinPeptide.py'>Protein-Peptide</a></li>")
print ("<li><a href=''>Contact</a></li>")
print ("</ul>")
print ("<div align='center'>")
print ("<h2> Below is shown the mapped binding sites for the sequences </h2>")


#capturing the entered pdb ids into list
for i in text_content:
        text_content1=text_content.replace(' ', '')
        l=text_content1.split(',')
# print l

foldernamer= ''.join(random.choices(string.ascii_letters, k=4))
os.mkdir('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/')

fileitem = form['filename']
fileattached = fileitem.value
with open('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/sequenceInputfile.fasta', 'w') as fout:
    fout.write(fileattached)

out_file='/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/seqaligned1.fasta'
muscle='/usr/local/bin/muscle'
clustalo='/usr/local/bin/clustalo'
pdbl=PDBList()
ppb=PPBuilder()
pdb_id=[] 
pdburl="https://files.rcsb.org/download/"
PDB="/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/PDB/"
pdb_seq_dict=OrderedDict()
pdb_seq_dict_numbering=OrderedDict()
pdb_seqfin=OrderedDict()
non_pdb_id=[]
non_pdb_seq_dict=OrderedDict()
non_pdb_seq_dictfin=OrderedDict()
pdbid_lig=OrderedDict()
check_aa=[]
NONHcheck_aa=[]
updated_indexed_bindingsite=OrderedDict()
updated_indexed_NONHbindingsite=OrderedDict()

Alignment_adjusted_indexed_Hbindingsite=OrderedDict()
Alignment_adjusted_indexed_NONHbindingsite=OrderedDict()

##dssp
dssp_urlpart1='https://mrs.cmbi.umcn.nl/search?db=all&q='
dssp_urlpart2='&count=3'
dssp_final_url1='https://mrs.cmbi.umcn.nl/entry?db=dssp&nr='
dssp_final_url2='&rq='

##dssp


aa_dict={'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

pdbsum_URL="http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode="
pdbsum_URL2="&template=links.html"
ebiurl="www.ebi.ac.uk"
pdbsum_dict=OrderedDict()
DSSP_SS=OrderedDict()
Prot_AA=OrderedDict()
aa_text_dict=OrderedDict()
PDB_code_ligand_dict ={}
for items_names in l:

    items_names_split= items_names.split('_')
    ligand_name= items_names_split[1]
    PDB_id = items_names_split[0].split(':')[0]
    PDB_code_ligand_dict[PDB_id]={}
    Chain_id = items_names_split[0].split(':')[1] 
    PDB_code_ligand_dict[PDB_id]['Chain']=Chain_id
    PDB_code_ligand_dict[PDB_id]['Ligand']=ligand_name
    pdbid_lig.setdefault(items_names_split[0],[]).append(ligand_name)


# fileitem = form['filename']

# # with open( 'Applications/XAMPP/xamppfiles/htdocs/MSALigMap/sequence1.fasta', 'w')  as f:
# if fileitem.filename:
    
#     # fileitem.file.write()
#     # print (fileitem.file.read(4096))
#     # strip leading path from file name to avoid
#     # directory traversal attacks
#     fn = os.path.basename(fileitem.filename)
#     print ('filename',fn)
#     print(fileitem.file.read())
#     ss = open('/Applications/XAMPP/xamppfiles/htdocs/MSALigMap/' + fn, 'wb')
#     ss.write('SMA')
#     for lines in fileitem.file.read(4096):
#         print (lines)
#     message = 'The file "' + fn + '" was uploaded successfully'
 
# else:
#     message = 'No file was uploaded'



def merge_pdb_nonpdb(dict1,dict2):
    combine=dict1.copy()
    combine.update(combine)
    return combine

with open('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/sequenceInputfile.fasta', 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        ids_pdb=record.id.split(':')[:1]
        if len(''.join(ids_pdb)) == 4:
            chain_name=record.id.split(':')[1:]
            ext=chain_name[0]
            filesset=pdbl.download_pdb_files(ids_pdb, obsolete=False, pdir=PDB, file_format="pdb", overwrite=False)
        if len(''.join(ids_pdb)) > 4:
            non_pdb_id.append(record.id)
            swis_seq=str(record.seq)
            non_pdb_seq_dict.setdefault(record.id,[]).append(record.seq)
            non_pdb_seq_dictfin.setdefault(record.id,[]).append(swis_seq)
            
        
        pdb_id.append(record.id)
    for data in os.listdir(PDB):
        paths='/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/PDB/'+ data
        idschange=paths.split("/")[1][3:]
        pdbid_dict=os.path.splitext(idschange)[0]+':'+ext
                
        with open(paths, 'r') as f:
            for line in f:
                lines=line.split()
                if lines[0]=="ATOM" and lines[2]=="CA" and lines[4]==ext  :
                    if len(lines[3])==3:
                        found= (aa_dict[key] for key in aa_dict.keys() if key == lines[3])
                        for item in found:
                            pdb_seq_dict.setdefault(pdbid_dict, []).append(item)
                            pdb_seq_dict_numbering.setdefault(pdbid_dict, []).append(lines[5])
                    elif lines[3][0]==ext:
                        found1= (aa_dict[key] for key in aa_dict.keys() if key == lines[3][1:])
                        for item1 in found1:
                            pdb_seq_dict.setdefault(pdbid_dict, []).append(item1)
                            pdb_seq_dict_numbering.setdefault(pdbid_dict, []).append(lines[5])
print("<br/>")
print("<br/>")



for keys,vals in pdb_seq_dict.items():
    aas_string=''.join(map(str,vals))
    pdb_seqfin.setdefault(keys,[]).append(aas_string)
 

combine=dict(list(non_pdb_seq_dictfin.items()) + list(pdb_seqfin.items()))
with open('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/trimmedfasta.fasta', 'w') as files:
    for seqids, seqn in combine.items():
        files.write( ">" + seqids)
        files.write("\n")
        files.write (seqn[0])
        files.write("\n") 

Omega_cline = subprocess.call([clustalo, '--infile', '/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/trimmedfasta.fasta', '--outfile', out_file])

with open('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/file_tab.txt', 'w') as files_tab:
    for record in SeqIO.parse(out_file, 'fasta'):
        tabform='{}\t{}'.format(record.description, record.seq)
        files_tab.write(tabform)
        files_tab.write('\n')

#Preparing PdbSum Url with selected PDB ids
for ids,lig in pdbid_lig.items():
    ids_splitslig=str(ids).split(':')[0]
    pdbsumurl=pdbsum_URL+ids_splitslig+pdbsum_URL2
    pdbsum_dict.setdefault('%s'%ids,[]).append(pdbsumurl)
#Extracting the Href links from PDBSum home page for the selected PDB ids and Ligands using BeautifulSoup

mydictcheck=OrderedDict()
pdbsum_peptide_url1= 'http://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode='
pdbsum_peptide_url2= '&template=interfaces.html&c=999'

for pdbpep,lid in pdbid_lig.items():
    pdbid_pep = pdbpep.split(':')[0]
    chain_first = pdbpep.split(':')[1]
    chain_second = lid[0]
    pdbsum_peptide_concaat_url = pdbsum_peptide_url1 + pdbid_pep + pdbsum_peptide_url2
    page = requests.get(pdbsum_peptide_concaat_url)
    
    chainurl1 = '/thornton-srv/databases/pdbsum/templates/gif/chain' + chain_first + '.jpg'
    chainurl2 = '/thornton-srv/databases/pdbsum/templates/gif/chain' + chain_second + '.jpg'
    temp_ebiurl = []
    soup = BeautifulSoup(page.content, 'html.parser')
    for h in soup.find_all('a', class_='menuClass' ):
        if chainurl1 in str(h) and chainurl2 in str(h):
            if not h['href'] in temp_ebiurl:
                temp_ebiurl.append(h['href'])
            for temp_ebiurl1 in temp_ebiurl:
                linkfin="http://" + ebiurl + temp_ebiurl1
                temp_ebiurl2 = []
                page2 = requests.get(linkfin)
                soup2 = BeautifulSoup(page2.content, 'html.parser')
                for h2 in soup2.find_all('a'):
                    search_1= '/thornton-srv/databases/cgi-bin/pdbsum/GetIface.pl?pdb='
                    search_2 = 'List of interactions'
                    if search_1 in str(h2) and search_2 in str(h2):
                        
                        if not h2['href'] in temp_ebiurl2:
                            temp_ebiurl2.append(h2['href'])
                            pep_pep_interactionlink ="http://www.ebi.ac.uk" + h2['href']
                            mydictcheck[pdbpep]=pep_pep_interactionlink

#selecting common ligand atoms that are hydrogen bonded in selected PDB structures
H_printing = False
H_atoms_commoncomp={}
for pdbsumids, pdbsumlink in mydictcheck.items():
    H_links_sel1=str(pdbsumlink).replace("&amp;", "&")
    pdbsumids_chain=str(pdbsumids).split(":")[1]

    weblink=requests.get(H_links_sel1)
    for H_atomlines in weblink.iter_lines():
        H_atomlines1=H_atomlines.strip()
        
        if H_atomlines1.startswith(b'Hydrogen bonds'):
            H_printing = True
        elif H_atomlines1.startswith(b'Non-bonded contacts'):
            H_printing = False
        if H_printing:
            if H_atomlines1.startswith((b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9')):
                H_atomlines2=H_atomlines1.split()
                H_atomlines2 = [w1.decode('utf8').replace('b', '') for w1 in H_atomlines2]
                H_atm_sel=H_atomlines2[3]
                if pdbsumids_chain==H_atomlines2[5]:
                    concat_aa_no=H_atm_sel+"_"+H_atomlines2[4]
                    if concat_aa_no not in check_aa:
                        check_aa.append(concat_aa_no)
                        found_aa= (aa_dict[key] for key in aa_dict.keys() if key == H_atomlines2[3])
                        for items_aa in found_aa:
                            H_atm_sel=items_aa+"_"+H_atomlines2[4]
                            H_atoms_commoncomp.setdefault(pdbsumids,[]).append(H_atm_sel)
    check_aa=[]


NONH_printing = False
NONHatoms_commoncomp={}

for NONHpdbsumids,NONHpdbsumlink in mydictcheck.items():
    NONHlinks_sel1=str(NONHpdbsumlink).replace("&amp;", "&")
    NONHpdbsumids_chain=str(NONHpdbsumids).split(":")[1]
    weblink1=requests.get(NONHlinks_sel1)
    for NONHatomlines in weblink1.iter_lines():
        NONHatomlines1=NONHatomlines.strip()
        if NONHatomlines1.startswith(b"Non-bonded contacts"):
            NONH_printing = True
            
        elif NONHatomlines1.startswith(b"Hydrogen bonds"):
            NONH_printing = False
        if NONH_printing:                
            if NONHatomlines1.startswith((b'0', b'1', b'2', b'3', b'4', b'5', b'6', b'7', b'8', b'9')):
                
                NONHatomlines2=NONHatomlines1.split()
                NONHatomlines2 = [w.decode('utf8').replace('b', '') for w in NONHatomlines2]
                NONHatm_sel=NONHatomlines2[3]
                if NONHatomlines2[5]==NONHpdbsumids_chain:
                   NONHconcat_aa_no=NONHatm_sel+"_"+ NONHatomlines2[4]
                   if NONHconcat_aa_no not in NONHcheck_aa:
                       NONHcheck_aa.append(NONHconcat_aa_no)
                       found_aa= (aa_dict[key] for key in aa_dict.keys() if key == NONHatomlines2[3])
                       for items_aa in found_aa:
                           NONH_atm_sel=items_aa+"_"+NONHatomlines2[4]
                           NONHatoms_commoncomp.setdefault(NONHpdbsumids,[]).append(NONH_atm_sel)
    NONHcheck_aa=[]                    

for ids,aa_numb1 in H_atoms_commoncomp.items():
    aa_numb= [w3.replace('u', '') for w3 in aa_numb1]
    for ids1, numb in pdb_seq_dict_numbering.items():
        if ids==ids1.upper():
            for j in aa_numb:
                numb_split=j.split('_')
                updated_indexed_bindingsite.setdefault(ids,[]).append(numb.index(numb_split[1]))

for ids_1,aa_numb1 in NONHatoms_commoncomp.items():
    for ids2, numb2 in pdb_seq_dict_numbering.items():
        if ids_1==ids2.upper():
            for j1 in aa_numb1:
                numb_split1=j1.split('_')
                updated_indexed_NONHbindingsite.setdefault(ids_1,[]).append(numb2.index(numb_split1[1]))


with open("/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/file_tab.txt",'r') as files:
    for line in files:
        
        line1=line.split()[1:]
        SS_list=[]
        pdbids_chain=line.split()[:1]
        
        pdbids_chain_split=''.join(pdbids_chain).split(':')
        if len(pdbids_chain_split) > 1:
            
            check_ids=''.join(pdbids_chain).split(':')[0]#[1:]
            
            pdb_chain=''.join(pdbids_chain).split(':')[1]
            dssp_urlpart12_concat= dssp_urlpart1+check_ids+dssp_urlpart2
            
            dssp_allDBpage = requests.get(dssp_urlpart12_concat)
            dsspsource = BeautifulSoup(dssp_allDBpage.content, 'html.parser')
            for MRSallDB_link in dsspsource.find_all('a'):
                for dssplink in MRSallDB_link:
                    
                    if dssplink==(check_ids.upper()):
                        
                        dssp_href= str(MRSallDB_link).split(';')
                        
                        if dssp_href[0]=='<a href="entry?db=dssp&amp':
                            nr_number=dssp_href[1].split('=')[1].split('&')[0]
                            dssp_finalurl_concatenate= dssp_final_url1 + nr_number+ dssp_final_url2+check_ids.upper()
                            dsspwebpagelink=requests.get(dssp_finalurl_concatenate)
                            for dssppagelines in dsspwebpagelink.iter_lines():
                                lines_in_dssp = dssppagelines
                                
                                if not lines_in_dssp.startswith('<') and len(lines_in_dssp)>3:
                                    linesdssp_part1=lines_in_dssp.split()
                                    if len(linesdssp_part1)>5 and linesdssp_part1[2]==pdb_chain:
                                        if str(linesdssp_part1[4]).startswith(('H','B','E','G','I','T', 'S')):
                                            aa_SS=linesdssp_part1[4]
                                            SS_list.append(aa_SS)
                                        else:
                                            aa_NoSS='-'
                                            SS_list.append(aa_NoSS)
                            #Secondary structure adjustment according to its sequence alignment
                            if len(pdbids_chain_split) > 1:
                                for seq_residues in line1:
                                    found2=[res_letter for res_letter in list(seq_residues)]
                                sscount=0
                                sscount1=0
                                num=0
                                SS_list_alignment_adjusted=[]
                                for  char,ssi in zip(found2,itertools.cycle(SS_list)): 
                                    
                                    if char!='-'  :
                                        SS_list_alignment_adjusted.append(SS_list[num])
                                        num+=1
                                        sscount+=1
                                    else:
                                        sscount1+=1
                                        SS_list_alignment_adjusted.append('-')
                                for t in range(0,len(''.join(SS_list_alignment_adjusted)),60):
                                    SS_temp=[]
                                    SS_temp.append(''.join(SS_list_alignment_adjusted)[t:t+60])
                                    DSSP_SS.setdefault(check_ids, []).append(''.join(SS_temp))
            PDBseqID_print=">"+''.join(check_ids)+'\t'+'\t'
        else:
            check_ids=''.join(pdbids_chain).split(':')[0]
            NonPDBseqID_print=''.join(check_ids)+'\t'+'\t'
        
        aa_text=list(str(line1)[2:])
        aa_text=aa_text[:-2]
        aa_text_dict[check_ids]=aa_text
        for aa_t in range(0,len(''.join(aa_text)),60):
            AA_temp=[]
            AA_temp.append(''.join(aa_text)[aa_t:aa_t+60])
            if check_ids.startswith('>'):
                check_ids=''.join(check_ids.split('>'))
                Prot_AA.setdefault(check_ids, []).append(''.join(AA_temp))
            else:
                Prot_AA.setdefault(check_ids, []).append(''.join(AA_temp))     

####updating index number for binding site based on alignment for the given sequences

for  Hkey1,Hvalues1 in updated_indexed_bindingsite.items():  
    H_count =0
    H_count1=0
    Hkey1=str(Hkey1).split(':')[0]
    if Hkey1.lower() in aa_text_dict.keys():
        for H_X, H_i in enumerate(aa_text_dict[str(Hkey1.lower()).split(':')[0]]):
                if H_i=='-':
                    H_count+=1
                else:
                     if H_count1 in Hvalues1:
                         Alignment_adjusted_indexed_Hbindingsite.setdefault(Hkey1, []).append(H_X)
                     H_count1+=1



# print (Alignment_adjusted_indexed_Hbindingsite)
for  NONHkey1,NONHvalues1 in updated_indexed_NONHbindingsite.items():  
    NONH_count =0
    NONH_count1=0
    NONHkey1=str(NONHkey1).split(':')[0]
    if NONHkey1.lower() in aa_text_dict.keys():
        for NONH_X, NONH_i in enumerate(aa_text_dict[NONHkey1.lower()]):
                if NONH_i=='-':
                    NONH_count+=1
                else:
                     if NONH_count1 in NONHvalues1:
                         Alignment_adjusted_indexed_NONHbindingsite.setdefault(NONHkey1, []).append(NONH_X)
                     NONH_count1+=1
prot_Idsname=list(Prot_AA.keys())
digits1=len(max(prot_Idsname))
f1 = '{0:>%d}: ' % (digits1)


print("<b>Color coded secondary structure elements </b> ")
print('<br/>')
structurecode="<span style='background-color:#8B008B'><font color='white'>B = beta-bridge residue</font></span>, <span style='background-color:#FFFF00'>E = extended strand (in beta ladder)</span>, <span style='background-color:#CD5C5C'><font color='white'>G = 3/10-helix</font></span>, <span style='background-color:#FF0000'><font color='white'>H= alpha-helix</font></span>,    <span style='background-color:#FA8072'>I = Pi helix</span>, <span style='background-color: #00FF00'> S = bend </span>, <span style='background-color: #008000'><font color='white'>T = H-bonded turn</font></span>"
print(structurecode)
print('<br/>')
print('<br/>')
print("<table style= 'border: 0'>")# align='center'>")
##holding the keys of  Prot_AA dictionary
SeqID_ref= list(sorted(Prot_AA.keys()))[0]

for itemlength in range(len(Prot_AA[SeqID_ref])):
    
    for dictkeys in sorted(Prot_AA.keys()):
        
        gapcount=0
        resindex=0
        
        if dictkeys in DSSP_SS:
            
            print("<tr style= 'border: 0'>")
            print("<td style= 'border: 0;width=5%'>")
            DSS_ID_print1=dictkeys+'_Sec_Str'+'\t'
            print(DSS_ID_print1)
            print("</td>")    
        if itemlength==0:
            
            count=0
            resindex=0
        elif  itemlength>0:
            count=itemlength*60
            resindex=itemlength*60
        
        if dictkeys in DSSP_SS.keys():
            
            ss=DSSP_SS[dictkeys][itemlength]
            print("<td style= 'border: 0;width=80%'>")
            for sec_str in ss:
                
                
                if sec_str=='H':
                    helix= "<span style='background-color:#FF0000;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(helix)
                elif sec_str=='E':
                    strand="<span style='background-color: #FFFF00;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(strand)
                elif sec_str=='S':
                    bend="<span style='background-color: #00FF00;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(bend)  
                elif sec_str=='T':
                    turn="<span style='background-color: #008000;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(turn)
                elif sec_str=='G':
                    helix_310="<span style='background-color: #CD5C5C;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(helix_310)
                elif sec_str=='I':
                    Pi_helix="<span style='background-color: #FA8072;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(Pi_helix)
                elif sec_str=='B':
                    beta_bridge="<span style='background-color: #8B008B;font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(beta_bridge)
                else:
                    other_SS="<span style='font-family:monospace;font-size:18px'>%s</span>"%sec_str
                    print(other_SS)
        print("</td>")
        print("</tr>")
        print("<tr style= 'border: 0'>")
        print("<td style= 'border: 0';width=5%'>")
        dictkeys2=dictkeys+'\t'
        dictkeys1="<text-align='left'>%s" %dictkeys2
        print(dictkeys1)
        print("</td>")
        ##print secondary structure and sequences
        print("<td style= 'border: 0'; width=80%'>")
        for i in Prot_AA[dictkeys][itemlength]:
            
            
            if i=='-':
                res_size= "<span style='font-family:monospace;font-size:18px'>%s</span>" %i
                print(res_size)
                count+=1
                gapcount+=1
            elif dictkeys.upper() in Alignment_adjusted_indexed_Hbindingsite.keys() or dictkeys.upper() in Alignment_adjusted_indexed_NONHbindingsite.keys():
                if count in Alignment_adjusted_indexed_Hbindingsite[dictkeys.upper()]: #and not Alignment_adjusted_indexed_NONHbindingsite[dictkeys]:
                    if i.startswith(('A','I','L','M','V')):
                        bolded="<b><font color='#FF1493'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    if i.startswith(('F','W', 'Y')):
                        bolded="<b><font color='#FF8C00'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded) 
                    if i.startswith(('K','R', 'H')):
                        bolded="<b><font color='red'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    if i.startswith(('E','D')):
                        bolded="<b><font color='#006400'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    if i.startswith(('N','Q','S','T')) :
                        bolded="<b><font color='blue'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    if i.startswith(('G','P')):
                        bolded="<b><font color='#800080'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    if i.startswith(('C')):
                        bolded="<b><font color='yellow'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</b>'
                        print(bolded)
                    #resindex+=1
                    count+=1
                elif  count in Alignment_adjusted_indexed_NONHbindingsite[dictkeys.upper()]: #and not Alignment_adjusted_indexed_Hbindingsite[dictkeys]:
                    if i.startswith(('A','I','L','M','V')):
                        underlined="<u><font color='#FF1493'><span style='font color: #FF1493;font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    if i.startswith(('F','W', 'Y')):
                        underlined="<u><font color='#FF8C00'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined) 
                    if i.startswith(('K','R', 'H')):
                        underlined="<u><font color='red'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    if i.startswith(('E','D')):
                        underlined="<u><font color='#006400'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    if i.startswith(('N','Q','S','T')):
                        underlined="<u><font color='blue'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    if i.startswith(('G','P')):
                        underlined="<u><font color='#800080'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    if i.startswith(('C')):
                        underlined="<u><font color='yellow'><span style='font-family:monospace;font-size:18px'>%s</span></font>" %i+'</u>'
                        print(underlined)
                    #resindex+=1
                    count+=1

                    
                else:
                    res_sizePDB= "<span style='font-family:monospace;font-size:18px'>%s</span>" %i
                    print(res_sizePDB)
                    #resindex+=1
                    count+=1
            else:
                res_sizeNonPDB= "<span style='font-family:monospace;font-size:18px'>%s</span>" %i
                print(res_sizeNonPDB)
        print("</td>")
    #     print("</pre>")
        print("</tr>")



    print("<tr style= 'border: 0; heigth: 60px '>")
    print("<td style= 'border: 0; height: 60px'>")
    print("</td>")
    print("</tr>")
        
print("</table>")


## making of Weblog and binding residues alignment
# merginf og hydrogen and nonhydrogen bondong residues        
mergeH_H_NHddict= {key: value + Alignment_adjusted_indexed_NONHbindingsite[key] for key, value in Alignment_adjusted_indexed_Hbindingsite.items()}

#making a unique list for each key in dictionary by removing duplicate items in the list
unique_mergeddict_H_NH= {k:sorted(set(j),key=j.index) for k,j in mergeH_H_NHddict.items()}

#making a list of binding residue indexes (according to aligned positions) from unique_mergeddict
H_NHnew_list=[]
for H_NHkey in unique_mergeddict_H_NH.keys():
    for index_items in unique_mergeddict_H_NH[H_NHkey]:
        if index_items not in H_NHnew_list:
            H_NHnew_list.append(index_items)

#making a dictionary for creating a weblogo
weblogo_align=OrderedDict()
weblogo_gap_to_X=OrderedDict()

for pdbid in aa_text_dict.keys():
    for H_NHnum in sorted(H_NHnew_list):
        # print ("H_NHnum", H_NHnum)
        H_NHresidues=aa_text_dict[pdbid][H_NHnum]
        weblogo_align.setdefault(pdbid, []).append(H_NHresidues)# dictionary where gaps are not replaced with X
        if aa_text_dict[pdbid][H_NHnum]!='-':
            residuesnew=aa_text_dict[pdbid][H_NHnum]
            weblogo_gap_to_X.setdefault(pdbid, []).append(residuesnew)
        else:
            residuesX='X'
            weblogo_gap_to_X.setdefault(pdbid, []).append(residuesX)# dictionary where gaps are replaced with X
        

#### For hydrongen bonded interactions only
H_new_list=[]
for key in Alignment_adjusted_indexed_Hbindingsite.keys():
    for H_index_items in Alignment_adjusted_indexed_Hbindingsite[key]:
        if H_index_items not in H_new_list:
            H_new_list.append(H_index_items)
    
H_weblogo_align=OrderedDict()
H_weblogo_gap_to_X=OrderedDict()

for H_pdbid in aa_text_dict.keys():
    for H_num in sorted(H_new_list):
        H_residues=aa_text_dict[H_pdbid][H_num]
        H_weblogo_align.setdefault(H_pdbid, []).append(H_residues)
        if aa_text_dict[H_pdbid][H_num]!='-':
            H_residuesnew=aa_text_dict[H_pdbid][H_num]
            H_weblogo_gap_to_X.setdefault(H_pdbid, []).append(H_residuesnew)
        else:
            H_residuesX='X'
            H_weblogo_gap_to_X.setdefault(H_pdbid, []).append(H_residuesX)


print("<div align='center'><b>Alignment of hydrogen bonded interacting residues</b></div>")
print("<br/>")


H_identical_index=[]
for H_ind, H_elems in enumerate(zip(*H_weblogo_align.values())):
        
        if H_elems[1:] == H_elems[:-1]:
            H_identical_index.append(H_ind)


print("<table style= 'border: 0'>")
for H_dictkeys in sorted(H_weblogo_align.keys()):
    print("<tr style= 'border: 0'>")
    print("<td style= 'border: 0'>")
    H_ident_count=0
    if str(H_dictkeys).startswith('>'):
        H_dictkeys_split=str(H_dictkeys).split('>')

        print(str(H_dictkeys_split[1]))
        print("</td>")
        
        print("<td>")
        for H_res in H_weblogo_align[H_dictkeys]:
            if H_ident_count in H_identical_index:
                H_underlined="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white'>%s</SPAN>" %H_res
                print(H_underlined)
                H_ident_count+=1
            else:
                H_no_conse="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black'>%s</SPAN>" %H_res+'</b>'
                print(H_no_conse)
                H_ident_count+=1
    else:
        print(H_dictkeys)
        print("</td>")
        print("<td style= 'border: 0'>")
        for H_res in H_weblogo_align[H_dictkeys]:
            if H_ident_count in H_identical_index:
                H_underlined1="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white'>%s</SPAN>" %H_res
                print(H_underlined1)
                H_ident_count+=1
            else:
                H_no_conse1="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black'>%s</SPAN>" %H_res+'</b>'
                print(H_no_conse1)
                H_ident_count+=1
        
    print("</td>")
    print("</tr>")
print("</table>")


# #making H_weblogo
# H_protmotif=[]
# H_bonding='Hbonding'+'.svg'
# print ("H_weblogo_gap_to_X",H_weblogo_gap_to_X)
# for H_webid in H_weblogo_gap_to_X.keys():
#     H_var=''.join(H_weblogo_gap_to_X[H_webid])
#     H_protmotif.append(Seq(H_var, generic_protein))
#     H_seq = motifs.create(H_protmotif) 
#     print (H_seq)
# H_seq.weblogo(H_bonding,format='SVG',xaxis_label= 'Residues', show_errorbars= False, color_scheme= 'color_chemistry')
# Hsrc="<center><embed src='%s#page=1&view=FitH ' /></center>"%H_bonding
# print(Hsrc)
# print("<br/>")
# print("<br/>")


#### For Non_hydrongen bonded interactions only            
NH_new_list=[]
for NHkey in Alignment_adjusted_indexed_NONHbindingsite.keys():
    for NHindex_items in Alignment_adjusted_indexed_NONHbindingsite[NHkey]:
        if NHindex_items not in NH_new_list:
            NH_new_list.append(NHindex_items)
    
NH_weblogo_align=OrderedDict()
NH_weblogo_gap_to_X=OrderedDict()

for NH_pdbid in aa_text_dict.keys():
    for NH_num in sorted(NH_new_list):
        NH_residues=aa_text_dict[NH_pdbid][NH_num]
        NH_weblogo_align.setdefault(NH_pdbid, []).append(NH_residues)
        if aa_text_dict[NH_pdbid][NH_num]!='-':
            NH_residuesnew=aa_text_dict[NH_pdbid][NH_num]
            NH_weblogo_gap_to_X.setdefault(NH_pdbid, []).append(NH_residuesnew)
        else:
            NH_residuesX='X'
            NH_weblogo_gap_to_X.setdefault(NH_pdbid, []).append(NH_residuesX)


#print("Alignment of Non-Hydrogen Interaction Residues")
print("<br/>")
print("<br/>")
print("<div align='center'><b>Alignment  of non-bonded interacting residues</b></div>")
print("<br/>")

NH_identical_index=[]
for NH_ind, NH_elems in enumerate(zip(*NH_weblogo_align.values())):
        
        if NH_elems[1:] == NH_elems[:-1]:
            NH_identical_index.append(NH_ind)


#NHSeqID_ref= list(sorted(NH_weblogo_align.keys()))[0]
print("<table style= 'border: 0'>")
for NH_dictkeys in sorted(NH_weblogo_align.keys()):
    print("<tr style= 'border: 0'>")
    print("<td style= 'border: 0'>")
    NH_ident_count=0
    if str(NH_dictkeys).startswith('>'):
        NH_dictkeys_split=str(NH_dictkeys).split('>')

        print(str(NH_dictkeys_split[1]))
        print("</td>")
        
        print("<td style= 'border: 0'>")
        for NH_res in NH_weblogo_align[NH_dictkeys]:
            if NH_ident_count in NH_identical_index:
                NH_underlined="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white'>%s</SPAN>" %NH_res
                print(NH_underlined)
                NH_ident_count+=1
            else:
                NH_no_conse="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black'>%s</SPAN>" %NH_res+'</b>'
                print(NH_no_conse)
                NH_ident_count+=1
    else:
        print(NH_dictkeys)
        print("</td>")
        print("<td style= 'border: 0'>")
        for NH_res in NH_weblogo_align[NH_dictkeys]:
            if NH_ident_count in NH_identical_index:
                NH_underlined1="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white'>%s</SPAN>" %NH_res
                print(NH_underlined1)
                NH_ident_count+=1
            else:
                NH_no_conse1="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black'>%s</SPAN>" %NH_res+'</b>'
                print(NH_no_conse1)
                NH_ident_count+=1
        
    print("</td>")
    print("</tr>")
print("</table>")

# #making H_NH_weblogo
# NH_protmotif=[]
# NH_bonding='Nhbonding'+'.svg'
# print ("NH_weblogo_gap_to_X", NH_weblogo_gap_to_X)
# for NH_webid in NH_weblogo_gap_to_X.keys():
#     NH_var=''.join(NH_weblogo_gap_to_X[NH_webid])
#     NH_protmotif.append(Seq(NH_var, generic_protein))
#     NH_seq = motifs.create(NH_protmotif) 
# NH_seq.weblogo(NH_bonding,format='SVG',xaxis_label= 'Residues', show_errorbars= False, color_scheme= 'color_chemistry')
# NHsrc="<center><embed src='%s#page=1&view=FitH ' /></center>"%NH_bonding
# print(NHsrc)
# print("<br/>")
# print("<br/>")


#print("Alignment of Both Hydrogen and Non-hydrogen interaction Residues")
print("<br/>")
print("<br/>")
print("<div align='center'><b>Alignment  of hydrogen and non-bonded interacting residues</b></div>")
print("<br/>")

H_NH_identical_index=[]
for H_NH_ind, H_NH_elems in enumerate(zip(*weblogo_align.values())):
        
        if H_NH_elems[1:] == H_NH_elems[:-1]:
            H_NH_identical_index.append(H_NH_ind)


#H_NHSeqID_ref= list(sorted(weblogo_align.keys()))[0]
print("<table style= 'border: 0'>")
for H_NH_dictkeys in sorted(weblogo_align.keys()):
    print("<tr style= 'border: 0'>")
    print("<td style= 'border: 0'>")
    H_NH_ident_count=0
    if str(H_NH_dictkeys).startswith('>'):
        H_NH_dictkeys_split=str(H_NH_dictkeys).split('>')

        print(str(H_NH_dictkeys_split[1]))
        print("</td>")
        
        print("<td style= 'border: 0'>")
        for H_NH_res in weblogo_align[H_NH_dictkeys]:
            if H_NH_ident_count in H_NH_identical_index:
                H_NH_underlined="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white;font-family:monospace;font-size:18px'>%s</SPAN>" %H_NH_res
                print(H_NH_underlined)
                H_NH_ident_count+=1
            else:
                H_NH_no_conse="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black;font-family:monospace;font-size:18px'>%s</SPAN>" %H_NH_res+'</b>'
                print(H_NH_no_conse)
                H_NH_ident_count+=1
    else:
        print(H_NH_dictkeys)
        print("</td>")
        print("<td style= 'border: 0'>")
        for H_NH_res in weblogo_align[H_NH_dictkeys]:
            if H_NH_ident_count in H_NH_identical_index:
                H_NH_underlined1="<SPAN STYLE='background-color:red;font-family:monospace;font-size:18px; font-weight:bold; color:white;font-family:monospace;font-size:18px'>%s</SPAN>" %H_NH_res
                print(H_NH_underlined1)
                H_NH_ident_count+=1
            else:
                H_NH_no_conse1="<b><SPAN STYLE='background-color:powderblue;font-family:monospace;font-size:18px; font-weight:bold; color:black;font-family:monospace;font-size:18px'>%s</SPAN>" %H_NH_res+'</b>'
                print(H_NH_no_conse1)
                H_NH_ident_count+=1
        
    print("</td>")
    print("</tr>")
print("</table>")

# #making H_NH_weblogo
# protmotif=[]
# H_NONH='HNONH'+'.svg'
# for webid in weblogo_gap_to_X.keys():
#     var=''.join(weblogo_gap_to_X[webid])
#     protmotif.append(Seq(var, generic_protein))
#     seq = motifs.create(protmotif) 
#     print (seq)
# seq.weblogo(H_NONH,format='SVG',xaxis_label= 'Residues', show_errorbars= False, color_scheme= 'color_chemistry')
# H_NONHsrc="<center><embed src='%s#page=1&view=FitH' /></center>"%H_NONH
# print(H_NONHsrc)
# print("<br/>")
# print("<br/>")

os.remove(out_file) 
os.remove('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/trimmedfasta.fasta') 
shutil.rmtree('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/obsolete') 
shutil.rmtree('/opt/lampp/htdocs/MSALigMap/tmp/ProtPep/foldernamer/PDB')

print ("</body>")
print ("</html>")
