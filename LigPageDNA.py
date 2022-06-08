#!/usr/bin/python

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
import re
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


print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"

#HTML style details for web page
print "<style>"
print "ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}"
print "li{float:left;}"
print "li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}"
print "li a:hover { background-color: #111111;}"
print "table, th, td { border: 2px solid black;}"
print ".footer { position: absolute; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }"
print "</style>"
#Style ends here

print "<title>LiBiSCo</title>"
print "</head>"
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
print "<div align='center'>"
print "<img src='Pic.gif' align='middle' width='80%' height='200'"
print "</div>"
print "<body>"
print "<ul>"
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>"
print "<li><a href='ProteinDNA.py'>Protein-DNA</a></li>"
print "<li><a href=''>Contact</a></li>"
print "</ul>"
print "<div align='center'>"
print "<h2> Below is shown the mapped binding sites for the sequences </h2>" 
#print lig_content
#capturing the entered pdb ids into list
for i in text_content:
        text_content1=text_content.replace(' ', '')
        l=text_content1.split(',')
# print l


out_file='seqaligned1.fasta'
muscle='/usr/local/bin/muscle'
clustalo='/usr/local/bin/clustalo'
pdbl=PDBList()
ppb=PPBuilder()
pdb_id=[] 
pdburl="https://files.rcsb.org/download/"
PDB="PDB/"
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

##NUCPLOT
NUcplot_URL1 = 'https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPS.pl?pdbcode='
NUcplot_URL2 = '&psfile=nucplot.ps'





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

def merge_pdb_nonpdb(dict1,dict2):
    combine=dict1.copy()
    combine.update(combine)
    return combine
PDB_ids_list= []
Exp_dict = {}
with open('WRKY.fasta', 'r') as f:
    for record in SeqIO.parse(f, "fasta"):
        # print (record.id)
        ids_pdb=record.id.split(':')[:1]
        # print (ids_pdb)
        chain_name=record.id.split(':')#[1:]
        if len(chain_name) > 1:
            pdb_code=chain_name[0]
            ext=chain_name[1]
            PDB_ids_list.append(pdb_code.lower()+':'+ext)
            # print (record.id)
            # print (pdb_code)
            # print (ext)

            filesset=pdbl.download_pdb_files(ids_pdb, obsolete=False, pdir=PDB, file_format="pdb", overwrite=False)
        if len(chain_name) == 1:
            non_pdb_id.append(record.id)
            swis_seq=str(record.seq)
            non_pdb_seq_dict.setdefault(record.id,[]).append(record.seq)
            non_pdb_seq_dictfin.setdefault(record.id,[]).append(swis_seq)
            
        
        pdb_id.append(record.id)
    for data in os.listdir(PDB):
        paths='PDB/'+ data
        idschange=paths.split("/")[1][3:]
        pdbid_dict=os.path.splitext(idschange)[0]+':'+ext
        with open(paths, 'r') as f:
            root,extension= os.path.splitext(paths)
            if extension == ".ent":
                for line in f:
                    lines=line.split()
                    if lines[0]=='EXPDTA' and lines[2]=='NMR':#.startswith('EXPDTA    SOLUTION NMR'):
                        NMR_printing = True
                        Exp_dict[pdbid_dict] = 'NMR'
                    elif lines[0]=='EXPDTA' and lines[1]=='X-RAY':#.startswith('EXPDTA    SOLUTION NMR'):
                        Exp_dict[pdbid_dict] = 'X-RAY'
                    # print (Exp_dict)  
                    if 'NMR' in Exp_dict.values() :
                        
                        if lines[0]=='MODEL' and lines[1]=='1':#("MODEL        1"):
                            NMR_printing = True
                        elif lines[0]=='MODEL' and lines[1]!='1':
                            NMR_printing = False
                        if NMR_printing:
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
                    if 'X-RAY' in Exp_dict.values():
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

# print (PDB_ids_list)
# print ("pdb_seq_dict", pdb_seq_dict)
# print ("pdb_seq_dict_numbering", pdb_seq_dict_numbering)
# print('<br/>')
# print('<br/>')
for keys,vals in pdb_seq_dict.items():
    aas_string=''.join(map(str,vals))
    pdb_seqfin.setdefault(keys,[]).append(aas_string)

# print (pdb_seqfin)

combine=dict(list(non_pdb_seq_dictfin.items()) + list(pdb_seqfin.items()))
# print (combine)

with open('trimmedfasta.fasta', 'w') as files:
    for seqids, seqn in combine.items():
        files.write( ">" + seqids)
        files.write("\n")
        files.write (seqn[0])
        files.write("\n")
        
#muscle_cline = MuscleCommandline(input="trimmedfasta.fasta", out=out_file)
Omega_out = subprocess.call([clustalo, '--infile', 'trimmedfasta.fasta','--outfile',  out_file])
seq1 = SeqIO.parse(out_file, 'fasta')
SeqIO.write(seq1, "file_tabDNA.fasta", "tab")
record_seq_dict = SeqIO.to_dict(SeqIO.parse(out_file, "fasta"))
# Sequence_aligned_df = pd.DataFrame(record_seq_dict)
# PDB_Sequence_aligned_df = Sequence_aligned_df[PDB_ids_list]


## Downloading NUCplot
aa = '''
(ala
(arg
(asn
(asp
(cys
(glu
(gln
(gly
(his
(ile
(leu
(lys
(met
(phe
(pro
(ser
(thr
(trp
(tyr
(val
'''.split()


aa = [each_string.title() for each_string in aa]

##NUCPLOT
DNA_binding_res_list2 = []
DNA_binding_res_list1 = OrderedDict()
with open("file_tabDNA.fasta",'r') as files:
    for line in files:
        pdbcode1=line.split()[0]
        # print (pdbcode1)
        pdbcode1_split=pdbcode1.split(':')
        if len(pdbcode1_split) >1:
            url = NUcplot_URL1 + pdbcode1_split[0] +NUcplot_URL2
            # print (url)

            data = requests.get(url, allow_redirects=True)
            nuclplotfilename = pdbcode1_split[0] + '_Nucplot.txt'
            open(nuclplotfilename, 'wb').write(data.content)



            with open(nuclplotfilename, 'r') as file1:
                for lines in file1:
                    if lines.startswith(tuple(aa)):
                        # print (lines)
                        lines_split_aminoacid= lines.split()[0].split('(', 1)[1].split('*')[0]
                        # print (lines_split_aminoacid)
                        if not lines_split_aminoacid in DNA_binding_res_list2:
                            DNA_binding_res_list2.append(lines_split_aminoacid.replace('))', ')'))
                            DNA_binding_res_list1.setdefault(pdbcode1, []).append(lines_split_aminoacid.replace('))', ')'))
                    aa1= [aa_string.replace('(','') for aa_string in aa]
                    t = filter(lambda x:x in lines, aa1)
                    for ts in t:
                        # print (lines)
                        if not lines.startswith(tuple(aa)):
                            if not lines.startswith('%'):
                                # print (lines)
                                lines_split2_aa = lines.split()[1].split('*')[0]
                                # print (lines_split2_aa)
                                if not lines_split2_aa in DNA_binding_res_list2:
                                    DNA_binding_res_list2.append(lines_split2_aa.replace('))', ')'))
                                    DNA_binding_res_list1.setdefault(pdbcode1, []).append(lines_split2_aa.replace('))', ')'))
            
                
# print (DNA_binding_res_list1)
# print (len(DNA_binding_res_list1))
# print('<br/>')
# print('<br/>')

DNA_binding_res_singlecode_dict={}
for PDBidKey, aavalue in DNA_binding_res_list1.items():
    aa_temp =[]
    # print (PDBidKey, aavalue)
    for aavalue1 in aavalue:
        # print (aavalue1)
        first_slit_chain = aavalue1.split('(')
        # print (first_slit_chain)
        aaminoacid = re.split(r'(\d+)', first_slit_chain[0])[0].upper()
        aaminoacidposition = re.split(r'(\d+)', first_slit_chain[0])[1]                     
        # print (aaminoacid, aaminoacidposition)
        if aaminoacid in aa_dict:
                
            aminoacid_position_value= aa_dict[aaminoacid] + '_' + aaminoacidposition
            # print (aminoacid_position_value)
            if not aminoacid_position_value in aa_temp:
                aa_temp.append(aminoacid_position_value)
                DNA_binding_res_singlecode_dict.setdefault(PDBidKey, []).append(aminoacid_position_value)
            
        # print (re.split('([^a-zA-Z0-9])', aavalue1))
# print (DNA_binding_res_singlecode_dict)

# print('<br/>')
# print('<br/>')
for ids,aa_numb in DNA_binding_res_singlecode_dict.items():
    for ids1, numb in pdb_seq_dict_numbering.items():
        if ids==ids1:
        #if ids==str(ids1).split(':')[0]:
            for j in aa_numb:
                numb_split=j.split('_')
                updated_indexed_bindingsite.setdefault(ids,[]).append(numb.index(numb_split[1]))
# print (updated_indexed_bindingsite)


####################
def Convert(string):
    list1=[]
    list1[:0]=string
    return list1
lines_DSSP_list=[]
lines_DSSP_dict ={}
always_print = False
with open("file_tabDNA.fasta",'r') as files:
    for line in files:
        # print (line)
        line1=line.split()[1:]
        SS_list=[]
        pdbids_chain=line.split()[:1]
        
        pdbids_chain_split=''.join(pdbids_chain).split(':')
        if len(pdbids_chain_split) > 1:
            
            check_ids=''.join(pdbids_chain).split(':')[0]#[1:]
            
            pdb_chain=''.join(pdbids_chain).split(':')[1]
            # print(pdb_chain)
            dssp_urlpart12_concat= dssp_urlpart1+check_ids+dssp_urlpart2
            dssp_allDBpage = requests.get(dssp_urlpart12_concat)
            dsspsource = BeautifulSoup(dssp_allDBpage.content, 'html.parser')
            for MRSallDB_link in dsspsource.find_all('a'):
                for dssplink in MRSallDB_link:
                    # print (dssplink)
                    if dssplink==(check_ids.upper()):
                        
                        dssp_href= str(MRSallDB_link).split(';')
                        
                        
                        if dssp_href[0]=='<a href="entry?db=dssp&amp':
                            
                            nr_number=dssp_href[1].split('=')[1].split('&')[0]
                            dssp_finalurl_concatenate= dssp_final_url1 + nr_number+ dssp_final_url2+check_ids.upper()
                            dsspwebpagelink=requests.get(dssp_finalurl_concatenate)
                            for dssppagelines in dsspwebpagelink.iter_lines():
                                lines_in_dssp = dssppagelines
                                # lines_in_dssp=str(dssppagelines.decode('utf8')).strip()
                                # print (lines_in_dssp)
                                if not lines_in_dssp.startswith('<') and len(lines_in_dssp)>3:
                                    linesdssp_part1=lines_in_dssp.split()
                                    if len(linesdssp_part1)>5 and linesdssp_part1[2]==pdb_chain:
                                        if str(linesdssp_part1[4]).startswith(('H','B','E','G','I','T', 'S')):
                                            # print(linesdssp_part1)
                                            if len(linesdssp_part1[4])==1:
                                                aa_SS=linesdssp_part1[4]
                                                SS_list.append(aa_SS)
                                            else:
                                                aa_NoSS='-'
                                                SS_list.append(aa_NoSS)
                                        else:
                                            aa_NoSS='-'
                                            SS_list.append(aa_NoSS)
                            # print ("SS_list",SS_list)
                            #Secondary structure adjustment according to its sequence alignment
                            if len(pdbids_chain_split) > 1:
                                # print (line1)
                                for seq_residues in line1:
                                    found2=[res_letter for res_letter in list(seq_residues)]
                                    # print ("found2", found2)
                                sscount=0
                                sscount1=0
                                num=0
                                SS_list_alignment_adjusted=[]
                                for  char,ssi in zip(found2,itertools.cycle(SS_list)): 
                                    
                                    if char!='-'  :
                                        # print ("SS_list[num]", SS_list[num], num)
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
        # print (line1)
        aa_text1=list(str(line1[0]))
        aa_text=aa_text1#[:-2]
        # print (check_ids, aa_text)
        aa_text_dict[check_ids]=aa_text
        for aa_t in range(0,len(''.join(aa_text)),60):
            AA_temp=[]
            AA_temp.append(''.join(aa_text)[aa_t:aa_t+60])
            # print (check_ids,'--->',AA_temp[0])
            if check_ids.startswith('>'):
                check_ids=''.join(check_ids.split('>'))
                # print (check_ids,'--->',AA_temp[0])
                Prot_AA.setdefault(check_ids, []).append(AA_temp[0])
            else:
                # print (check_ids,'--->',AA_temp[0])
                Prot_AA.setdefault(check_ids, []).append(AA_temp[0])
            # print ("Prot_AA",Prot_AA)

####updating index number for DNA binding site based on alignment for the given sequences
print('<br/>')
print('<br/>')
for  Hkey1,Hvalues1 in updated_indexed_bindingsite.items():  
    # print (Hkey1,Hvalues1)
    H_count =0
    H_count1=0
    # print ("Hkey1",Hkey1)
    Hkey1=str(Hkey1).split(':')[0]
    # print (aa_text_dict.keys())
    if Hkey1 in aa_text_dict.keys():
        for H_X, H_i in enumerate(aa_text_dict[str(Hkey1).split(':')[0]]):
                if H_i=='-':
                    H_count+=1
                else:
                     if H_count1 in Hvalues1:
                         Alignment_adjusted_indexed_Hbindingsite.setdefault(Hkey1, []).append(H_X)
                     H_count1+=1

# print (Alignment_adjusted_indexed_Hbindingsite)

prot_Idsname=list(Prot_AA.keys())
digits1=len(max(prot_Idsname))
f1 = '{0:>%d}: ' % (digits1)


print("<b>Color coded secondary structure elements </b> ")
print('<br/>')
structurecode="<span style='background-color:#8B008B'><font color='white'>B = beta-bridge residue</font></span>, <span style='background-color:#FFFF00'>E = extended strand (in beta ladder)</span>, <span style='background-color:#CD5C5C'><font color='white'>G = 3/10-helix</font></span>, <span style='background-color:#FF0000'><font color='white'>H= alpha-helix</font></span>,    <span style='background-color:#FA8072'>I = Pi helix</span>, <span style='background-color: #00FF00'> S = bend </span>, <span style='background-color: #008000'><font color='white'>T = H-bonded turn</font></span>"
print(structurecode)
print('<br/>')
print('<br/>')
print("<table style= 'border: 0'>")
##holding the keys of  Prot_AA dictionary
SeqID_ref= list(sorted(Prot_AA.keys()))[0]

for itemlength in range(len(Prot_AA[SeqID_ref])):
    for dictkeys in sorted(Prot_AA.keys()):
        gapcount=0
        resindex=0
        
        if dictkeys in DSSP_SS:
            print("<tr style= 'border: 0'>")
            # print("<td")
            # print("</td>")
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
            # print("<pre>")
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
            # print("</pre>")
            print("</tr>")
        print("<tr style= 'border: 0'>")
        print("<td style= 'border: 0';width=5%'>")
        dictkeys2=dictkeys+'\t'
        dictkeys1="<text-align='left'>%s" %dictkeys2
        print(dictkeys1)
        print("</td>")
        ##print secondary structure and sequences
        print("<td style= 'border: 0'; width=80%'>")
        # print("<pre>")
        for i in Prot_AA[dictkeys][itemlength]:
            
            if i=='-':
                res_size= "<span style='font-family:monospace;font-size:18px'>%s</span>" %i
                print(res_size)
                count+=1
                gapcount+=1
            elif dictkeys in Alignment_adjusted_indexed_Hbindingsite.keys() or dictkeys in Alignment_adjusted_indexed_NONHbindingsite.keys():
                if count in Alignment_adjusted_indexed_Hbindingsite[dictkeys]: #and not Alignment_adjusted_indexed_NONHbindingsite[dictkeys]:
                    
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
        # print("</pre>")
        print("</tr>")

    print("<tr style= 'border: 0; heigth: 60px '>")
    print("<td style= 'border: 0; height: 60px'>")
    print("</td>")
    print("</tr>")

        
                
print("</table>")

print('<br/>')

#making a list of binding residue indexes (according to aligned positions) from unique_mergeddict
H_NHnew_list=[]
for H_NHkey in Alignment_adjusted_indexed_Hbindingsite.keys():
    for index_items in Alignment_adjusted_indexed_Hbindingsite[H_NHkey]:
        if index_items not in H_NHnew_list:
            H_NHnew_list.append(index_items)
# print ("H_NHnew_list",H_NHnew_list)

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


# making H_NH_weblogo
# protmotif=[]
# H_NONH='HNONH'+'.svg'
# for webid in weblogo_gap_to_X.keys():
#     var=''.join(weblogo_gap_to_X[webid])
#     protmotif.append(Seq(var, IUPAC.extended_protein))
# seq = motifs.create(protmotif) 
# # print (seq)
# seq.weblogo(H_NONH,format='SVG',xaxis_label= 'Residues', show_errorbars= False, color_scheme= 'color_chemistry')
# H_NONHsrc="<center><embed src='%s#page=1&view=FitH ' /></center>"%H_NONH
# print(H_NONHsrc)
# print("<br/>")
# print("<br/>")

#print("Alignment of DNA interaction Residues")
print("<div align='center'><b>Alignment  of DNA binding residues</b></div>")
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
        # print("<pre>")
        for H_NH_res in weblogo_align[H_NH_dictkeys]:
            if H_NH_ident_count in H_NH_identical_index:
                H_NH_underlined="<SPAN STYLE='background-color:red; font-weight:bold; color:white'>%s</SPAN>" %H_NH_res
                print(H_NH_underlined)
                H_NH_ident_count+=1
            else:
                H_NH_no_conse="<b><SPAN STYLE='background-color:powderblue; font-weight:bold; color:black'>%s</SPAN>" %H_NH_res+'</b>'
                print(H_NH_no_conse)
                H_NH_ident_count+=1
    else:
        #boxing.write("<tr>")
        #boxing.write("<td>")
        print(H_NH_dictkeys)
        print("</td>")
        print("<td style= 'border: 0'>")
        # print("<pre>")
        for H_NH_res in weblogo_align[H_NH_dictkeys]:
            if H_NH_ident_count in H_NH_identical_index:
                H_NH_underlined1="<SPAN STYLE='background-color:red; font-weight:bold; color:white'>%s</SPAN>" %H_NH_res
                print(H_NH_underlined1)
                H_NH_ident_count+=1
            else:
                H_NH_no_conse1="<b><SPAN STYLE='background-color:powderblue; font-weight:bold; color:black'>%s</SPAN>" %H_NH_res+'</b>'
                print(H_NH_no_conse1)
                H_NH_ident_count+=1
        
    # print("</pre>")
    print("</td>")
    print("</tr>")
print("</table>") 

os.remove(out_file) 
os.remove('trimmedfasta.fasta') 
shutil.rmtree('obsolete') 
shutil.rmtree('PDB')



print "</body>"
print "</html>"
