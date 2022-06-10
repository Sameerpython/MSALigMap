#!/usr/bin/env python3.9

# Import modules for CGI handling 
import cgi, cgitb 
# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from fields
if form.getvalue('subject'):
   subject = form.getvalue('subject')
else:
   subject = "Not set"


print("Content-type:text/html\r\n\r\n")
print("<html>")
print("<head>")

print("<style>")
#header styling
"""print
.images{
    width: 100%;
    height: 200px;
    background-image: url('Pic.gif') ;
    background-size: cover;
    padding-left: 80px;
    
}

"""
#########

print("div.container { width:100%; border: 1px solid grey;}")
print("p.hello {font-size:10px; color:blue}")
print("ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}")
print("li{float:left;}")
print("li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}")
print("li a:hover { background-color: #111111;}")
print(".footer { position: fixed; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }")
print(" #container1{ width:1000px; height=1500px; line-height:1.6; font-size:22px;color:black}")

##CSS
#print

print("</style>")

print("<title>MSALigMap</title>")
print("</head>")
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
#print "<div align='center'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "</div>"

#print "<div id='container3' style='position:relative;'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "<h1 style='position:absolute; top:100px; left:20px;'>LiBisCo </h1>  "
#print "<h1 style='position:absolute; top:100px; left:20px;'> Ligand Binding Site Comparison </h1>   " 
#print "</div>"

print("<div class='images'>")
#print "<div class='main'>"
#print "<style='font-size:120px'>LiBiSCo</font>","</br>"
#print "<font color='red'> Li </font>gand <font color='red'>Bi</font>nding <font color='red'>S</font>ite <font color='red'>Co</font>mparison"
print("  </div>")
#print "</div>"



print("<body>")
print("<ul>")
print("<li><a href='HomePage.py'>Home</a></li>")
print("<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>")
print("<li><a href='ProteinDNA.py'>Protein-DNA</a></li>")
print("<li><a href='ProteinPeptide.py'>Protein-Peptide</a></li>")
print("<li><a href='Contact.py'>Contact</a></li>")
print("</ul>")
#Body size Determination


print("<div align='center'>")
print("<h1 font-size='30px'> Comparing the Binding Residues for the Selected Protein Sequences</h1>")
#print "</br>"
print("<div id='container1'>")
print("<p id='hello' align='justify' >Proteins binds to ligands (small molecules, DNA, RNA, peptide and proteins) to perform all kinds of important and essential cellular processes. They bind via a network of weak, noncovalent intermolecular interactions such as hydrogen bonding, hydrophobic and electrostatic interactions. Thus, binding of substrate is required for many proteins to function properly. Ligands are recognized in the binding pockets of the proteins through amino acids of specific properties. With increasing number of genomes having now sequenced, computational approach is the method of choice for characterizing proteins using homology-based information. Homology based approach is widely used for functional annotating sequences for identifying ligand, peptide and DNA binding sites etc from proteins with known function that are bound with ligand molecules. The MSALigMap is a tool that can be used for annotating protein sequences based on homology using protein structural information. User should submit protein sequences (homologous) and homologous protein structures bound to ligands of interest or DNA to map the binding sites for the uncharacterized protein sequences.</p>")



print("</div>")
print("</div>")
#print "<div class='footer'>"
#print "<p><img src='University.jpg' float='left'  width='40' height='40'<a href='https://bioenv.gu.se/english/staff?userId=xarohe'>Prof. Henrik Aronsson Group</a>, Department of Biological and Environmental Sciences, University of Gothenburg, Sweden. </p>"
#print "</div>"

print("</body>")
print("</html>")
