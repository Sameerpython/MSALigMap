#!/usr/bin/python

# Import modules for CGI handling 
import cgi, cgitb 

# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from fields
if form.getvalue('subject'):
   subject = form.getvalue('subject')
else:
   subject = "Not set"


print "Content-type:text/html\r\n\r\n"
print "<html>"
print "<head>"

print "<style>"
#header styling
print"""
.images{
    width: 100%;
    height: 200px;
    background-image: url('Pic.gif') ;
    background-size: cover;
    padding-left: 80px;
    
}

"""
#########

print "div.container { width:100%; border: 1px solid grey;}"
print "ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}"
print "li{float:left;}"
print "li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}"
print "li a:hover { background-color: #111111;}"
print ".footer { position: fixed; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }"
print " #container1{ width:1000px; height=1500px; line-height:1.6;}"

#style for divindg into 2 columns

print "* {box-sizing: border-box;}"
print ".column {float: left;width: 50%;padding: 10px;height: 300px;}"
print ".row:after {content: "";display: table;clear: both;}"


print "</style>"

print "<title>LiBiSCo</title>"
print "</head>"
#print "<h1>CoFact<style=color:blue;>Comp</style></h1>"
#print "<div align='center'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "</div>"

#print "<div id='container3' style='position:relative;'>"
#print "<img src='Title_image1.png' align='middle' width='1000' height='200'"
#print "<h1 style='position:absolute; top:100px; left:20px;'>LiBisCo </h1>  "
#print "<h1 style='position:absolute; top:100px; left:20px;'> Ligand Binding Site Comparison </h1>   " 
#print "</div>"

print "<div class='images'>"
#print "<div class='main'>"
#print "<style='font-size:120px'>LiBiSCo</font>","</br>"
#print "<font color='red'> Li </font>gand <font color='red'>Bi</font>nding <font color='red'>S</font>ite <font color='red'>Co</font>mparison"
print"  </div>"
#print "</div>"



print "<body>"
print "<ul>"
print "<li><a href='HomePage.py'>Home</a></li>"
print "<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>"
print "<li><a href='ProteinDNA.py'>Protein-DNA</a></li>"
print "<li><a href='ProteinPeptide.py'>Protein-Peptide</a></li>"
print "<li><a href='Contact.py'>Contact</a></li>"
print "</ul>"
#Body size Determination


print "<div align='center'>"
print "<h2> Comparing the Binding Residues for the Selected Protein Structures</h2>"
#print "</br>"
print "<div id='container1'>"


print "<h2> Protein - Ligand Binding Site Mapping</h1>"

print "<p> Upload protein sequences in fasta format:</p>"
print  "<form enctype='multipart/form-data' action='LigPage.py' method = 'post' >"
print  "<input type='file' id='myFile' name='filename'>"

print "<p> Enter PDB ids and Ligand ids (eg: 3WXB:A_NDP):</p>"
print  "<textarea rows='2' cols='10' name = 'textcontent' cols = '10' rows = '10'>"
print  "</textarea>"

print " <p><input type = 'submit' value = 'Upload' /></p>"

print  " </form>"
print "</div>"




print "</div>"

#print "<div class='footer'>"
#print "<p><img src='University.jpg' float='left'  width='40' height='40'<a href='https://bioenv.gu.se/english/staff?userId=xarohe'>Prof. Henrik Aronsson Group</a>, Department of Biological and Environmental Sciences, University of Gothenburg, Sweden. </p>"
#print "</div>"

print "</body>"
print "</html>"
