#!/usr/bin/env python3

# Import modules for CGI handling 
import cgi, cgitb 
# Create instance of FieldStorage
form = cgi.FieldStorage()

# Get data from fields
if form.getvalue('subject'):
   subject = form.getvalue('subject')
else:
   subject = "Not set"


print ("Content-type:text/html\r\n\r\n")
print ("<html>")
print ("<head>")

print ("<style>")
#header styling
print ("""
.images{
    width: 100%;
    height: 200px;
    background-image: url('Pic.gif') ;
    background-size: cover;
    padding-left: 80px;
    
}

""")
#########

print ("div.container { width:100%; border: 1px solid grey;}")
print ("ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}")
print ("li{float:left;}")
print ("li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}")
print ("li a:hover { background-color: #111111;}")
print (".footer { position: fixed; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }")
print (" #container1{ width:1000px; height=1500px; line-height:1.6;}")

#style for divindg into 2 columns

print ("* {box-sizing: border-box;}")
print (".column {float: left;width: 50%;padding: 10px;height: 300px;}")
print (".row:after {content: "";display: table;clear: both;}")


print ("</style>")

print ("<title>MSALigMap</title>")
print ("</head>")


print ("<div class='images'>")
print (" </div>")



print ("<body>")
print ("<ul>")
print ("<li><a href='HomePage.py'>Home</a></li>")
print ("<li><a href='ProteinLigand.py'>Protein-Ligand</a></li>")
print ("<li><a href='ProteinDNA.py'>Protein-DNA</a></li>")
print ("<li><a href='ProteinPeptide.py'>Protein-Peptide</a></li>")
print ("<li><a href='Contact.py'>Contact</a></li>")
print ("</ul>")


print ("<div align='center'>")
print ("<h2> Comparing the Binding Residues for the Selected Protein Sequences</h2>")
print ("<div id='container1'>")



print ("<h2> Protein - DNA binding Site Mapping</h1>")


print ("<p> Upload protein sequences in fasta format:</p>")
print  ("<form action='LigPageDNA.py' method = 'post' target = '_blank'>")
print  ("<input type='file' id='myFile' name='filename'>")

print (" <p><input type = 'submit' value = 'Submit' /></p>")

print  (" </form>")
print ("</div>")




print ("</div>")


print ("</body>")
print ("</html>")

