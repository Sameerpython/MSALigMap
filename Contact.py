#!/usr/bin/env python3



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
#container1{
display:flex;
flex-direction:column;
align-item:flex-start;
border:1px solid black;
padding-top:20px;
-webkit-box-shadow: 5px 5px 15px 5px #000000; 
box-shadow: 5px 5px 15px 5px #000000;

}
.container4{
margin-bottom:10px;
padding_bottom:10px;
           }
#container3{
margin-top:70px;
}

""")
#########

print ("div.container { width:100%; border: 1px solid grey;}")
print ("ul{list-style-type: none;margin: 0;padding: 0; overflow: hidden;background-color: #333333;}")
print ("li{float:left;}")
print ("li a {display: block;color: white;text-align: center;padding: 16px;font-size:20px; text-decoration: none;}")
print ("li a:hover { background-color: #111111;}")
print (".footer { position: fixed; left: 0; bottom: 0; width: 100%; height:60px;  background-color: #808080; color: white; text-align: center; }")
print (" #container1{ width:900px; height=1500px; line-height:0;}")

#style for divindg into 2 columns

print ("* {box-sizing: border-box;}")
print (".column {float: left;width: 50%;padding: 10px;height: 300px;}")
print (".row:after {content: "";display: table;clear: both;}")


print ("</style>")

print ("<title>LiBiSCo</title>")
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

#Body size Determination


print ("<div align='center'>")
print ("<h1> Please send any comments or suggestions about the MSALigMap server to:</h1>")
print ("<div id='container1'>")



print("<div  >")

print ("<IMG SRC='HA.jpg' ALIGN='top' width=120px height=120px  style='border-radius:15px';>")
print("<div>")


print("<div class='container4'>")
print ("<p><b>Professor Henrik Aronsson</b></p>")
print  ("<p>Head of Department</p>")
print ("<p>Department of Biological & Environmental Sciences</p>")
print ("<p>Email: henrik.aronsson@bioenv.gu.se</p>")
print("<div>")



print("<div  id='container3'>")
print ("<IMG SRC='SH3.jpeg' ALIGN='top' width=120px height=120px style='border-radius:15px';>")
print("<div>")

print("<div class='container4'>")
print ("<p><b>Sameer Hassan</b></p>")
print  ("<p>Postdoc Researcher</p>")
print ("<p>Department of Biosciences and Nutrition</p>")
print ("<p>Email: sameer.hassan@ki.se</p>")
print ("</div>")

print("<div id='container3'>")

print ("<IMG SRC='MT.jpg' ALIGN='top' width=120px height=120px style='border-radius:15px';>")
print("<div>")
print("<div class='container4'>")
print ("<p><b>Mats Topel</b></p>")
print  ("<p>Researcher</p>")
print ("<p>Department of Marine Sciences</p>")
print ("<p>Email: mats.topel@marine.gu.se</p>")
print ("</div>")

print("<div id='container3'>")
print ("<IMG SRC='SH3.jpeg' ALIGN='top' width=120px height=120px  style='border-radius:15px';>")
print("<div>")
print("<div class='container4'>")
print ("<p><b>Haleemath Sameena Sameer</b></p>")
print  ("<p>Master Student, Department of Computer and Systems Sciences</p>")
print ("<p>Stockholm University</p>")
print ("<p>Email: sameenah.tm@gmail.com </p>")
print ("</div>")

print ("</div>")
print ("</div>")

print ("</body>")
print ("</html>")
