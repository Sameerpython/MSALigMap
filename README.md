MSALigMap (Multiple Sequence Alignment Ligand Mapping) is a tool for mapping active site amino acid residues that bind selected ligands or DNA on to target protein sequences of interest. To run ligand (small molecule) or DNA binding site analysis, the user needs to specify Lig or DNA respectively.
MSALigMap is implemented in Python (3) and is run from the command line. To be able to use the tool one will need:

1. Python 3 or later 
2. Installation of Clustal Omega and keep the ClustalO.exe in the same working directory

Example:
Two input files are supplied to MSALigMap.py: one containing multiple sequences of interest in fasta format and a second file in a tab-delimited format containing PDB id with chain information and ligand codes ( e.g. 3WXB:A	NDP).

Below is an example shown how to run MSALigMap:
To run Ligand or DNA analysis, the user needs to specify Lig or DNA as mentioned below.

a. Ligand binding site analysis:

Usage: MSALigMap <Lig> <Input.fasta> <Ligand.txt>
Input.fasta: The input.fasta file to be provided is a fasta file containing protein sequences that includes sequences of PDB structures and non-PDB sequences
Ligand.txt: The input file is a text file containing the PDB code with the chain id and ligand code (as provided in PDB) in tab-delimited format. 

Thus, the format of PDB ids, followed by chain id and ligand names in ligand.txt should be as follows:

3WXB:A	NDP
3O26:A	NDP

The PDB code is followed by ':' and thereafter comes the chain id. The ligand name (code as provided in PDB) is separated by a tab.

In the sequence file, the sequence name for the PDB structure should be the same as provided in the ligand.txt file (see example below).

>3WXB:A
MGSSHHHHHHSSGLVPRGSHMGELRVRSVLVTGANRGIGLGFVQHLLALSNPPEWVFATCRDPKGQRAQELQKLASKHPNLVIVPLEVTDPASIKAAAASVGERLKGSGLNLLINNAGIARANTIDNETLKDMSEVYTTNTIAPLLLSQAFLPMLKKAAQENPGSGLSCSKAAIINISSTAGSIQDLYLWQYGQALSYRCSKAALNMLTRCQSMGYREHGIFCVALHPGWVKTDMGGTLEDKSRVTVDESVGGMLKVLSNLSEKDSGAFLNWEGKVMAW


b. DNA binding site analysis:
  
Usage: MSALigMap <DNA> <input.fasta> 
Input.fasta: The input.fasta file to be provided is a fasta file containing protein sequences that includes sequences of PDB structures and non-PDB sequences. For the PDB sequence, the chain id should be mentioned following the PDB code (6KKS:A)
>6KKS:A
GNNEYKKGLWTVEEDKILMDYVKAHGKGHWNRIAKKTGLKRCGKSCRLRWMNYLSPNVKR
GNFTEQEEDLIIRLHKLLGNRWSLIAKRVPGRTDNQVKNYWNTHLSKKLGIKDQKTKQS
