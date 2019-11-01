MSALigMap (Multiple Sequence Alignment Ligand Mapping) is a tool for mapping active site amino acid residues that bind selected ligands on to target protein sequences of interest.
MSALigMap is implemented in Python (2.7) and is run from the command line. To be able to use the tool one will need:

1. Python 2.7 or later 
2. Installation of Clustal Omega

Example:
Two input files are supplied to MSALigMap.py: one containing multiple sequences of interest in fasta format and a second file in a tab-delimited format containing PDB id with chain information and ligand codes ( e.g. 3WXB:A	NDP).

Below is an example shown how to run MSALigMap:

Usage: MSALigMap <input.fasta> <ligand.txt>
Input: The first input file to be provided is a fasta file containing protein sequences that includes sequences of PDB structures and non-PDB sequences
Ligand: The second input file is a text file containing the PDB code with the chain id and ligand code (as provided in PDB) in tab-delimited format. 

Thus, the format of PDB ids, followed by chain id and ligand names in ligand.txt should be as follows:

3WXB:A	NDP
3O26:A	NDP

The PDB code is followed by ':' and thereafter comes the chain id. The ligand name (code as provided in PDB) is separated by a tab.

In the sequence file, the sequence name for the PDB structure should be the same as provided in the ligand.txt file (see example below).

>3WXB:A
MGSSHHHHHHSSGLVPRGSHMGELRVRSVLVTGANRGIGLGFVQHLLALSNPPEWVFATCRDPKGQRAQELQKLASKHPNLVIVPLEVTDPASIKAAAASVGERLKGSGLNLLINNAGIARANTIDNETLKDMSEVYTTNTIAPLLLSQAFLPMLKKAAQENPGSGLSCSKAAIINISSTAGSIQDLYLWQYGQALSYRCSKAALNMLTRCQSMGYREHGIFCVALHPGWVKTDMGGTLEDKSRVTVDESVGGMLKVLSNLSEKDSGAFLNWEGKVMAW



